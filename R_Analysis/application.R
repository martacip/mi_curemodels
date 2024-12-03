library(dplyr)
library(gtsummary)
library(Hmisc)
library(plyr)
library(knitr)
library(curephEM)
library(ggplot2)
library(patchwork)
library(mice)
library(survival)
library(survminer)

source('functions.R')

load("case_study.RData")
#pnr: patient ID
#start: entry date
#dor: date of randomization
#end: end of follow-up 
#fut: time of follow-up
#stat: status
#trt: allocated treatment (0 conventional, 1 DI)
#age: age (in days)
#sex: sex (0 male, 1 female)
#hist: histological response (0 poor, 1 good)

data = Dhist_reduced %>% select(pnr, start, dor, end, fut, stat, trt, age, sex, hist)
data$fut = (as.numeric(data$fut))/365.25 #follow-up time in years
data$age = (as.numeric(data$age))/365.25 #age in years

data$trt = as.factor(data$trt)
levels(data$trt) = c("Reg-C", "Reg-DI")
label(data$trt) = "Allocated treatment"

data$sex = as.factor(data$sex)
levels(data$sex) = c("Male", "Female")
data$sex = ordered(data$sex, levels = c("Female", "Male"))
label(data$sex) = "Sex"

label(data$age) = "Age"

data$hist = as.factor(data$hist)
levels(data$hist) = c("Poor", "Good")
label(data$hist) = "Histological response"





###Complete-case analysis 
data_cc = data[complete.cases(data),]
model = cureph(Surv.cure(fut,stat) ~ hist + sex, formula2 = ~ hist + trt, data = data_cc)
summary(model)


#################################################################################


###Multiple imputation
data_miss = data
colnames(data_miss)[5] = "Y"
colnames(data_miss)[6] = "status"
colnames(data_miss)[7] = "Z"
data_miss$Z = ifelse(data_miss$Z=="Reg-C",0,1)
colnames(data_miss)[9] = "X"
data_miss$X = ifelse(data_miss$X=="Female",0,1)
colnames(data_miss)[10] = "W"
data_miss$W = ifelse(data_miss$W=="Poor",0,ifelse(data_miss$W=="Good",1,NA))

### STEP 0: Initialization
alpha_X0 = model[["coefficients"]][["logistic"]][["(Intercept)"]]
alpha_W = model[["coefficients"]][["logistic"]][["histGood"]]
alpha_X = model[["coefficients"]][["logistic"]][["sex.L"]]
beta_W = model[["coefficients"]][["cox"]][["histGood"]]
beta_Z = model[["coefficients"]][["cox"]][["trtReg-DI"]]

fit_cox = coxph(Surv(Y, status) ~ W + Z, data = data_miss, method = "breslow")
risk = basehaz(fit_cox, centered = FALSE)
colnames(risk)[2] = "Y"
data_miss = left_join(data_miss, risk)
colnames(data_miss)[11] = "H0" 


#randomly filling missing values of W
missing_indices = which(is.na(data_miss$W))
for (i in missing_indices) {
  data_miss$W[i] = runif(1)
}


result_list_exact = list()
result_list_approx = list()


###Imputation
K = 10
M = 10
W_distribution = "Bernoulli"

for (k in 1:K) {
  for (m in 1:M) {
    
    ### STEP 1: Estimate the baseline cumulative hazard rate H0(Y)
    
    # uncure probability
    data_miss$X0 = 1
    data_miss$pi = as.vector(p_func(c(alpha_X0, alpha_W, alpha_X), 
                                    cbind(data_miss$X0, data_miss$W, data_miss$X)))
    P = data_miss$pi
    
    # survival probability  
    n = 429
    data_miss$G_expect = rbinom(n, size = 1, prob = P)
    uncured = data_miss[data_miss$G_expect==1,]
    uncured$surv = surv_uncured_func(uncured, beta_W, beta_Z)
    data_miss$surv = 1
    data_miss[data_miss$G_expect==1,]$surv = uncured$surv
    
    # estimate of the posterior expectation Q
    plateau = max(data_miss$Y[data_miss$status==1])
    data_miss$G = ifelse(data_miss$status==1, 1, ifelse(data_miss$Y > plateau, 0, NA))
    data_miss$Q = Q_func(data_miss, plateau)
    
    # estimate of the baseline cumulative hazard
    data_miss$H0 = H0_func(data_miss, c(beta_W, beta_Z))
    
    
    
    
    ### STEP 2: Draw estimates for model parameters alpha0, alpha, and beta
    
    # Fit logistic regression model
    cure_param = cure_func(data_miss)
    
    # Fit Cox model
    cox_param = cox_uncured_func(data_miss)
    
    # Draw values of alpha and beta from Normal distribution
    parameters = lin_pred_func(data_miss, cure_param, cox_param)
    alpha_X0 = parameters$alpha_X0 
    alpha_W = parameters$alpha_W
    alpha_X = parameters$alpha_X
    beta_W = parameters$beta_W
    beta_Z = parameters$beta_Z
    data_miss$lp_cure = parameters$lp_cure
    data_miss$lp_surv = parameters$lp_surv
    
    
    
    ### STEP 3: Impute cure status G
    data_miss$G = impute_G_func(data_miss)
    
    
    
    ### STEP 4: Impute missing covariate W
    first_iter = (m==1)
    
    # Exact conditional distribution
    exact = imputation_exact_func(data_miss,missing_indices,alpha_X0,alpha_W,alpha_X,beta_W,beta_Z,
                                  first_iter)
    
    
    # Approximated conditional distribution
    approx = imputation_approx_func(data_miss, missing_indices, k)
    
    
  }
  
  result_list_exact[[k]] = exact
  result_list_approx[[k]] = approx
  
}


#exact conditional distribution
real_parameters = c(1,1,1,1,1)
est_exact = mi_smcure_func(result_list_exact)
rubin_exact = rubin_rules_func(est_exact[[1]], est_exact[[2]])
names(rubin_exact[[1]]) = c("alpha_X0", "alpha_W", "alpha_X", "beta_W", "beta_Z")
names(rubin_exact[[2]]) = c("alpha_X0", "alpha_W", "alpha_X", "beta_W", "beta_Z")
coef_exact_df = data.frame(exact_alpha0 = rubin_exact[[1]][["alpha_X0"]],
                           exact_alphaW = rubin_exact[[1]][["alpha_W"]],
                           exact_alphaX = rubin_exact[[1]][["alpha_X"]],
                           exact_betaW = rubin_exact[[1]][["beta_W"]],
                           exact_betaZ = rubin_exact[[1]][["beta_Z"]])
sd_exact_df = data.frame(exact_alpha0 = rubin_exact[[2]][["alpha_X0"]],
                         exact_alphaW = rubin_exact[[2]][["alpha_W"]],
                         exact_alphaX = rubin_exact[[2]][["alpha_X"]],
                         exact_betaW = rubin_exact[[2]][["beta_W"]],
                         exact_betaZ = rubin_exact[[2]][["beta_Z"]])


#approximated conditional distribution
est_approx = mi_smcure_func(result_list_approx)
rubin_approx = rubin_rules_func(est_approx[[1]], est_approx[[2]])
names(rubin_approx[[1]]) = c("alpha_X0", "alpha_W", "alpha_X", "beta_W", "beta_Z")
names(rubin_approx[[2]]) = c("alpha_X0", "alpha_W", "alpha_X", "beta_W", "beta_Z")
coef_approx_df = data.frame(approx_alpha0 = rubin_approx[[1]][["alpha_X0"]],
                            approx_alphaW = rubin_approx[[1]][["alpha_W"]],
                            approx_alphaX = rubin_approx[[1]][["alpha_X"]],
                            approx_betaW = rubin_approx[[1]][["beta_W"]],
                            approx_betaZ = rubin_approx[[1]][["beta_Z"]])
sd_approx_df = data.frame(approx_alpha0 = rubin_approx[[2]][["alpha_X0"]],
                          approx_alphaW = rubin_approx[[2]][["alpha_W"]],
                          approx_alphaX = rubin_approx[[2]][["alpha_X"]],
                          approx_betaW = rubin_approx[[2]][["beta_W"]],
                          approx_betaZ = rubin_approx[[2]][["beta_Z"]])







###################################################
### SATURATED MODEL

#Complete-case
data_cc$interaction = ifelse(data_cc$hist=="Good" & data_cc$trt=="Reg-DI", 1, 0)
model_bis = cureph(Surv.cure(fut,stat) ~ sex + hist + trt + interaction,
                   formula2 = ~ sex + hist + trt + interaction, 
                   data = data_cc)
summary(model_bis)


#Multiple imputation
data_miss = data
colnames(data_miss)[5] = "Y"
colnames(data_miss)[6] = "status"
colnames(data_miss)[7] = "Z"
data_miss$Z = ifelse(data_miss$Z=="Reg-C",0,1)
colnames(data_miss)[9] = "X"
data_miss$X = ifelse(data_miss$X=="Female",0,1)
colnames(data_miss)[10] = "W"
data_miss$W = ifelse(data_miss$W=="Poor",0,ifelse(data_miss$W=="Good",1,NA))
data_miss$R = data_miss$W*data_miss$Z


### STEP 0: Initialization
alpha_X0 = model_bis[["coefficients"]][["logistic"]][["(Intercept)"]]
alpha_W = model_bis[["coefficients"]][["logistic"]][["histGood"]]
alpha_X = model_bis[["coefficients"]][["logistic"]][["sex.L"]]
alpha_Z = model_bis[["coefficients"]][["logistic"]][["trtReg-DI"]]
alpha_R = model_bis[["coefficients"]][["logistic"]][["interaction"]]
beta_W = model_bis[["coefficients"]][["cox"]][["histGood"]]
beta_X = model_bis[["coefficients"]][["cox"]][["sex.L"]]
beta_Z = model_bis[["coefficients"]][["cox"]][["trtReg-DI"]]
beta_R = model_bis[["coefficients"]][["cox"]][["interaction"]]

fit_cox = coxph(Surv(Y, status) ~ W + X + Z + R, data = data_miss, method = "breslow")
risk = basehaz(fit_cox, centered = FALSE)
colnames(risk)[2] = "Y"
data_miss = left_join(data_miss, risk)
colnames(data_miss)[12] = "H0" 


#randomly filling missing values of W
missing_indices = which(is.na(data_miss$W))
for (i in missing_indices) {
  data_miss$W[i] = runif(1)
}
data_miss$R = data_miss$W*data_miss$Z


result_list_exact = list()
result_list_approx = list()

source('functions_5.R')

###Imputation
K = 10
M = 10
W_distribution = "Bernoulli"

for (k in 1:K) {
  for (m in 1:M) {
    
    ### STEP 1: Estimate the baseline cumulative hazard rate H0(Y)
    
    # uncure probability
    data_miss$X0 = 1
    data_miss$pi = as.vector(p_func(c(alpha_X0, alpha_W, alpha_X, alpha_Z, alpha_R), 
                                    cbind(data_miss$X0, data_miss$W, data_miss$X, data_miss$Z, data_miss$R)))
    P = data_miss$pi
    
    # survival probability  
    n = 429
    data_miss$G_expect = rbinom(n, size = 1, prob = P)
    uncured = data_miss[data_miss$G_expect==1,]
    uncured$surv = surv_uncured_func(uncured, beta_W, beta_X, beta_Z, beta_R)
    data_miss$surv = 1
    data_miss[data_miss$G_expect==1,]$surv = uncured$surv
    
    # estimate of the posterior expectation Q
    plateau = max(data_miss$Y[data_miss$status==1])
    data_miss$G = ifelse(data_miss$status==1, 1, ifelse(data_miss$Y > plateau, 0, NA))
    data_miss$Q = Q_func(data_miss, plateau)
    
    # estimate of the baseline cumulative hazard
    data_miss$H0 = H0_func(data_miss, c(beta_W, beta_X, beta_Z, beta_R))
    
    
    
    
    ### STEP 2: Draw estimates for model parameters alpha0, alpha, and beta
    
    # Fit logistic regression model
    cure_param = cure_func(data_miss)
    
    # Fit Cox model
    cox_param = cox_uncured_func(data_miss)
    
    # Draw values of alpha and beta from Normal distribution
    parameters = lin_pred_func(data_miss, cure_param, cox_param)
    alpha_X0 = parameters$alpha_X0 
    alpha_W = parameters$alpha_W
    alpha_X = parameters$alpha_X
    alpha_Z = parameters$alpha_Z
    alpha_R = parameters$alpha_R
    beta_W = parameters$beta_W
    beta_X = parameters$beta_X
    beta_Z = parameters$beta_Z
    beta_R = parameters$beta_R
    data_miss$lp_cure = parameters$lp_cure
    data_miss$lp_surv = parameters$lp_surv
    
    
    
    ### STEP 3: Impute cure status G
    data_miss$G = impute_G_func(data_miss)
    
    
    
    ### STEP 4: Impute missing covariate W
    first_iter = (m==1)
    
    # Exact conditional distribution
    exact = imputation_exact_func(data_miss,missing_indices,alpha_X0,alpha_W,alpha_X,alpha_Z,alpha_R,
                                  beta_W,beta_X,beta_Z,beta_R,
                                  first_iter)
    
    
    # Approximated conditional distribution
    approx = imputation_approx_func(data_miss, missing_indices, k)
    
    
  }
  
  result_list_exact[[k]] = exact
  result_list_approx[[k]] = approx
  
}


#exact conditional distribution
real_parameters = c(1,1,1,1,1,1,1,1,1)
est_exact = mi_smcure_func(result_list_exact)
rubin_exact = rubin_rules_func(est_exact[[1]], est_exact[[2]])
names(rubin_exact[[1]]) = c("alpha_X0", "alpha_W", "alpha_X", "alpha_Z", "alpha_R", "beta_W", "beta_X", "beta_Z", "beta_R")
names(rubin_exact[[2]]) = c("alpha_X0", "alpha_W", "alpha_X", "alpha_Z", "alpha_R", "beta_W", "beta_X", "beta_Z", "beta_R")
coef_exact_df = data.frame(exact_alpha0 = rubin_exact[[1]][["alpha_X0"]],
                           exact_alphaW = rubin_exact[[1]][["alpha_W"]],
                           exact_alphaX = rubin_exact[[1]][["alpha_X"]],
                           exact_alphaZ = rubin_exact[[1]][["alpha_Z"]],
                           exact_alphaR = rubin_exact[[1]][["alpha_R"]],
                           exact_betaW = rubin_exact[[1]][["beta_W"]],
                           exact_betaX = rubin_exact[[1]][["beta_X"]],
                           exact_betaZ = rubin_exact[[1]][["beta_Z"]],
                           exact_betaR = rubin_exact[[1]][["beta_R"]])
sd_exact_df = data.frame(exact_alpha0 = rubin_exact[[2]][["alpha_X0"]],
                         exact_alphaW = rubin_exact[[2]][["alpha_W"]],
                         exact_alphaX = rubin_exact[[2]][["alpha_X"]],
                         exact_alphaZ = rubin_exact[[2]][["alpha_Z"]],
                         exact_alphaR = rubin_exact[[2]][["alpha_R"]],
                         exact_betaW = rubin_exact[[2]][["beta_W"]],
                         exact_betaX = rubin_exact[[2]][["beta_X"]],
                         exact_betaZ = rubin_exact[[2]][["beta_Z"]],
                         exact_betaR = rubin_exact[[2]][["beta_R"]])


#approximated conditional distribution
est_approx = mi_smcure_func(result_list_approx)
rubin_approx = rubin_rules_func(est_approx[[1]], est_approx[[2]])
names(rubin_approx[[1]]) = c("alpha_X0", "alpha_W", "alpha_X", "alpha_Z", "alpha_R", "beta_W", "beta_X", "beta_Z", "beta_R")
names(rubin_approx[[2]]) = c("alpha_X0", "alpha_W", "alpha_X", "alpha_Z", "alpha_R", "beta_W", "beta_X", "beta_Z", "beta_R")
coef_approx_df = data.frame(approx_alpha0 = rubin_approx[[1]][["alpha_X0"]],
                            approx_alphaW = rubin_approx[[1]][["alpha_W"]],
                            approx_alphaX = rubin_approx[[1]][["alpha_X"]],
                            approx_alphaZ = rubin_approx[[1]][["alpha_Z"]],
                            approx_alphaR = rubin_approx[[1]][["alpha_R"]],
                            approx_betaW = rubin_approx[[1]][["beta_W"]],
                            approx_betaX = rubin_approx[[1]][["beta_X"]],
                            approx_betaZ = rubin_approx[[1]][["beta_Z"]],
                            approx_betaR = rubin_approx[[1]][["beta_R"]])
sd_approx_df = data.frame(approx_alpha0 = rubin_approx[[2]][["alpha_X0"]],
                          approx_alphaW = rubin_approx[[2]][["alpha_W"]],
                          approx_alphaX = rubin_approx[[2]][["alpha_X"]],
                          approx_alphaZ = rubin_approx[[2]][["alpha_Z"]],
                          approx_alphaR = rubin_approx[[2]][["alpha_R"]],
                          approx_betaW = rubin_approx[[2]][["beta_W"]],
                          approx_betaX = rubin_approx[[2]][["beta_X"]],
                          approx_betaZ = rubin_approx[[2]][["beta_Z"]],
                          approx_betaR = rubin_approx[[2]][["beta_R"]])
