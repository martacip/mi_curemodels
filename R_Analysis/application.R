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



### Baseline characteristics
summary_1 = data %>% select(trt, sex, age, hist) %>%
  tbl_summary(
    type = list(all_continuous() ~ "continuous"),
    statistic = all_continuous() ~ "{median} ({min}, {max})",
    missing = "ifany"
  ) %>%
  bold_labels()%>%
  italicize_levels() %>%
  add_stat_label(
    label = list(all_categorical() ~ "n (%)", all_continuous() ~ "median (range)"))

summary_1 %>% as_flex_table()


dati = data
jpeg("pfs.jpg ", width = 9,  height = 6,  units = 'in',  res = 600)
ggsurvplot(
  fit = survfit(Surv(fut, stat) ~ 1, data = dati), 
  xlab = "Years since end of treatment", 
  ylab = "PFS probability",
  palette = "#DE1A1A",
  conf.int.fill = "#DE1A1A",
  xlim = c(0, 10),
  break.time.by = 2,
  risk.table = T,
  risk.table.height = 0.3,
  surv.scale="percent",
  conf.int.style="ribbon",conf.int = T,
  # legend = "none",
  title = "Osteosarcoma progression-free survival",
  tables.theme = theme_survminer(font.main = 0.01,font.x =11, font.y =11, font.tickslab = 11))
dev.off()


###Complete-case analysis 
data_cc = data[complete.cases(data),]
model = cureph(Surv.cure(fut,stat) ~ hist + sex, formula2 = ~ hist + trt, data = data_cc)
summary(model)

coef_data = data.frame(
  Variable = c("Intercept", "Hist Resp", "Sex"),
  OR = c(2.562, 0.293, 1.581),
  Lower = c(1.6450, 0.1653, 1.0938),
  Upper = c(3.9900, 0.5194, 2.3843)
)

coef_data$Variable = factor(coef_data$Variable, levels = rev(coef_data$Variable))

colors = c("Intercept" = "#b8b8b8", "Hist Resp" = "#1a80bb", 
           "Sex" = "#ea801c")

include_in_legend = c("Hist Resp", "Sex")

forest_plot1 = ggplot(coef_data, aes(x = OR, y = Variable)) +
  geom_point(aes(color = Variable), size = 5) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = colors, 
                     breaks = include_in_legend,
                     labels = include_in_legend) +
  labs(
    title = "Incidence",
    x = "Odds Ratio (OR)",
    y = NULL,
    caption = "CI: 95%"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 20, face = "bold"),  # Title size
    axis.title.x = element_text(size = 16),  # X-axis label size
    axis.text.x = element_text(size = 14),   # X-axis tick labels size
    axis.text.y = element_text(size = 14),   # Y-axis tick labels size
    legend.text = element_text(size = 14),   # Legend text size
    legend.title = element_text(size = 16),  # Legend title size
    plot.caption = element_text(size = 12)   # Caption size
  ) +
  xlim(0, 4)

forest_plot1 = forest_plot1 + 
  scale_y_discrete(labels = c(
    "Intercept" = "Intercept",
    "Hist Resp" = "Good vs Poor",
    "Sex" = "Male vs Female"))


print(forest_plot1)

coef_data = data.frame(
  Variable = c("Hist Resp", "Treatment"),
  HR = c(0.8678, 0.8151),
  Lower = c(0.5566, 0.5548),
  Upper = c(1.353, 1.197)
)

coef_data$Variable = factor(coef_data$Variable, levels = rev(coef_data$Variable))

colors = c("Hist Resp" = "#298c8c", "Treatment" = "#a00000")

forest_plot2 = ggplot(coef_data, aes(x = HR, y = Variable)) +
  geom_point(aes(color = Variable), size = 5) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = colors) +
  labs(
    title = "Latency",
    x = "Hazard Ratio (HR)",
    y = NULL,
    caption = "CI: 95%"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 20, face = "bold"),  # Title size
    axis.title.x = element_text(size = 16),  # X-axis label size
    axis.text.x = element_text(size = 14),   # X-axis tick labels size
    axis.text.y = element_text(size = 14),   # Y-axis tick labels size
    legend.text = element_text(size = 14),   # Legend text size
    legend.title = element_text(size = 16),  # Legend title size
    plot.caption = element_text(size = 12)   # Caption size
  ) +
  guides(color = guide_legend(reverse = TRUE)) +
  xlim(0.5, 1.5)


print(forest_plot2)

combined_plot = forest_plot1 + forest_plot2 + plot_layout(ncol = 2)

ggsave("plot.png", combined_plot, width = 11, height = 5, dpi = 300)


# pen = cox_cure_net(~ trt + sex + hist,
#                    ~ trt + sex + hist,
#                    time = fut, event = stat, data = data_cc,
#                    surv_alpha = 1, cure_alpha = 1)
# lapply(coef(pen), function(x) x[x != 0])

# ### BRIER SCORE
# n = 379
# data_cc = data_cc[sample(1:379, 379, replace = FALSE),]
# 
# X = as.matrix(data_cc[, c(9,10)])
# Z = as.matrix(data_cc[, 7])
# 
# Y = data_cc$fut  # follow-up times
# Status = data_cc$stat       # censoring indicator
# 
# T.max = as.numeric(max(Y[Status == 1]))  # cure threshold
# 
# dataset = data.frame(Y, Status, Z, X)
# 
# B = rep(-1, 379)
# B[Status == 1] = 1  
# B[Status == 0 & Y > T.max] = 0 
# 
# weights = rep(-1, 379)
# weights[B == -1] = 0
# 
# G = survfit(Surv(time = dataset$Y, event = 1 - dataset$Status) ~ 1, data = dataset)
# 
# G_surv_fun = function(t) {
#   r = t
#   for (j in 1:length(t)) {
#     ind = 1
#     while (G$time[ind] <= t[j] & ind <= length(G$time)) {
#       ind = ind + 1
#     }
#     if (ind == 1) {
#       r[j] = 1
#     }
#     if (ind == length(G$time)) {
#       r[j] = G$surv[ind]
#     }
#     if (ind > 1 & ind < length(G$time)) {
#       r[j] = G$surv[ind - 1]
#     }
#   }
#   return(r)
# }
# G_est = G_surv_fun(pmin(dataset$Y, rep(T.max, 379)))
# 
# weights[B != -1] = 1 / (n * G_est[B != -1])
# 
# K = 10
# folds_1 = cut(seq(1, n), breaks = K, labels = FALSE)  # 10-fold cross-validation
# 
# Brier_run = 1:10*NA
# 
# for (i in 1:10) {
#   folds = sample(folds_1,379,replace = FALSE)
#   Brier_val = 1:K
#   k = 0
#   
#   while (k < K) {
#     k = k + 1
#     train_data = which(folds != k)
#     
#     # Train set
#     X = as.matrix(data_cc[train_data, c(9, 10)])
#     Z = as.matrix(data_cc[train_data, 7])
#     
#     Y = data_cc$fut[train_data]  # follow-up times
#     Status = data_cc$stat[train_data]        # censoring indicator
#     
#     dataset = data.frame(Y, Status, Z, X)
#     dataset$Z = as.numeric(as.factor(dataset$Z))-1
#     dataset$sex = as.numeric(as.factor(dataset$sex))-1
#     dataset$hist = as.numeric(as.factor(dataset$hist))-1
#     
#     # Weights computation
#     beta.init = coxph(Surv(Y,Status) ~ hist + sex, subset = Status!=0,
#                       data = dataset, method = "breslow")$coef 
#     gamma.init = glm(Status ~ hist + Z, family = binomial(link=logit), data=dataset)$coefficients
#     
#     # Logistic/Cox mixture cure model
#     LCMM = smcure(Surv(Y,Status) ~ hist + sex, cureform = ~ hist + Z, model = "ph", 
#                   data = dataset, Var = FALSE) 
#     
#     # Test set
#     test_data = which(folds == k)
#     
#     X = as.matrix(data_cc[test_data, c(9, 10)])
#     Z = as.matrix(data_cc[test_data, 7])
#     
#     Y = data_cc$fut[test_data]  # follow-up times
#     Status = data_cc$stat[test_data]        # censoring indicator
#     
#     dataset = data.frame(Y, Status, Z, X)
#     dataset$Z = as.numeric(as.factor(dataset$Z))-1
#     dataset$sex = as.numeric(as.factor(dataset$sex))-1
#     dataset$hist = as.numeric(as.factor(dataset$hist))-1
#     
#     s = sort(LCMM$s, decreasing = TRUE)
#     Time_train = sort(LCMM$Time)
#     
#     surv_fun = function(t) {
#       r = t
#       for (j in 1:length(t)) {
#         ind = 1
#         while (Time_train[ind] <= t[j] & ind <= length(Time_train)) {
#           ind = ind + 1
#         }
#         if (ind == 1) {
#           r[j] = 1
#         }
#         if (ind == length(Time_train)) {
#           r[j] = s[ind]
#         }
#         if (ind > 1 & ind < length(Time_train)) {
#           r[j] = s[ind - 1]
#         }
#       }
#       return(r)
#     }
#     S_test = surv_fun(dataset$Y)
#     
#     Floor = sqrt(.Machine$double.eps)
#     
#     # Logistic/Cox mixture cure model prediction
#     ScoreLC = LCMM$b[1] + dataset$hist * LCMM$b[2] + dataset$Z * LCMM$b[3] 
#     dataset$predLC <- (exp(ScoreLC) / (1 + exp(ScoreLC)))  # Probability of being uncured    
#     
#     Brier_val[k] <- sum(weights[folds == k] * ((B[folds == k] - dataset$predLC)^2))
#   }   
#   
#   Brier_run[i] <- sum(Brier_val)
# }
# 
# 
# mean(Brier_run)
# sd(Brier_run)
# 
# 
# 
# 
# 
# 
# model_bis = cureph(Surv.cure(fut,stat) ~ hist + sex + trt + age, 
#                    formula2 = ~ hist + trt + sex + age, data = data_cc)
# summary(model_bis)
# 
# 
# data_cc = data_cc[sample(1:379, 379, replace = FALSE),]
# 
# X = as.matrix(data_cc[, c(9,10,7,8)])
# 
# Y = data_cc$fut  # follow-up times
# Status = data_cc$stat       # censoring indicator
# 
# T.max = as.numeric(max(Y[Status == 1]))  # cure threshold
# 
# dataset = data.frame(Y, Status, X)
# 
# B = rep(-1, 379)
# B[Status == 1] = 1  
# B[Status == 0 & Y > T.max] = 0 
# 
# weights = rep(-1, 379)
# weights[B == -1] = 0
# 
# G = survfit(Surv(time = dataset$Y, event = 1 - dataset$Status) ~ 1, data = dataset)
# 
# G_surv_fun = function(t) {
#   r = t
#   for (j in 1:length(t)) {
#     ind = 1
#     while (G$time[ind] <= t[j] & ind <= length(G$time)) {
#       ind = ind + 1
#     }
#     if (ind == 1) {
#       r[j] = 1
#     }
#     if (ind == length(G$time)) {
#       r[j] = G$surv[ind]
#     }
#     if (ind > 1 & ind < length(G$time)) {
#       r[j] = G$surv[ind - 1]
#     }
#   }
#   return(r)
# }
# G_est = G_surv_fun(pmin(dataset$Y, rep(T.max, 379)))
# 
# weights[B != -1] = 1 / (n * G_est[B != -1])
# 
# K = 10
# folds_1 = cut(seq(1, n), breaks = K, labels = FALSE)  # 10-fold cross-validation
# 
# Brier_run = 1:50*NA
# 
# for (i in 1:50) {
#   folds = sample(folds_1,379,replace = FALSE)
#   Brier_val = 1:K
#   k = 0
#   
#   while (k < K) {
#     k = k + 1
#     train_data = which(folds != k)
#     
#     # Train set
#     X = as.matrix(data_cc[train_data, c(9, 10, 7, 8)])
#     
#     Y = data_cc$fut[train_data]  # follow-up times
#     Status = data_cc$stat[train_data]        # censoring indicator
#     
#     dataset = data.frame(Y, Status, X)
#     dataset$trt = as.numeric(as.factor(dataset$trt))-1
#     dataset$sex = as.numeric(as.factor(dataset$sex))-1
#     dataset$hist = as.numeric(as.factor(dataset$hist))-1
#     dataset$age = as.numeric(dataset$age)
#     
#     # # Weights computation
#     # beta.init = coxph(Surv(Y,Status) ~ hist + sex + trt + age, subset = Status!=0,
#     #                   data = dataset, method = "breslow")$coef 
#     # gamma.init = glm(Status ~ hist + trt + sex + age, family = binomial(link=logit), data=dataset)$coefficients
#     
#     # Logistic/Cox mixture cure model
#     LCMM = smcure(Surv(Y,Status) ~ hist + sex + trt + age, cureform = ~ hist + trt + sex + age, model = "ph", 
#                   data = dataset, Var = FALSE) 
#     
#     # Test set
#     test_data = which(folds == k)
#     
#     X = as.matrix(data_cc[test_data, c(9, 10, 7, 8)])
#     
#     Y = data_cc$fut[test_data]  # follow-up times
#     Status = data_cc$stat[test_data]        # censoring indicator
#     
#     dataset = data.frame(Y, Status, X)
#     dataset$trt = as.numeric(as.factor(dataset$trt))-1
#     dataset$sex = as.numeric(as.factor(dataset$sex))-1
#     dataset$hist = as.numeric(as.factor(dataset$hist))-1
#     dataset$age = as.numeric(dataset$age)
#     
#     s = sort(LCMM$s, decreasing = TRUE)
#     Time_train = sort(LCMM$Time)
#     
#     surv_fun = function(t) {
#       r = t
#       for (j in 1:length(t)) {
#         ind = 1
#         while (Time_train[ind] <= t[j] & ind <= length(Time_train)) {
#           ind = ind + 1
#         }
#         if (ind == 1) {
#           r[j] = 1
#         }
#         if (ind == length(Time_train)) {
#           r[j] = s[ind]
#         }
#         if (ind > 1 & ind < length(Time_train)) {
#           r[j] = s[ind - 1]
#         }
#       }
#       return(r)
#     }
#     S_test = surv_fun(dataset$Y)
#     
#     Floor = sqrt(.Machine$double.eps)
#     
#     # Logistic/Cox mixture cure model prediction
#     ScoreLC = LCMM$b[1] + dataset$hist * LCMM$b[2] + dataset$trt * LCMM$b[3] + dataset$sex * LCMM$b[4] +
#       dataset$age * LCMM$b[5]
#     dataset$predLC <- (exp(ScoreLC) / (1 + exp(ScoreLC)))  # Probability of being uncured    
#     
#     Brier_val[k] <- sum(weights[folds == k] * ((B[folds == k] - dataset$predLC)^2))
#   }   
#   
#   Brier_run[i] <- sum(Brier_val)
# }
# 
# 
# mean(Brier_run)
# sd(Brier_run)





#################################################################################

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



coef_data = data.frame(
  Variable = c("Intercept", "Hist Resp", "Sex"),
  OR = c(1.853401, 0.2929949, 1.910775),
  Lower = c(1.143909, 0.1652835, 1.13517),
  Upper = c(3.002944, 0.5193865, 2.216311)
)

coef_data$Variable = factor(coef_data$Variable, levels = rev(coef_data$Variable))

colors = c("Intercept" = "#b8b8b8", "Hist Resp" = "#1a80bb", 
           "Sex" = "#ea801c")

include_in_legend = c("Hist Resp", "Sex")

forest_plot1 = ggplot(coef_data, aes(x = OR, y = Variable)) +
  geom_point(aes(color = Variable), size = 5) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = colors, 
                     breaks = include_in_legend,
                     labels = include_in_legend) +
  labs(
    title = "Incidence",
    x = "Odds Ratio (OR)",
    y = NULL,
    caption = "CI: 95%"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 20, face = "bold"),  # Title size
    axis.title.x = element_text(size = 16),  # X-axis label size
    axis.text.x = element_text(size = 14),   # X-axis tick labels size
    axis.text.y = element_text(size = 14),   # Y-axis tick labels size
    legend.text = element_text(size = 14),   # Legend text size
    legend.title = element_text(size = 16),  # Legend title size
    plot.caption = element_text(size = 12)   # Caption size
  ) +
  xlim(0, 4)

forest_plot1 = forest_plot1 + 
  scale_y_discrete(labels = c(
    "Intercept" = "Intercept",
    "Hist Resp" = "Good vs Poor",
    "Sex" = "Male vs Female"))

print(forest_plot1)

coef_data = data.frame(
  Variable = c("Hist Resp", "Treatment"),
  HR = c(0.8677649, 0.8150901),
  Lower = c(0.5566117, 0.5548125),
  Upper = c(1.352857, 1.137471)
)

coef_data$Variable = factor(coef_data$Variable, levels = rev(coef_data$Variable))

colors = c("Hist Resp" = "#298c8c", "Treatment" = "#a00000")

forest_plot2 = ggplot(coef_data, aes(x = HR, y = Variable)) +
  geom_point(aes(color = Variable), size = 5) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = colors) +
  labs(
    title = "Latency",
    x = "Hazard Ratio (HR)",
    y = NULL,
    caption = "CI: 95%"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 20, face = "bold"),  # Title size
    axis.title.x = element_text(size = 16),  # X-axis label size
    axis.text.x = element_text(size = 14),   # X-axis tick labels size
    axis.text.y = element_text(size = 14),   # Y-axis tick labels size
    legend.text = element_text(size = 14),   # Legend text size
    legend.title = element_text(size = 16),  # Legend title size
    plot.caption = element_text(size = 12)   # Caption size
  ) +
  guides(color = guide_legend(reverse = TRUE))+
  xlim(0.5, 1.5)

print(forest_plot2)

combined_plot2 = forest_plot1 + forest_plot2 + plot_layout(ncol = 2)

ggsave("plot2.png", combined_plot2, width = 11, height = 5, dpi = 300)











coef_data = data.frame(
  Variable = c("Intercept", "Hist Resp", "Sex"),
  OR = c(1.842783, 0.3377349, 1.78545),
  Lower = c(1.187122, 0.1990126, 1.10987),
  Upper = c(2.86057, 0.573154, 2.372254)
)

coef_data$Variable = factor(coef_data$Variable, levels = rev(coef_data$Variable))

colors = c("Intercept" = "#b8b8b8", "Hist Resp" = "#1a80bb", 
           "Sex" = "#ea801c")

include_in_legend = c("Hist Resp", "Sex")

forest_plot1 = ggplot(coef_data, aes(x = OR, y = Variable)) +
  geom_point(aes(color = Variable), size = 5) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = colors, 
                     breaks = include_in_legend,
                     labels = include_in_legend) +
  labs(
    title = "Incidence",
    x = "Odds Ratio (OR)",
    y = NULL,
    caption = "CI: 95%"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 20, face = "bold"),  # Title size
    axis.title.x = element_text(size = 16),  # X-axis label size
    axis.text.x = element_text(size = 14),   # X-axis tick labels size
    axis.text.y = element_text(size = 14),   # Y-axis tick labels size
    legend.text = element_text(size = 14),   # Legend text size
    legend.title = element_text(size = 16),  # Legend title size
    plot.caption = element_text(size = 12)   # Caption size
  ) +
  xlim(0, 4)

forest_plot1 = forest_plot1 + 
  scale_y_discrete(labels = c(
    "Intercept" = "Intercept",
    "Hist Resp" = "Good vs Poor",
    "Sex" = "Male vs Female"))

print(forest_plot1)

coef_data = data.frame(
  Variable = c("Hist Resp", "Treatment"),
  HR = c(0.9274517, 0.8004891),
  Lower = c(0.6185817, 0.5730739),
  Upper = c(1.390547, 1.11815)
)

coef_data$Variable = factor(coef_data$Variable, levels = rev(coef_data$Variable))

colors = c("Hist Resp" = "#298c8c", "Treatment" = "#a00000")

forest_plot2 = ggplot(coef_data, aes(x = HR, y = Variable)) +
  geom_point(aes(color = Variable), size = 5) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = colors) +
  labs(
    title = "Latency",
    x = "Hazard Ratio (HR)",
    y = NULL,
    caption = "CI: 95%"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 20, face = "bold"),  # Title size
    axis.title.x = element_text(size = 16),  # X-axis label size
    axis.text.x = element_text(size = 14),   # X-axis tick labels size
    axis.text.y = element_text(size = 14),   # Y-axis tick labels size
    legend.text = element_text(size = 14),   # Legend text size
    legend.title = element_text(size = 16),  # Legend title size
    plot.caption = element_text(size = 12)   # Caption size
  ) +
  guides(color = guide_legend(reverse = TRUE)) +
  xlim(0.5, 1.5)

print(forest_plot2)

combined_plot3 = forest_plot1 + forest_plot2 + plot_layout(ncol = 2)

ggsave("plot3.png", combined_plot3, width = 11, height = 5, dpi = 300)




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
