args = commandArgs(trailingOnly = TRUE)

library(smcure)
library(curephEM)
library(mice)
library(survival)
library(survminer)
library(dplyr)
library(extRemes)
require(doParallel)
require(foreach)


####################
# initial settings #
####################

scenario = "3"

#import functions
source('functions.R')

#import parameters
source('parameters_Marta.R')




comp_coef = function(n, M, K, W_distribution, seed=NULL){
  if (is.null(seed))  seed = sample(1:10000, 1)
  set.seed(seed)
  
  
  #############
  # variables #
  #############
  
  # XZW = mvtnorm::rmvnorm(n, c(0.5, 0.5, 0.5), matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1), nrow=3, ncol=3))
  # X = rep(0,n)
  # X[XZW[,1]>0.5] = 1
  # Z = rep(0,n)
  # Z[XZW[,2]>0.5] = 1
  # 
  # if (W_distribution=="Normal"){
  #   W = XZW[,3]}
  # else if (W_distribution=="Bernoulli"){
  #   W = rep(0,n)
  #   W[XZW[,3]>0.5] = 1
  # }
  
  X = rbinom(n, 1, 0.5)
  Z = rbinom(n, 1, 0.5)

  if (W_distribution=="Normal"){
    W = rnorm(n, 0.5)}
  else if (W_distribution=="Bernoulli"){
    W = rbinom(n, 1, 0.5)
  }
  
  
  
  
  #############
  # incidence #
  #############
  
  X0 = rep(1, n) #intercept
  lp_incidence = cbind(X0, W, X) #incidence design matrix
  
  P = as.vector(p_func(alpha_values, lp_incidence)) #uncure probability
  
  #generation of cure status 
  G_expect = rbinom(n, size = 1, prob = P) #G=1 if uncured 
  
  
  
  ###########
  # latency #
  ###########
  
  #generate survival times using a Weibull model
  lp_latency = cbind(W, Z)
  
  Y = rep(0, n)
  
  u = runif(n)
  time = rep(0,n)
  time = as.numeric((-log(u) / (lambda*exp(lp_latency%*%beta_values)))^(1/rho))
  # time1 = rep(0, n)
  # time1 = as.numeric((-log(u) / (lambda*exp(-rho*lp_latency%*%gamma_values)))^(1/rho))
  # time2 = rep(0, n)
  # time2 = exp(mi + lp_latency%*%gamma_values + sigma*revd(n))
  
  time[time > tau_0] = tau_0 #truncation of survival times 
  # time1[time1 > tau_0] = tau_0
  # time2[time2 > tau_0] = tau_0
  # summary(time)
  # summary(time1)
  # summary(time2)
  
  time[which(G_expect==0)] = 20000
  Y = time   
  
  
  #censoring   
  C = as.numeric(rexp(n, cens_rate))
  C[C > tau_1] = tau_1  #truncate censoring times
  
  
  #follow-up time
  TgreatC = which(time > C)
  Y[TgreatC] = C[TgreatC]
  status = as.numeric(time <= C)
  
  
  # #checking survival time generation
  # summary(coxph(Surv(Y, status) ~ W + Z))
  # 
  # #checking cured fraction generation
  # summary(glm(status ~ W + X))
  
  
  
  #observed uncured indicator
  dataKM = as.data.frame(cbind(Y, status))
  fitKM = survfit(Surv(Y, status) ~ 1, data = dataKM)
  plateau = fitKM$time[tail(which((diff(fitKM$surv) < 0) == T), 1) + 1] #to find the elbow of the survival curve
  G = ifelse(status==1, 1, ifelse(Y > plateau, 0, NA))
  
  
  
  #summary statistics on cure and censoring rates
  nobs_plateau = sum(Y > plateau)
  infomat = matrix(0, nrow = 1, ncol = 6)
  colnames(infomat) = c("n", "Cure rate", "Cens.rate", "Cens.level","% obs in plateau", "% cured in plateau")
  infomat[1, 1] = n                                     # Sample size
  infomat[1, 2] = round(length(which(G_expect == 0)) / n * 100, 3)       # Cure rate
  infomat[1, 3] = cens_rate
  infomat[1, 4] = round(sum(1 - status) / n * 100, 3)    # Censoring level
  infomat[1, 5] = round(nobs_plateau / n * 100, 3)      # % of obs. in plateau
  infomat[1, 6] = round(sum(G_expect[which(Y > plateau)] == 0) / nobs_plateau * 100, 3)
  
  
  #create data frame for the analysis with no missing values
  data_full = data.frame(Y, status, W, X0, X, Z, G)
  
  
  ##missing mechanisms
  W_mcar_15 = rep(0, n)
  W_mcar_15[sample(1:n, 0.15*n)] = NA
  
  # MCAR with 30% missing values
  W_mcar_30 = rep(0, n)
  W_mcar_30[sample(1:n, 0.3*n)] = NA
  
  # MAR with logit form and 30% missing values 
  #Calculate response propensity:
  logit_prob = 0.3 - 0.5*X - 0.5*Z  #response model
  rp = exp(logit_prob) / (exp(logit_prob) + 1) # Suppress values between 0 and 1 via inverse-logit
  
  # rp can be seen as probability to respond in y. 
  # See literature about reponse propensity for more details
  
  # Create missings based on rp
  W_mar = rbinom(n, 1, rp)
  
  data = data_full
  data_mar = ampute(data[,-7], prop = 0.3, patterns = c(1,1,0,1,1,1), mech = "MAR")
  data_mar_bis = data_mar$amp
  
  if (miss_mech=="MCAR_15"){  
    W[is.na(W_mcar_15)] = NA} 
  else if (miss_mech=="MCAR_30"){  
    W[is.na(W_mcar_30)] = NA}
  else if (miss_mech=="MAR"){  
    W = data_mar_bis$W}
  
  
  
  ##create a dataframe for the MI and CCA
  data_miss = data.frame(Y, status, W, X0, X, Z, G)
  
  data_cca = data_miss[!is.na(data_miss$W),]

  
  
  
  ##############
  # imputation #
  ##############
  
  ### STEP 0: Initialization
  
  #estimate of alpha and beta on the full dataset
  # fit_cure_sm = smcure(Surv(Y,status) ~ W + Z, cureform = ~ W + X,
  #                   data = data_full[,-dim(data_full)[2]], model="ph")
  fit_cure = cureph(Surv.cure(Y,status) ~ W + X, formula2 = ~ W + Z,
                    data = data_full[,-dim(data_full)[2]])
  alpha_X0 = fit_cure[["coefficients"]][["logistic"]][["(Intercept)"]]
  alpha_W = fit_cure[["coefficients"]][["logistic"]][["W"]]
  alpha_X = fit_cure[["coefficients"]][["logistic"]][["X"]]
  beta_W = fit_cure[["coefficients"]][["cox"]][["W"]]
  beta_Z = fit_cure[["coefficients"]][["cox"]][["Z"]]
  
  #estimate of H0 on the full data 
  fit_cox = coxph(Surv(Y, status) ~ W + Z, data = data_full, method = "breslow")
  risk = basehaz(fit_cox, centered = FALSE)
  colnames(risk)[2] = "Y"
  data_miss = left_join(data_miss, risk)
  colnames(data_miss)[8] = "H0" 
  
  
  #randomly filling missing values of W
  missing_indices = which(is.na(data_miss$W))
  for (i in missing_indices) {
    data_miss$W[i] = runif(1)
  }
  
  
  
  
  result_list_exact = list()
  result_list_approx = list()
  
  
  for (k in 1:K) {
    for (m in 1:M) {
      
      ### STEP 1: Estimate the baseline cumulative hazard rate H0(Y)
      
      # uncure probability
      data_miss$pi = as.vector(p_func(c(alpha_X0, alpha_W, alpha_X), cbind(data_miss$X0, data_miss$W, 
                                                                           data_miss$X)))
      
      # survival probability  
      uncured = data_miss[G_expect==1,]
      uncured$surv = surv_uncured_func(uncured, beta_W, beta_Z)
      data_miss$surv = 1
      data_miss[G_expect==1,]$surv = uncured$surv
      
      # estimate of the posterior expectation Q
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
  
  
  
  
  
  ##############
  # cure model #
  ##############
  
  #no missing values
  # cm_all = smcure(Surv(Y,status) ~ W + Z, cureform = ~ W + X, data = data_full, model="ph")
  # coef_all_df = data.frame(full_alpha0 = cm_all[["b"]][["(Intercept)"]],
  #                          full_alphaW = cm_all[["b"]][["Z[, -1]W"]],
  #                          full_alphaX = cm_all[["b"]][["Z[, -1]X"]],
  #                          full_betaW = cm_all[["beta"]][["X[, -1]W"]],
  #                          full_betaZ = cm_all[["beta"]][["X[, -1]Z"]])
  # sd_all_df = data.frame(full_alpha0 = cm_all[["b_sd"]][[1]],
  #                        full_alphaW = cm_all[["b_sd"]][[2]],
  #                        full_alphaX = cm_all[["b_sd"]][[3]],
  #                        full_betaW = cm_all[["beta_sd"]][[1]],
  #                        full_betaZ = cm_all[["beta_sd"]][[2]])
  cm_all = cureph(Surv.cure(Y,status) ~ W + X, formula2 = ~ W + Z,
                    data = data_full[,-dim(data_full)[2]])
  coef_all_df = data.frame(full_alpha0 = cm_all[["coefficients"]][["logistic"]][["(Intercept)"]],
                           full_alphaW = cm_all[["coefficients"]][["logistic"]][["W"]],
                           full_alphaX = cm_all[["coefficients"]][["logistic"]][["X"]],
                           full_betaW = cm_all[["coefficients"]][["cox"]][["W"]],
                           full_betaZ = cm_all[["coefficients"]][["cox"]][["Z"]])
  reslog = summary(cm_all)[["coefficients"]][["logistic"]]
  rescox = summary(cm_all)[["coefficients"]][["cox"]]
  sd_all_df = data.frame(full_alpha0 = reslog[1,3],
                         full_alphaW = reslog[2,3],
                         full_alphaX = reslog[3,3],
                         full_betaW = rescox[1,3],
                         full_betaZ = rescox[2,3])

    
  
  #complete-case analysis
  # cm_cc = smcure(Surv(Y,status) ~ W + Z, cureform = ~ W + X, data = data_cca, model="ph")
  # coef_cca_df = data.frame(cc_alpha0 = cm_cc[["b"]][["(Intercept)"]],
  #                          cc_alphaW = cm_cc[["b"]][["Z[, -1]W"]],
  #                          cc_alphaX = cm_cc[["b"]][["Z[, -1]X"]],
  #                          cc_betaW = cm_cc[["beta"]][["X[, -1]W"]],
  #                          cc_betaZ = cm_cc[["beta"]][["X[, -1]Z"]])
  # sd_cca_df = data.frame(cc_alpha0 = cm_cc[["b_sd"]][[1]],
  #                        cc_alphaW = cm_cc[["b_sd"]][[2]],
  #                        cc_alphaX = cm_cc[["b_sd"]][[3]],
  #                        cc_betaW = cm_cc[["beta_sd"]][[1]],
  #                        cc_betaZ = cm_cc[["beta_sd"]][[2]])
  cm_cc = cureph(Surv.cure(Y,status) ~ W + X, formula2 = ~ W + Z,
                  data = data_cca)
  coef_cca_df = data.frame(cc_alpha0 = cm_cc[["coefficients"]][["logistic"]][["(Intercept)"]],
                           cc_alphaW = cm_cc[["coefficients"]][["logistic"]][["W"]],
                           cc_alphaX = cm_cc[["coefficients"]][["logistic"]][["X"]],
                           cc_betaW = cm_cc[["coefficients"]][["cox"]][["W"]],
                           cc_betaZ = cm_cc[["coefficients"]][["cox"]][["Z"]])
  reslog = summary(cm_cc)[["coefficients"]][["logistic"]]
  rescox = summary(cm_cc)[["coefficients"]][["cox"]]
  sd_cca_df = data.frame(cc_alpha0 = reslog[1,3],
                         cc_alphaW = reslog[2,3],
                         cc_alphaX = reslog[3,3],
                         cc_betaW = rescox[1,3],
                         cc_betaZ = rescox[2,3])
  
  
  #exact conditional distribution
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
  
  
  
  
  #create a matrix with all the regression coefficients
  return(list(coef = cbind(seed, coef_all_df, coef_cca_df, coef_exact_df, coef_approx_df),
              sd = cbind(seed, sd_all_df, sd_cca_df, sd_exact_df, sd_approx_df),
              infomat = cbind(seed, infomat)))
  
}



#I = 1; n = 500; M = 1; K = 1; W_distribution = "Normal"
I = as.numeric(args[1]) #number of data sets
n = as.numeric(args[2]) #number of observations
M = as.numeric(args[3]) #number of iterations 
K = as.numeric(args[4]) #number of imputed data sets
W_distribution = args[5] #type of distribution (Bernoulli or Normal)


# ncores = detectCores() - 2
cl = makeCluster(120)
registerDoParallel(cl)

coef_mat = foreach(i=1:I, .export = c(".GlobalEnv"), .combine = "rbind", .packages = c("mice","smcure","survival","survminer","dplyr","extRemes","curephEM")) %dopar%{
  comp_coef(n, M, K, W_distribution)
}

stopCluster(cl)


save.image(file = paste0("simulation_", W_distribution, "_", scenario, ".RData"))
