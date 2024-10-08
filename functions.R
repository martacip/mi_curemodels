### Functions to be used in the imputation algorithm


#uncure probability 
p_func = function(alpha,x){
  xalpha = as.numeric(x %*% alpha)
  r = (1 + exp(-xalpha)) ^ (-1) #logistic
  return(r)}



#estimate the survival of the uncured
surv_uncured_func = function(df, beta_W, beta_Z){
  df$Surv = exp(-(df$H0*exp(beta_W*df$W + beta_Z*df$Z)))
  return(df$Surv)
}



#estimate Q 
Q_func = function(df, plateau){
  df$Q[df$status==1] = 1
  df$Q[df$status==0] = (df$pi[df$status==0]*df$surv[df$status==0])/
    (1-df$pi[df$status==0]+df$pi[df$status==0]*df$surv[df$status==0])
  df$Q[df$Y>plateau]=0
  return(df$Q)
}



#estimate H0
H0_func = function(df, beta){
  df$H0Y = 0
  #loop through each row of df
  for(i in 1:length(df$Y)){
    H0t = 0
    event_times = df$Y[which(df$status==1)]
    for(tj in event_times[event_times<=df$Y[i]]){
      # number events at time tj
      Dj = sum(event_times==tj)
      # i in Rj (risk set at time tj (df$Y>=tj))
      Rj=df[df$Y>=tj,]
      r = Rj$Q * exp((cbind(Rj$W,Rj$Z))%*%beta)
      denominator = sum(r)
      H0t = H0t+Dj/denominator
    }
    df$H0Y[i] = H0t
  }
  return(df$H0Y)
}



#estimate parameters of cure (pi=P(Gi=1|Xi))
cure_func = function(df){
  pi = glm(G ~ W + X, family = binomial(link='logit'), data = df)
  cure_param = list(coefs = coef(pi), vcov_matrix = vcov(pi))
  return(cure_param)
}



#estimate parameters of survival of the uncured
cox_uncured_func = function(df){
  df_uncured = df[df$G==1,]
  surv_object = Surv(time = df_uncured$Y, event = df_uncured$status)
  coxph_model = coxph(surv_object ~ W + Z, data = df_uncured)
  cox_param = list(coefs = coef(coxph_model), vcov_matrix = vcov(coxph_model))
  return(cox_param)
}



#draw alphas and and add columns with linear predictors
lin_pred_func = function(df, cure_param, cox_param){
  alpha_star = as.data.frame(mvtnorm::rmvnorm(1, cure_param$coefs, cure_param$vcov_matrix))
  beta_star = as.data.frame(mvtnorm::rmvnorm(1, cox_param$coefs, cox_param$vcov_matrix))
  
  colnames(alpha_star) = paste("alpha", colnames(alpha_star), sep = "_")
  colnames(alpha_star)[colnames(alpha_star) == "alpha_(Intercept)"] = "alpha_X0"
  
  colnames(beta_star) = paste("beta", colnames(beta_star), sep="_")
  
  alpha_X0 = alpha_star$alpha_X0
  alpha_W = alpha_star$alpha_W
  alpha_X = alpha_star$alpha_X
  
  beta_W = beta_star$beta_W
  beta_Z = beta_star$beta_Z
  
  df$lp_cure = alpha_X0 + df$W*alpha_W + df$X*alpha_X
  df$lp_surv = df$W*beta_W + df$Z*beta_Z
  
  return(list(alpha_X0 = alpha_X0, 
              alpha_W = alpha_W,
              alpha_X = alpha_X,
              beta_W = beta_W,
              beta_Z = beta_Z,
              lp_cure = df$lp_cure,
              lp_surv = df$lp_surv))
}



#function to impute G
impute_G_func = function(df){
  #logistic function
  r = -df$H0*exp(df$lp_surv) + df$lp_cure
  r = 1/(1+exp(-r))
  #0 or 1 result
  G_imputed = as.numeric(c(runif(n)<=r))
  G_miss_index = is.na(df$G)
  df$G[G_miss_index] = G_imputed[G_miss_index]
  return(df$G)
}



#function to impute W using exact conditional distribution
imputation_exact_func = function(df,missing_indices,alpha_X0,alpha_W,alpha_X,beta_W,beta_Z,first_iter){
  
  if (first_iter) {
    #regression uses only observed values in the first iteration
    df_complete = df[-missing_indices,]
  } else {
    df_complete = df
  }
  
  
  #regression on W using X and Z
  if (W_distribution=="Bernoulli"){
    reg = glm(W ~ X + Z, family = binomial(link='logit'), data = df_complete) 
  }
  else if (W_distribution=="Normal"){
    reg = lm(W ~ X + Z, data = df_complete)
  }
  
  coefs = coef(reg)
  vcov_matrix = vcov(reg)
  
  
  #draw thetas from multivariate normal distribution
  theta_star = as.data.frame(mvtnorm::rmvnorm(1, coefs, vcov_matrix))
    
  colnames(theta_star) = paste("theta", colnames(theta_star), sep = "_")
  colnames(theta_star)[colnames(theta_star) == "theta_(Intercept)"] = "theta_X0"
    
  theta_X0 = theta_star$theta_X0
  theta_X = theta_star$theta_X
  theta_Z = theta_star$theta_Z
    
  if (W_distribution=="Bernoulli"){
      #logistic function
      r = (theta_X0 + df$X*theta_X + df$Z*theta_Z) + #mu_i
        (df$G*df$status*beta_W) - #first term
        (df$G*df$H0*(exp(beta_W)-1)*exp(df$Z*beta_Z)) + #second term
        (alpha_W*df$G) + #third term
        log(1 + exp(alpha_X0 + df$X*alpha_X)) - #fourth term
        log(1 + exp(alpha_X0 + df$W*alpha_W + df$X*alpha_X))
      R = 1/(1+exp(-r))
      
      #0 or 1 result
      W_imputed = as.numeric(c(runif(n)<=R))
    }
    
    else if (W_distribution=="Normal"){
      #kernel of the distribution
      mu = theta_X0 + df$X*theta_X + df$Z*theta_Z
      residuals = residuals(reg)
      sigma = var(residuals)
        
      target = function(W, X, Z, Y, status, G, alpha_X0, alpha_W, alpha_X, beta_W, beta_Z, H0, mu, sigma){
        kernel = ((exp(alpha_X0+X*alpha_X+W*alpha_W) / (1 + exp(alpha_X0+X*alpha_X+W*alpha_W)))^(G)) *
        ((exp((Z*beta_Z+W*beta_W)*status) * exp(-H0 * exp(Z*beta_Z+W*beta_W)))^(G)) *
        ((1 / (1 + exp(alpha_X0 + X*alpha_X + W*alpha_W)))^(1 - G)) *
        (exp((-(W - mu)^2) / (2*sigma^2)))
        return(kernel)
      }
      
      metropolis_step = function(W, X, Z, Y, status, G, alpha_X0, alpha_W, alpha_X, beta_W, beta_Z, H0, mu, sigma, sd){ 
        #the sd of the proposal distribution is a tuning parameter
        proposed_W = rnorm(1, mean = W, sd = sd)
        accept_prob = min(1, target(proposed_W, X, Z, Y, status, G, alpha_X0, alpha_W, alpha_X, beta_W, beta_Z, H0, mu, sigma) / target(W, X, Z, Y, status, G, alpha_X0, alpha_W, alpha_X, beta_W, beta_Z, H0, mu, sigma))
        u = runif(1)
        if (u <= accept_prob){
          value = proposed_W
          accepted = TRUE
        } else {
          value = W
          accepted = FALSE
        }
        out = data.frame(value = value, accepted = accepted)
        return(out)
      }
      
      
      metropolis_sampler = function(initial_value, X, Z, Y, status, G, alpha_X0, alpha_W, alpha_X, beta_W, beta_Z, H0, mu, sigma, n=500, sd = 1, burnin = 0, lag = 1){
        
        results = list()
        current_state = initial_value
        for (i in 1:burnin) {
          out = metropolis_step(current_state, X, Z, Y, status, G, alpha_X0, alpha_W, alpha_X, beta_W, beta_Z, H0, mu, sigma, sd)
          current_state = out$value
        }
        for (i in 1:n) {
          for (j in 1:lag) {
            out = metropolis_step(current_state, X, Z, Y, status, G, alpha_X0, alpha_W, alpha_X, beta_W, beta_Z, H0, mu, sigma, sd)
            current_state = out$value
          }
          results[[i]] = out
        }
        results = do.call(rbind, results)
        return(results)
      }
        
      out = metropolis_sampler(initial_value = df$W, df$X, df$Z, df$Y, df$status, df$G, alpha_X0, alpha_W, alpha_X, beta_W, beta_Z, df$H0, mu, sigma, burnin = 1000, lag = 100)
       
      W_imputed = out$value 
      }
    
    #put imputed W into original df and create list of df
    df$W[missing_indices] = W_imputed[missing_indices]
    
    return(df)
  }






#function to impute W using mice (approximated imputation for the CPH cure model)
imputation_approx_func = function(df, missing_indices, k){
  df_imp = data.frame(df)
  
  #select variables for imputation
  df_imp$G_status = df_imp$G*df_imp$status
  df_imp$G_H0 = df_imp$G*df_imp$H0
  df_imp$G_H0_X_Z = df_imp$G*df_imp$H0*df_imp$X*df_imp$Z
  df_imp = df_imp[ , c("W", "X", "Z", "G", "G_status", "G_H0", "G_H0_X_Z")]
  

  if (W_distribution=="Bernoulli"){
    #W as a factor for imputation
    df_imp$W[missing_indices] = NA
    df_imp$W = as.factor(df_imp$W)}
  else if (W_distribution=="Normal"){
    df_imp$W[missing_indices] = NA}
  
  #impute using MICE
  if (W_distribution=="Bernoulli"){
    imp = mice(df_imp, m = K, print=FALSE, method = "logreg", maxit = 1)}
  else if (W_distribution=="Normal"){
    imp = mice(df_imp, m = K, print=FALSE, method = "pmm", maxit = 1)}
  
  #put imputed W into original df and create list of df
  df_imputed = complete(imp, k)
  if (W_distribution=="Bernoulli"){
    df_imputed$W = as.numeric(levels(df_imputed$W))[df_imputed$W]}
    df$W[missing_indices] = df_imputed$W[missing_indices]
    
  return(df)
}




#compute the coefficients and standard errors of the cure ph model
mi_smcure_func = function(data_list){
  coef_mi_matrix = matrix(nrow = 0, ncol = length(real_parameters))
  se_mi_matrix = matrix(nrow = 0, ncol = length(real_parameters))
  
  for(j in 1:length(data_list)){
    df = data_list[[j]]
    
    #cure model for multiple imputation case
    log = capture.output({cm_mi = cureph(Surv.cure(Y,status) ~ W + X, formula2 = ~ W + Z,
                                 data = df)
    #standard error
    se_cm_mi = c(summary(cm_mi)[["coefficients"]][["logistic"]][,3],
                 summary(cm_mi)[["coefficients"]][["cox"]][,3])})
    
    coef_mi_matrix = rbind(coef_mi_matrix,c(summary(cm_mi)[["coefficients"]][["logistic"]][,1],
                       summary(cm_mi)[["coefficients"]][["cox"]][,1]))
    se_mi_matrix = rbind(se_mi_matrix,se_cm_mi)
  }
  return(list(coef_mi_matrix, se_mi_matrix))
}




#Rubin's Rules for estimators and standard errors
rubin_rules_func = function(coef_mi_matrix, se_mi_matrix){
  #coefficients
  U =  colMeans(coef_mi_matrix)
  
  #variance / standard error
  V = colMeans(se_mi_matrix^2)
  
  B_aux = sweep(coef_mi_matrix, 2, U, '-')^2
  B = colSums(B_aux)/(K-1)
  
  S = V+(1+1/K)*B
  sd = sqrt(S)
  
  return(list(U, sd))
}
