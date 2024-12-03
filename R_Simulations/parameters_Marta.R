#parameters

if (scenario=='1'){  
  #incidence
  alpha_values= c(1, -1, 0.5)
  #latency
  gamma_values = c(0.138, 0)
  rho = 1.45
  lambda = 0.25 
  beta_values = -gamma_values*rho
  mi = -log(lambda)/rho
  sigma = 1/rho
  #others
  cens_rate = 0.08 #censoring rate
  tau_0 = 8 #truncation of survival times
  tau_1 = 10 #end of study
  miss_mech = "MCAR_15" 
  sd = 1
}

if (scenario=='2'){  
  #incidence
  alpha_values = c(0.1,0.5,0.5)
  #latency
  gamma_values = c(-0.345,-0.345)
  rho = 1.45
  lambda = 0.25
  beta_values = -gamma_values*rho
  mi = -log(lambda)/rho
  sigma = 1/rho
  #others
  cens_rate = 0.1
  tau_0 = 8
  tau_1 = 10
  miss_mech = "MCAR_30"
  sd = 1
}

if (scenario=='3'){  
  #incidence
  alpha_values = c(0.1,0.5,0.5)
  #latency
  gamma_values = c(0,-0.345)
  rho = 1.45 
  lambda = 0.25
  beta_values = -gamma_values*rho
  mi = -log(lambda)/rho
  sigma = 1/rho
  #others
  cens_rate = 0.1
  tau_0 = 8
  tau_1 = 10
  miss_mech = "MAR"
  sd = 1
}


real_parameters = c(alpha_values,beta_values)
