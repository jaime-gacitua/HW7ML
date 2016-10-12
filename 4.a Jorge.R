  rm(list=ls())
  require(geoR)
  # 
  m = 8
  A = c(62,60,63,59,NA,NA,NA,NA)
  B = c(63,67,71,67,68,68,NA,NA)
  C = c(68,66,71,67,68,68,NA,NA)
  D = c(56,62,60,61,63,64,63,59)
  J = 4
  dos_n = 200
  n0 = 100
  iterations = dos_n + n0
  
  data = rbind(A,B,C,D)
  
  # Vector of theta_j, theta_j_hat and V_theta_j
  theta_j = rep(1,4)
  theta_j_hat = rep(1,4)
  V_theta_j = rep(1,4)
  
  n = length(data[!is.na(data)])
  
  # Getting starting points for thetas
  for (i in 1:J){
    theta_j[i] <- mean(na.omit(data[i,]))
  }
  # Starting point for mu
  mu = mean(theta_j)
  
  for (k in 1:100){
  # 3
  temp_sum = 0
    for (j in 1:J) {
      temp_sum = temp_sum + sum((na.omit(data[j,])-theta_j[j])^2)
    }
  sig_hat_sq = temp_sum*(1/n)
  sig_sq = rinvchisq(1, df = n, scale = sig_hat_sq)
  
  # 4 
  
  tao_hat_sq = (1/(J-1))*sum((theta_j-mu)^2)
  tao_sq = rinvchisq(1, df = (J-1), scale = tao_hat_sq)
  
  # 1 MAYBE THETA J MUST BE INSIDE OTHER LOOP
  for (i in 1:J){
    theta_j_hat[i] = ((1/tao_sq)*mu+(length(na.omit(data[i,]))/sig_sq)*mean(na.omit(data[i,])))/
      ((1/tao_sq)+(length(na.omit(data[i,]))/sig_sq))
    V_theta_j[i] = 1/((1/tao_sq)+(length(na.omit(data[i,]))/sig_sq))
    theta_j[i] = rnorm(1,mean = theta_j_hat[i], sd = sqrt(V_theta_j[i]))
  }
  
  for (i in 1:J){
    
  }
  
  # 2
  
  mu_hat = mean(theta_j)
  mu = rnorm(1,mean = mu_hat,sd = sqrt(tao_sq/J))
  }
