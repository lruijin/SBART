library(datasets)
air <- datasets::airquality
sum(is.na(air$Ozone))
air <- na.omit(air)

library(Rcpp)
sourceCpp("soft_bart.cpp")
source("SoftBart.R")
source("data_func.R")
id_list = split(1:nrow(air) , 1:5)

X = data.matrix(air[,c("Solar.R","Wind","Temp")])
X = apply(X,2,function(x) (x - min(x))/(max(x) - min(x)))
X = X/ncol(X)
theta_width = 0.1
update_theta_width = 0.08
expTrue = 1
Eta <- Eta_theta <- matrix(NA, 5,3)
MAE <- rep(NA,5)
for(i in 1:5){
  idx = id_list[[i]]
  
  X_test = X[idx,]
  X_train = X[-idx,]
  Y_train= log(air[-idx,"Ozone"])
  Y_test= log(air[idx,"Ozone"])
  acc.rate = 0.6
  while(acc.rate > 0.4|acc.rate < 0.15){
    fit = simbart2(X = X_train, Y = Y_train, X_test= X_test, 
                   hypers = Hypers(X_train, Y_train, num_tree = 10, 
                                   sigma_hat = NULL, sim = T, sparse = F,
                                   prq = rep(0.01,2), M2 = 20),
                   opts = Opts(num_burn = 4000, num_save = 2000, num_print = 500,
                               num_thin = 1, num_update_theta = 400, expTrue = expTrue,
                               theta_width = theta_width, # for beta distribution 5,
                               update_theta_width = update_theta_width, # 0.2 if sparse = F, P = 5
                               update_sigma_mu = FALSE,
                               update_s = FALSE,
                               update_alpha = FALSE,
                               update_beta = FALSE,
                               update_gamma = FALSE,
                               update_tau = TRUE,
                               update_tau_mean = FALSE,
                               update_sigma = TRUE))
    acc.rate = sum(diff(fit$theta[,1])!=0)/nrow(fit$theta)
  }
  
  Eta_theta[i,] <- theta2eta(apply(fit$theta,2,median))
  Eta[i,] <- apply(fit$eta,2,median)
  sum(apply(fit$eta, 2, median)^2)
  mae_simBart_test = mae(fit$y_hat_test_mean, Y_test)
  MAE[i] <- mae_simBart_test
  
}
median(MAE)
Eta
rowSums(Eta^2)
