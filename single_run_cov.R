library(Rcpp)
sourceCpp("soft_bart.cpp")
source("SoftBart.R")
source("data_func.R")

single_run_cov <- function(seed, eta, sigma, n_train, n_test, theta_width = 0.2,
                           load.data = F, update_theta_width = 0.2, sparse = F, 
                           prq = rep(0.01,length(eta)-1), expTrue = 3, prob = c(0.5, 0.5),
                           sim_cov = F){
  P = length(eta)
  theta = eta2theta(eta)
  set.seed(1234 + seed)
  
  prob.name = as.character(round(prob[1]*100))
  if(load.data){
    load(paste0("data_cov/Prob",prob.name, "_P", P, "_sigma",
                sub("\\.", "", paste0(sigma*10)),
                "_n", n_train, "/sample", seed, ".Rdata"))
  }else{
    sim_data <- gen_pool_data(n_train, n_test, P, eta, sigma, prob)
    # sim_data.plot <- data.frame(x = sim_data$X %*% eta, 
    #                             y = sim_data$Y,
    #                             ytrue = sim_data$mu,
    #                             v = factor(sim_data$V,labels = c("Group 1", "Group 2")))
    # library(ggplot2)
    # ggplot(sim_data.plot, aes(x = x, group = v)) +
    #   geom_line(aes(y = ytrue, col = v)) +
    #   geom_point(aes(y = y, col = v))+
    #   theme_classic() +
    #   labs(y = "Response",
    #        x = "Single index",
    #        col = "Covariate") + 
    #   theme(legend.position = "bottom",
    #         axis.title = element_text(size = 20),
    #         axis.text = element_text(size = 15),
    #         legend.title = element_text(size = 20),
    #         legend.text = element_text(size = 15))
    save(sim_data, file = paste0("data_cov/Prob",prob.name, "_P", P, "_sigma",
                                 sub("\\.", "", paste0(sigma*10)),
                                 "_n", n_train, "/sample", seed, ".Rdata"))
  }
  eta = sim_data$eta
  sigma = sim_data$sigma
  P = length(eta)
  theta = eta2theta(eta)
  
  acc_tmp = 0.6
  if(!sim_cov){
    sim_data$V = NULL
    sim_data$V_test = NULL
  }
  for(j in 1:5){
    cat("Try: ", j, "\n")
    tStart = Sys.time()
    fit = simbart2(X = sim_data$X, V = sim_data$V, Y = sim_data$Y, X_test = sim_data$X_test,
                   V_test = sim_data$V_test,
                   hypers = Hypers(sim_data$X, sim_data$V, sim_data$Y, num_tree = 10, 
                                   sigma_hat = NULL, sim = T, sparse = sparse,
                                   prq = prq, M2 = 15),
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
    tEnd = Sys.time()
    acc.rate = sum(diff(fit$theta[,1])!=0)/nrow(fit$theta)
    # par(mfrow = c(2,1))
    # plot(fit$theta[,1],type = "l")
    # plot(fit$eta[,1],type = "l")
    # sum(colMeans(fit$eta)^2)
    if(acc.rate > 0.15 & acc.rate < 0.4){
      fit_tmp = fit
      acc_tmp = acc.rate
      break
    }else{
      if(acc.rate < 0.15){
        theta_width = theta_width / (1+update_theta_width)
      }else if(acc.rate > 0.4){
        theta_width = theta_width * (1+update_theta_width)
      }
      if(max(abs(acc.rate - c(0.15, 0.4))) < max(abs(acc_tmp-c(0.15,0.4)))){
        fit_tmp = fit
        acc_tmp <- acc.rate
      }
    }
  }
  
  eta.true = eta
  theta.true = theta
  acc.rate <- acc_tmp
  fit = fit_tmp
  Time = as.numeric(tEnd - tStart)
  
  theta.sample = fit$theta
  eta.sample = fit$eta
  
  yhat_train = fit$y_hat_train
  yhat_test = fit$y_hat_test
  
  ytrue_train = sim_data$mu
  ytrue_test = sim_data$mu_test
  
  ytrain <- sim_data$Y
  ytest <- sim_data$Y_test
  
  if(sim_cov){
    save("Time", "eta.true", "theta.true", "acc.rate", 
         "theta.sample", "eta.sample", "theta_width",
         "yhat_train", "yhat_test", 
         "ytrue_train", "ytrue_test","ytrain","ytest",
         file = paste0("results_cov/Prob",prob.name, "_P", P, "_sigma",
                       sub("\\.", "", paste0(sigma*10)),
                       "_n", n_train, "/result", seed, ".Rdata"))
  }else{
    save("Time", "eta.true", "theta.true", "acc.rate", 
         "theta.sample", "eta.sample", "theta_width",
         "yhat_train", "yhat_test", 
         "ytrue_train", "ytrue_test","ytrain","ytest",
         file = paste0("results_standard1/Prob",prob.name, "_P", P, "_sigma",
                       sub("\\.", "", paste0(sigma*10)),
                       "_n", n_train, "/result", seed, ".Rdata")) 
  }
  return(theta_width)
  
}
