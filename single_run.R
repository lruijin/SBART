library(Rcpp)
sourceCpp("soft_bart.cpp")
source("SoftBart.R")
source("data_func.R")

single_run <- function(seed, eta, sigma, func, n_train, n_test){
    P = length(eta)
    theta = eta2theta(eta)
    set.seed(1234 + seed)
    sim_data <- gen_data(n_train, n_test, P, func, eta, sigma)
    
    save(sim_data, file = paste0("data/", func,"_P", P, "_sigma",
                                 sub("\\.", "", paste0(sigma/10)),
                                 "_n", n_train, "/sample", i, ".Rdata"))
    
    tStart = Sys.time()
    
    fit = simbart2(X = sim_data$X, Y = sim_data$Y, X_test = sim_data$X_test, 
                   hypers = Hypers(sim_data$X, sim_data$Y, num_tree = 1, 
                                   sigma_hat = NULL, sim = T),
                   opts = Opts(num_burn = 2000, num_save = 3000, 
                               theta_width = 0.3,
                               update_sigma_mu = FALSE,
                               update_s = FALSE,
                               update_alpha = FALSE,
                               update_beta = FALSE,
                               update_gamma = FALSE,
                               update_tau = TRUE,
                               update_tau_mean = FALSE,
                               update_sigma = TRUE))
    tEnd = Sys.time()
    
    Time = as.numeric(tEnd - tStart)
    theta.CI.lower = apply(fit$theta, 2, quantile, 0.025)
    theta.CI.upper = apply(fit$theta, 2, quantile, 0.975)
    eta.CI.lower = apply(fit$eta, 2, quantile,0.025)
    eta.CI.upper = apply(fit$eta, 2, quantile,0.975)
    eta.est = colMeans(fit$eta)
    sse_eta = sum((eta.est - eta)^2)
    theta.est = colMeans(fit$theta)
    sse_theta = sum((theta.est - theta)^2)
    acc.rate = sum(diff(fit$theta[,1])!=0)/nrow(fit$theta)
    rmse_simBart_train = rmse(fit$y_hat_train_mean, sim_data$mu)
    rmse_simBart_test = rmse(fit$y_hat_test_mean, sim_data$mu_test)
    mae_simBart_train = mae(fit$y_hat_train_mean, sim_data$mu)
    mae_simBart_test = mae(fit$y_hat_test_mean, sim_data$mu_test)
    
    save("Time", "eta", "theta", "theta.CI.lower", "theta.CI.upper", 
         "eta.CI.lower", "eta.CI.upper", 
         "eta.est", "sse_eta", 
         "theta.est", "sse_theta",
         "acc.rate", "rmse_simBart_train", "rmse_simBart_test", 
         "mae_simBart_train", "mae_simBart_test", 
         file = paste0("results/", func,"_P", P, "_sigma",
                       sub("\\.", "", paste0(sigma/10)),
                       "_n", n_train, "/result", i, ".Rdata"))
}