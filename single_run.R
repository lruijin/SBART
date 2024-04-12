library(Rcpp)
sourceCpp("soft_bart.cpp")
source("SoftBart.R")
source("data_func.R")

single_run <- function(seed, eta, sigma, func, n_train, n_test, theta_width = 0.2,
                       load.data = F, update_theta_width = 0.2, sparse = F, 
                       prq = rep(0.01,P-1), expTrue = 3){
    P = length(eta)
    theta = eta2theta(eta)
    set.seed(1234 + seed)
    
    
    if(load.data){
        load(paste0("data/", func,"_P", P, "_sigma",
                    sub("\\.", "", paste0(sigma*10)),
                    "_n", n_train, "/sample", seed, ".Rdata"))
    }else{
        sim_data <- gen_data(n_train, n_test, P, func, eta, sigma)
        save(sim_data, file = paste0("data/", func,"_P", P, "_sigma",
                                     sub("\\.", "", paste0(sigma*10)),
                                     "_n", n_train, "/sample", seed, ".Rdata"))
    }
    
    for(j in 1:5){
        acc_tmp = 0.6
        cat("Try: ", j, "\n")
        tStart = Sys.time()
        fit = simbart2(X = sim_data$X, Y = sim_data$Y, X_test = sim_data$X_test, 
                       hypers = Hypers(sim_data$X, sim_data$Y, num_tree = 10, 
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
        par(mfrow = c(2,1))
        plot(fit$theta[,1],type = "l")
        plot(fit$eta[,1],type = "l")
        sum(colMeans(fit$eta)^2)
        if(acc.rate > 0.15 & acc.rate < 0.4){
            fit_tmp = fit
            acc_tmp = acc.rate
            break
        }else{
            if(max(abs(acc.rate - c(0.15, 0.4))) < max(abs(acc_tmp-c(0.15,0.4)))){
                fit_tmp = fit
                acc_tmp <- acc.rate
            }
        }
    }
    
    acc.rate <- acc_tmp
    fit = fit_tmp
    Time = as.numeric(tEnd - tStart)
    theta.CI.lower = apply(fit$theta, 2, quantile, 0.025)
    theta.CI.upper = apply(fit$theta, 2, quantile, 0.975)
    theta.est = apply(fit$theta,2,mean)
    
    eta.est = theta2eta(theta.est)
    eta.est = colMeans(fit$eta)
    sse_eta = sum((eta.est - eta)^2)
    
    eta.CI.lower = apply(fit$eta, 2, quantile, 0.025)
    eta.CI.upper = apply(fit$eta, 2, quantile, 0.975)
    sse_theta = sum((theta.est - theta)^2)
    # acc.rate = sum(diff(fit$theta[,1])!=0)/nrow(fit$theta)
    rmse_simBart_train = rmse(fit$y_hat_train_mean, sim_data$mu)
    rmse_simBart_test = rmse(fit$y_hat_test_mean, sim_data$mu_test)
    mae_simBart_train = mae(fit$y_hat_train_mean, sim_data$mu)
    mae_simBart_test = mae(fit$y_hat_test_mean, sim_data$mu_test)
    
    if(sparse){
        sel <- colMeans(fit$delta)
        theta4eta = colMeans(fit$theta)
        theta4eta[sel < 0.5] = pi/2
        eta_hat = theta2eta(theta4eta)
        sse_eta_hat = sum((eta_hat - eta)^2)
        save("Time", "eta", "theta", "theta.CI.lower", "theta.CI.upper", 
             "eta.CI.lower", "eta.CI.upper", 
             "eta.est", "sse_eta", 
             "theta.est", "sse_theta",
             "acc.rate", "rmse_simBart_train", "rmse_simBart_test", 
             "mae_simBart_train", "mae_simBart_test", "sel",
             file = paste0("results_sparse/", func,"_P", P, "_sigma",
                           sub("\\.", "", paste0(sigma*10)),
                           "_n", n_train, "/result", seed, ".Rdata"))
    }else{
        save("Time", "eta", "theta", "theta.CI.lower", "theta.CI.upper", 
             "eta.CI.lower", "eta.CI.upper", 
             "eta.est", "sse_eta", 
             "theta.est", "sse_theta",
             "acc.rate", "rmse_simBart_train", "rmse_simBart_test", 
             "mae_simBart_train", "mae_simBart_test", 
             file = paste0("results/", func,"_P", P, "_sigma",
                           sub("\\.", "", paste0(sigma*10)),
                           "_n", n_train, "/result", seed, ".Rdata")) 
    }
    
}
