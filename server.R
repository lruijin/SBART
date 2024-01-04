sourceCpp("soft_bart.cpp")
source("SoftBart.R")
source("data_func.R")

eta = rep(1/sqrt(3),3)
eta = eta/sum(eta^2)
sigma = 0.01
func = "f2"
n_train = 80
n_test = 20
P = length(eta)


set.seed(1234)
sim_data <- gen_data(n_train, n_test, P, func, eta, sigma)

# simbart
tStart = Sys.time()

fit = simbart2(X = sim_data$X, Y = sim_data$Y, X_test = sim_data$X_test, 
                hypers = Hypers(sim_data$X, sim_data$Y, num_tree = 1, 
                                sigma_hat = NULL, sim = T),
                opts = Opts(num_burn = 2000, num_save = 3000, 
                            theta_width = 0.12,
                            update_sigma_mu = FALSE,
                            update_s = FALSE,
                            update_alpha = FALSE,
                            update_beta = FALSE,
                            update_gamma = FALSE,
                            update_tau = TRUE,
                            update_tau_mean = FALSE,
                            update_sigma = TRUE))
tEnd = Sys.time()

# plot(fit$y_hat_train_mean, sim_data$mu)
# plot(fit$y_hat_test_mean, sim_data$mu_test)
Time = as.numeric(tEnd - tStart)
CI.lower = apply(fit$eta, 2, quantile, 0.025)
CI.upper = apply(fit$eta, 2, quantile, 0.975)
acc.rate = sum(diff(fit$eta[,1])!=0)/nrow(fit$eta)
rmse_simBart_train = rmse(fit$y_hat_train_mean, sim_data$mu)
rmse_simBart_test = rmse(fit$y_hat_test_mean, sim_data$mu_test)
mae_simBart_train = mae(fit$y_hat_train_mean, sim_data$mu)
mae_simBaet_test = mae(fit$y_hat_test_mean, sim_data$mu_test)
