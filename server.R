source("data_func.R")
eta = c(0.61547, 0.7854)
eta = eta/sum(eta^2)
sigma = 0.01
func = "f2"
n_train = 80
P = length(eta)

set.seed(1234)

sim_data <- gen_data(n_train, 20, P, func, eta, sigma)

# simbart
tStart = Sys.time()
fit = simbart2(X = sim_data$X, Y = sim_data$Y, X_test = sim_data$X[1,], 
                hypers = Hypers(sim_data$X, sim_data$Y, num_tree = 50, 
                                sigma_hat = NULL, sim = F),
                opts = Opts(num_burn = 0, num_save = 50, 
                            update_sigma_mu = FALSE,
                            update_s = FALSE,
                            update_alpha = FALSE,
                            update_beta = FALSE,
                            update_gamma = FALSE,
                            update_tau = TRUE,
                            update_tau_mean = FALSE,
                            update_sigma = TRUE))
tEnd = Sys.time()

plot(fit2$y_hat_train_mean, sim_data$mu)
plot(fit2$y_hat_test_mean, sim_data$mu_test)
Time = as.numeric(tEnd - tStart)
CI.lower = apply(fit2$eta, 2, quantile, 0.025)
CI.upper = apply(fit2$eta, 2, quantile, 0.975)
acc.rate = sum(diff(fit2$eta[,1])!=0)/nrow(fit2$eta)
rmse_simBart_train = rmse(fit2$y_hat_train_mean, sim_data$mu)
rmse_simBart_test = rmse(fit2$y_hat_test_mean, sim_data$mu_test)
mae_simBart_train = mae(fit2$y_hat_train_mean, sim_data$mu)
mae_simBaet_test = mae(fit2$y_hat_test_mean, sim_data$mu_test)
