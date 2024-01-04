# library(devtools)
# install_github("theodds/SoftBART")
library(SoftBart)
library(Rcpp)
sourceCpp("soft_bart.cpp")
source("SoftBart.R")
# data used in Park et al (2005) and Dhara et al (2020).
data("airquality")
air = airquality[complete.cases(airquality),]
air = air[,c("Ozone","Solar.R","Wind","Temp")]
X = air[,c("Solar.R","Wind","Temp")]
Y = air$Ozone
X_test = air[,c("Solar.R","Wind","Temp")]

fit_air = simbart2(X = X, Y = Y, X_test = X_test, 
                hypers = Hypers(X, Y, num_tree = 50, sigma_hat = NULL, sim = T),
                opts = Opts(num_burn = 500, num_save = 500, 
                            theta_width = 0.3,
                            update_sigma_mu = FALSE,
                            update_s = FALSE,
                            update_alpha = FALSE,
                            update_beta = FALSE,
                            update_gamma = FALSE,
                            update_tau = TRUE,
                            update_tau_mean = FALSE,
                            update_sigma = TRUE))




str(fit_air)
plot(fit_air)
plot(fit_air$theta[,1], type = "l")
sum(diff(fit_air$theta[,1])!=0)/500

set.seed(1234)

f_fried <- function(x) 10 * sin(pi * x) # * x[,2] + 20 * (x[,3] - 0.5)^2 #+  10 * x[,4] + 5 * x[,5]
gen_data <- function(n_train, n_test, P, theta, sigma) {
  X <- matrix(runif(n_train * P), nrow = n_train)
  
  mu <- f_fried(X %*% matrix(theta,ncol = 1))
  X_test <- matrix(runif(n_test * P), nrow = n_test)
  mu_test <- f_fried(X_test%*% matrix(theta,ncol = 1))
  Y <- mu + sigma * rnorm(n_train)
  Y_test <- mu_test + sigma * rnorm(n_test)
  return(list(X = X, Y = Y, mu = mu, X_test = X_test, Y_test = Y_test, mu_test = mu_test))
}

eta = c(0.61547, 0.7854)
set.seed(1234)

sim_data <- gen_data(50, 100, 2, "f2", eta, 0.01)

# simbart
fit2 = simbart2(X = sim_data$X, Y = sim_data$Y, X_test = sim_data$X_test, 
               hypers = Hypers(sim_data$X, sim_data$Y, num_tree = 50, sigma_hat = NULL, sim = T),
               opts = Opts(num_burn = 2000, num_save = 3000, 
                           theta_width = 0.11,
                           update_sigma_mu = FALSE,
                           update_s = FALSE,
                           update_alpha = FALSE,
                           update_beta = FALSE,
                           update_gamma = FALSE,
                           update_tau = TRUE,
                           update_tau_mean = FALSE,
                           update_sigma = TRUE))



sum(diff(fit2$eta[,1])!=0)/500


par(mfrow = c(1,3))
plot(fit2$eta[,1],type = "l")
abline(h = eta[1], col = 2)

plot(fit2$eta[,2],type = "l")
abline(h = eta[2], col = 2)

plot(fit2$eta[,3],type = "l")
abline(h = 0, col = 2)

fit2$eta[1,]

rmse <- function(x,y) sqrt(mean((x-y)^2))
rmse(fit2$y_hat_test_mean, sim_data$mu_test)
rmse(fit2$y_hat_train_mean, sim_data$mu)
str(fit2)
plot(fit2)


fit <- softbart(X = sim_data$X, Y = sim_data$Y, X_test = sim_data$X_test,
                hypers = Hypers(sim_data$X, sim_data$Y, num_tree = 50, temperature = 1),
                opts = Opts(num_burn = 500, num_save = 500, update_tau = TRUE))
rmse(fit$y_hat_test_mean, sim_data$mu_test)
rmse(fit$y_hat_train_mean, sim_data$mu)
plot(fit)
