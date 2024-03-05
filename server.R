source("single_run.R")

eta = c(rep(1/sqrt(3),3),rep(0,2))
eta = eta/sum(eta^2)
sigma = 0.3
func = "f1"
n_train = 100
n_test = 20
theta_width = 0.15

for(i in 1:100){
    single_run(i, eta, sigma, func, n_train, n_test, theta_width = theta_width)
}




