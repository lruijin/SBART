source("single_run.R")

eta = c(rep(1/sqrt(3),3),0,0)
eta = eta/sum(eta^2)
sigma = 0.1
func = "f2"
n_train = 50
n_test = 20

for(i in 1:100){
    single_run(i, eta, sigma, func, n_train, n_test)
}




