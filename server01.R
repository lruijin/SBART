source("single_run.R")

eta = c(rep(1/sqrt(3),3),rep(0,2))
eta = eta/sum(eta^2)
sigma = 0.01
func = "f1"
n_train = 100
n_test = 20
theta_width = 0.2# 0.2 for f1
load.data = T
sparse = T
update_theta_width = 0.1 #0.1 for f1
prq = rep(0.02, 4) # rep(0.01, 4) for f1
expTrue = 2

for(i in 1:200){
    single_run(i, eta, sigma, func, n_train, n_test, theta_width = theta_width,
               load.data = load.data, update_theta_width = update_theta_width, 
               sparse = sparse, prq  = prq)
}




