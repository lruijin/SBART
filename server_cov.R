setwd('/storage1/fs1/r.lu/Active/Projects/simBart_cov')
source("single_run_cov.R")

args <- commandArgs(trailingOnly = T)
prob1 <- as.numeric(args[1])
if(!prob1 %in% c(0.1,0.3,0.5)){
  stop("The current simulation settings only allow for prob1 = 0.1, 0.3 and 0.5")
}
sigma <- as.numeric(args[2])
if(!sigma %in% c(0.01,0.1,0.3,0.5)){
  stop("The current simulation seetings only allow sigma = 0.01, 0.1, 0.3 and 0.5")
}
n_train <- as.integer(args[3])
if(!n_train %in% c(50,100)){
  stop("The training dataset is only allow n = 50 and 100")
}


eta = rep(1/sqrt(3),3)
eta = eta/sum(eta^2)
# sigma = 0.01
prob = c(prob1, 1-prob1)
# n_train = 100
n_test = 20
sparse = F

update_theta_width = 0.1 #0.1 for f1
# prq = rep(0.02, null_p + 2) # rep(0.01, 4) for f1
# expTrue = 2
theta_width = 0.16
load.data = F
for(i in 1:200){
  update = single_run_cov(i, eta, sigma, n_train, n_test, theta_width = theta_width,
                      load.data = load.data, update_theta_width = update_theta_width,
                      sparse = sparse, prob = prob, sim_cov = F)
  theta_width = update
}

load.data = T
theta_width = 0.1
for(i in 1:200){
  update = single_run_cov(i, eta, sigma, n_train, n_test, theta_width = theta_width,
                          load.data = load.data, update_theta_width = update_theta_width,
                          sparse = sparse, prob = prob, sim_cov = T)
  theta_width = update
}



