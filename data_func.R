f1 <- function(x){
  return(x + x^2)
}

f2 <- function(x){
  return(exp(x))
}

f3 <- function(x){
  return(0.1*x + (sin(0.5*x))^3)
}

f_pool <- function(x, v){
  res = x + (2*x)^2
  res[v == 1] = exp(x[v == 1])-1
  return(res)
}

# plot out the true function
# x.eg <- seq(0,1,len = 1000)
# plot(x.eg, f_pool(x.eg, rep(0,length(x.eg))), type = "l")
# lines(x.eg, f_pool(x.eg, rep(1, length(x.eg))), col = "red")

gen_data <- function(n_train, n_test, P, func, eta, sigma){
  X <- matrix(runif(n_train * P), nrow = n_train)
  mu <- eval(get(func)(X %*% matrix(eta, ncol = 1)))
  
  X_test <- matrix(runif(n_test * P), nrow = n_test)
  mu_test <-  eval(get(func)(X_test %*% matrix(eta,ncol = 1)))
  
  Y <- mu + sigma * rnorm(n_train)
  Y_test <- mu_test + sigma * rnorm(n_test)
  return(list(X = X, Y = Y, mu = mu, X_test = X_test, Y_test = Y_test, mu_test = mu_test,
              eta = eta, sigma = sigma))
}
# plot out data added to the true function

gen_pool_data <- function(n_train, n_test, P, eta, sigma, prob){
  func = "f_pool"
  V = matrix(sample(c(0,1), n_train, replace = T, prob = prob),ncol = 1)
  X <- matrix(runif(n_train * P), nrow = n_train)
  mu <- eval(get(func)(X %*% matrix(eta, ncol = 1), V))
  
  V_test = matrix(sample(c(0,1), n_test, replace = T, prob = prob),ncol = 1)
  X_test <- matrix(runif(n_test * P), nrow = n_test)
  mu_test <-  eval(get(func)(X_test %*% matrix(eta,ncol = 1), V_test))
  
  Y <- mu + sigma * rnorm(n_train)
  Y_test <- mu_test + sigma * rnorm(n_test)
  return(list(X = X, V = V, Y = Y, mu = mu, 
              X_test = X_test, V_test = V_test, Y_test = Y_test, mu_test = mu_test,
              eta = eta, sigma = sigma, prob = prob))
}

rmse <- function(x,y) sqrt(mean((x-y)^2))
nrmse <- function(x,y) sqrt(mean(((x-y)/y)^2))
mae <- function(x,y) mean(abs(x-y))
nmae <- function(x,y) mean(abs((x-y)/y))

# transform eta to theta
eta2theta <- function(eta){
  P = length(eta) - 1
  theta = eta[1:P]
  theta[1] = acos(eta[1])
  dem = sin(theta[1])
  for(i in 2: P){
    theta[i] = acos(eta[i] / dem)
    dem = dem * sin(theta[i])
  }
  return(theta)
}

theta2eta <- function(theta){
  eta = cos(theta[1])
  for(i in 2:length(theta)){
    eta = c(eta,
            cos(theta[i]) * prod(sin(theta[1:(i-1)])))
  }
  eta = c(eta, prod(sin(theta)))
  return(eta)
}
