f1 <- function(x){
  return(x + x^2)
}

f2 <- function(x){
  return(exp(x))
}

f3 <- function(x){
  return(0.1*x + (sin(0.5*x))^3)
}

gen_data <- function(n_train, n_test, P, func, eta, sigma){
  X <- matrix(runif(n_train * P), nrow = n_train)
  mu <- eval(get(func)(X %*% matrix(eta, ncol = 1)))
  
  X_test <- matrix(runif(n_test * P), nrow = n_test)
  mu_test <-  eval(get(func)(X_test %*% matrix(eta,ncol = 1)))
  
  Y <- mu + sigma * rnorm(n_train)
  Y_test <- mu_test + sigma * rnorm(n_test)
  return(list(X = X, Y = Y, mu = mu, X_test = X_test, Y_test = Y_test, mu_test = mu_test))
}

rmse <- function(x,y) sqrt(mean((x-y)^2))
mae <- function(x,y) mean(abs(x-y))

# transform eta to theta
eta2theta <- function(eta){
  P = length(eta) - 1
  theta = eta[1:P]
  theta[1] = acos(eta[1])
  dem = sin(theta[1])
  for(i in 2: P){
    theta[i] = acos(theta[i] / dem)
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
}
