### functions for the prior distributions of states and 
## rates specific to the lotka volterra demo example


## theta prior function
theta_prior = function(N){
  return(cbind(th1 = runif(N,-8,8), 
               th2 = runif(N,-8,8), 
               th3 = runif(N,-8,8)))
}

## x prior function
x_prior <- function(N){
  return(cbind(x1 = rpois(N,50), x2 = rpois(N,100)))
}

prior <- function(N){
  x = x_prior(N)
  th = theta_prior(N)
  ## return a list rather than a matrix ready for use with mclapply
  return(lapply(1:N, function(i){c(x[i,], th[i,])}))
}
