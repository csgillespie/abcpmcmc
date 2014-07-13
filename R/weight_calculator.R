#' Weight Calculator
#' 
#' A function to create a weight calclulator object for use in 
#' calculating particle weights in ABC SMC.
#' This example is specific to the Lotka Volterra model and priors
#' as a demonstration
#' @param weights A vector of weights that will be updated
#' @param samples A matrix of samples that have the associated weights
#' @param sigma A covariance matrix according to a Gaussian proposal that has been used 
#' @return A function which calculates weights for a matrix of proposals
#' @export
weightCalculator = function(weights,samples,sigma){
  w = weights/sum(weights)
  s = samples[,1:3]
  o = sigma
  return(function(proposals){
    k = densityMvNorm(s, mean=proposals, sigma = o, log = FALSE)
    temp = w*k
    weights = (1.0/4096.0)/apply(temp,1,sum)
    ## weights are 0 for parameter vectors outside prior support
    weights[which(abs(proposals[,1]) > 8 |
                    abs(proposals[,2]) > 8 | 
                    abs(proposals[,3]) > 8)] = 0
    return(weights)
  })
}