### an mcmc chain builder

################################################
# Public functions
################################################
#' MCMC Chain
#' 
#' A function which does the set up work for a particle MCMC chain, returning a function which can be used to run 
#' the chain given some initial parameter values, a number of iterations and an amount to thin samples by.
#' @importFrom mvtnorm rmvnorm
#' @importFrom smfsb pfMLLik
#' @param proposalVariance The covariance matrix to be used in the Gaussian random walk proposal
#' @param dataLikelihood A function which calculates the data likelihood
#' @param nParticles The number of particles to be used in the particle filter
#' @param xPrior A function which returns samples from the prior distribution on the states
#' @param data The observed data, as a matrix
#' @param simulator The model simulator to be used for advancing the state of the system given some rate parameters
#' @return A fucntion whihc runs a particle MCMC chain 
#' @export

mcmcChain <- function(proposalVariance, dataLikelihood, 
                      nParticles, xPrior, data, simulator){
  tune = proposalVariance
  ## build a particle filter which estimates the marginal log likelihood, see smfsb for details
  pf = pfMLLik(nParticles,  simx0=xPrior, 
               t0=0, stepFun=simulator, 
               dataLik=dataLikelihood, data=data)
  return(function(initial, nIters, nThin){
    currLl = 1e99
    currTh = initial
    message(currTh)
    ## set up a matrix to store accepted samples
    sample = matrix(0, nrow = nIters, ncol = length(initial))
    for(i in 1:nIters){
      for(j in 1:nThin){
        r = rmvnorm(1,sigma = tune)
        message(paste("current Theta = ",currTh,
                      "random element = ",r))
        propTh = currTh + r
        propLl = pf(exp(propTh))
        if(log(runif(1)) < (propLl - currLl)){
          currTh = propTh
          currLl = propLl
        }
      }
      sample[i,] = currTh
    }
    return(sample)
  })
}