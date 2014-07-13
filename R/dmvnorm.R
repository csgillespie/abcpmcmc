
##################################################################
# Public functions
##################################################################
#' Multivariate Normal Density 
#' 
#' Calculates mutlivariate normal densities. This is different to the dmvnorm function 
#' in the mvtnorm package in that it takes a matrix for both x and mean. It then calculates 
#' a vector of densities according to dmvnorm(x[i,],mean[i,],sigma,log = FALSE).
#' To aid computation the mahalanobis distances are calculated in parallel using mclapply.
#' @importFrom parallel mclapply
#' @param x A matrix of values
#' @param mean A matrix of means
#' @param sigma A covariance matrix
#' @param log Boolean for whether we want log densities or not
#' @details My own implementation of the multivariate normal density function
#' for increased efficiency for application in this package
#' because there are so many repeated calls to densitymvnorm using a 
#' given sigma matrix it made sense to have one that could take a matrix 
#' of means as well as a matrix of x's and treat them as 
#' paired, returning the density of x[1,] given mean[1,], x[2,] 
#' given mean[2,] also stops repeated inversions of the matrix sigma, and
#' calculates the densities in parallel
#' @return A vector of densities
#' @export
densityMvNorm = function (x, mean, sigma, log = FALSE) 
{
  ## takes a matrix of means rather than a single vector
  if (missing(sigma)) sigma = diag(ncol(x))
  
  if (NCOL(x) != NCOL(sigma)) {
    stop("x and sigma have non-conforming size")
  }
  if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
                   check.attributes = FALSE)) {
    stop("sigma must be a symmetric matrix")
  }
  if (NCOL(mean) != NROW(sigma)) {
    stop("mean and sigma have non-conforming size")
  }
  ## invert matrix before hand so only do this once
  prec = solve(sigma)
  means = lapply(1:dim(mean)[1],function(i){mean[i,]})
  distval = do.call(rbind, 
                     mclapply(means, 
                              mahalanobis,x = x, cov = prec,inverted = TRUE))
  logdet = sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
  logretval = -(ncol(x) * log(2 * pi) + logdet + distval)/2
  
  if (log) 
    return(logretval)
  else 
    return(exp(logretval))
}