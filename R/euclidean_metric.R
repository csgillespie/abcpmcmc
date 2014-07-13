## a simple euclidean metric function for use with the abc pmcmc hybrid demo

################################################################
# Public functions
################################################################
#' Eulidean Metric
#'
#' Creates a function which calculates the Euclidean distances between some simulated data and an observation time index 
#'
#' @param observed The set of observed data
#' @return A function which calculates the Euclidean distance given some simulated data and a row/time index
#' @export
metric <- function(observed){
  xobs = observed
  return(function(xsim, i){
    return(sum((xsim - xobs[i,])*(xsim-xobs[i,])))
  })
}
