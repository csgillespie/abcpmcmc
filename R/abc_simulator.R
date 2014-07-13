#' Abc Simulator wrapper
#'
#' Wraps up a model simulator function for use with ABC 
#' @inheritParams abcSmc
#' @return A function for building the ABC wrapped simulators
#' @export
abcSimulatorWrapper <- function(distance, simulator, times){
  deltas = c(0,diff(times))
  sim = simulator
  metric = distance
  
  ## returns a function which builds a simulatorwrapper object based on the 
  ## observation regime, simulator and metric function being constant
  ##, which itself returns a function
  return(function(epsilon){
    tolerance = epsilon
    ## returns a function which takes in some states and log rate parameters as a single vector and simulates
    return(function(pars){
      d = 0
      deltat = 0
      ## to further speed up computation, we dont simulate for 
      ## any parameter values whihc are proposed outside of the 
      ## region of prior support, specific to the LV model
      if(any(abs(pars[3:5]) > 8)) return()
      for(i in 1:length(deltas)){
        t = deltat
        deltat = deltat + deltas[i]
        ## the following is specific to the Lotka Volterra model with 2 species and 3 parameters
        x = sim(x0 = pars[1:2],t,deltat,th = exp(pars[3:5]))
        d = d + metric(x,i)
        ## check if the simulation has already exceeded the tolerance, if so return nothing and stop forward simulation
        if(sqrt(d) > tolerance)
          return()
      }
      return(c(pars[3:5], d = sqrt(d))) ## the subsetting of pars is specific to LV example
    })
  })
}
