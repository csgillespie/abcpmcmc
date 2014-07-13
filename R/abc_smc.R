#' ABC SMC
#'
#' A function which returns an ABC SMC algorithm function given an initial tolerance, a 
#' model simulator, a vector of observation times, a distance calculator function which contains 
#' the observed data information and a prior distribution for both the rates and states.
#' 
#' @importFrom parallel mclapply 
#' @importFrom mvtnorm rmvnorm
#' @param epsilon1 The initial tolerance value for the ABC SMC algorithm
#' @param simulator A model simulator
#' @param times A vector of observation times of the data
#' @param distance A metric function
#' @param prior A function which gives returns prior samples of the states and rate parameters
#' @return A function which runs the ABC SMC algorithm for a desired accepted sample size, over a desired number of populations for a given adaptive tolerance selected as the alpha ntile of the distribution of distances
#' @export
abcSmc = function(epsilon1, simulator, times, distance,prior){
  ## wraps up the simulator using the ABC wrapper function
  simfunction = abcSimulatorWrapper(distance,simulator,times)
  priorfunction = prior
  sim = simfunction(epsilon1)
  
  return(function(nSamples,nDist,alpha){
    message("running abc smc")
    for(population in 1:nDist){
      samples = 0
      batch = 10000
      ## empty matrix to store results
      out = c()
      
      while(samples < nSamples){
        ## simulate forward and keep any non null in a matrix
        out = rbind(do.call(rbind, Filter(Negate(is.null),mclapply(priorfunction(batch),sim))),out)
        ## if none are retained in a batch set samples = 0
        samples = ifelse(!is.null(nrow(out)),nrow(out),0)
        message(paste(population,"\t",samples,sep = ""))
      }
      ## the number of parameters is one less since kept samples have attached their distance
      pars = ncol(out) - 1
      
      message("updating weights")
      if(population == 1){
        weights = rep(1/samples,samples) ## normalised weights
      }else{
        ## construct a new weight calculator function based on current samples and weigths
        weightCalc = weightCalculator(weights,post,sigma)
        weights = weightCalc(out[,1:pars])
        weights = weights/sum(weights)
      }
      post = out
      ## work out desired new tolerance
      tol = quantile(out[,pars+1],alpha)
      ## the samples which beat new tolerance
      indicies = which(out[,pars+1] < tol)
      good = out[indicies,]
      goodweights = weights[indicies]
      goodweights = goodweights/sum(goodweights) ## normalise
      ## calculate tuning parameters for innovations
      message("Calculating new proposal covariance")
      sigma = matrix(0,nrow = pars, ncol = pars)
      
      for(i in 1:samples){
        for(j in 1:length(indicies)){
          v = good[j,1:pars] - out[i,1:pars]
          v = v %*% t(v)
          sigma = sigma + weights[i]*goodweights[j]*v
        }
      }
      
      ## if ess < nSamples/2 resample
      message("checking ESS")
      s = sum(weights)
      if(s*s/sum(weights*weights) < samples/2){
        message("resampling")
        index = sample(1:samples,nSamples,weights,replace = TRUE)
        weights = weights[index]
        post = out[index,]
      }
      if(population == 1){
        df = cbind(post[,1],population)
        colnames(df) <- c("th1","population")
        df = as.data.frame(df)
      }
      else{
        x = cbind(post[,1],population)
        colnames(x) <- c("th1","population")
        df = rbind(df,x)
      } 
      ## construct new abc simulator based on new tolerance
      sim = simfunction(tol)
      message(paste("New tolerance = ",tol,sep = ""))
      ## construct a new "prior" function based on proposals from current samples
      priorfunction = function(N){
        x = x_prior(N)
        i = sample(1:dim(post)[1],N,weights,replace=TRUE)
        th = post[i,pars] + rmvnorm(N,sigma = sigma)
        return(lapply(1:N,function(i){c(x[i,],th[i,])}))
      }
    }
    ## final sample
    i = sample(1:nrow(post),nSamples,weights,replace=TRUE)
    return(list("post" = post, "df" = df))
  })
}
