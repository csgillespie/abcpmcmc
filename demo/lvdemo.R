## necessary library for parallelisation
## set the number of cores available for use to the maximum detected by the system
library("parallel")
nCores = detectCores()
options(mc.cores = nCores)


## package associated with smfsb, contains the data set and a c implemented discretised gillespie simulator
## load in the reference data sets which are part of the smfsb package
library("smfsb")
data(LVdata)

lvdata = as.matrix(LVnoise10)
tvec = seq(0,30,2)

## theta prior function
theta_prior <- function(N){
  return(cbind(th1 = rnorm(N,1,2),th2 = rnorm(N,-5.3,2),th3 = rnorm(N,-0.5,2)))
}

## x prior function
x_prior <- function(N,t0 = 0){
  return(cbind(x1 = rpois(N,50),x2 = rpois(N,100)))
}

prior <- function(N){
  x = x_prior(N)
  th = theta_prior(N)
  ## return a list rather than a matrix ready for use with mclapply
  return(lapply(1:N,function(i){c(x[i,],th[i,])}))
}

## eucliden metric
euc = metric(LVnoise10)

# first tolerance
e1 = 750

## create the abcSmc runnner function
abcrunner = abcSmc(epsilon1 = e1, simulator = stepLVc, times = tvec, distance = euc, prior = prior)

## run the output for a few populations to keep 1000 samples
th = abcrunner(nSamples = 1000, nDist = 10, alpha = 0.3)

## set up the data likelihood function for the particle filter
dataLikelihood <- function(x,y,log=TRUE){
  ll=sum(dnorm(y,x,noiseSD,log=TRUE))
  if (log) return(ll)
  else return(exp(ll))
}

iters = 250
thin = 50
proposalVariance = var(th[,1:3])

chainRunner = mcmcChain(proposalVariance,dataLikelihood,100,x_prior,LVnoise10,stepLVc)

##sample from output of the ABC some intialisations
initialisations = th[sample(1:1000,nCores,replace = TRUE),1:3]
initialisations = lapply(1:nrow(initialisations),function(i){initialisations[i,]})

## run the chains in parallel
output = mclapply(initialisations,chainRunner,nIters = iters, nThin = thin)
