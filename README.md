abcpmcmc
========
[![Build Status](https://travis-ci.org/csgillespie/abcpmcmc.png?branch=master)](https://travis-ci.org/csgillespie/abcpmcmc)

This package illustrates example one of the paper

>Owen, JR, Gillespie, CS, Wilkinson, DJ, [Scalable
>  Inference for Markov Processes with Intractable Likelihoods](http://dx.doi.org/10.1007/s11222-014-9524-7), 
>  Statistics and Computing, 2014.

Source code for the other examples can be obtained by contacting Jamie Owen <j.r.owen@newcastle.ac.uk)>

Installation
------------

This package is only hosted on github. To install the package, 
```r
install.packages("devtools")
library(devtools)
install.packages(c("smfsb", "mvtnorm"))
install_github("abcpmcmc", "csgillespie")
```

Note Windows users have to first install [Rtools](http://cran.rstudio.com/bin/windows/Rtools/).


