library(dplyr)
library(MASS)
library(mnormt)
library(foreach)
library(doRNG)
library(metafor)

## true parameters setting ----------
s = c(15,50) ## #oiclif population studies
set = expand.grid(
  t.theta = c(-2), ## true theta
  t.tau = sqrt(c(0.1, 0.3, 0.7)), ## true tau  
  t.rho = c(0.8), ## for population data rho does not matter the estimates
  n.median = c(20,50,100), ## median number of total subjects,
  grp.r = c(1,2), ##, 2group ratio: treat:control
  pmax = 0.99,
  pmin = 0.2
) %>% arrange(t.theta,n.median)
set$ymin = 5
set$ymax = 15
set$nmin = ifelse(set$n.median==20,30,ifelse(set$n.median==50,50,500))
set$nmax = ifelse(set$n.median==20,60,ifelse(set$n.median==50,200,700))
set$p0 = ifelse(set$n.median==20,0.2,ifelse(set$n.median==50,0.1,0.002))