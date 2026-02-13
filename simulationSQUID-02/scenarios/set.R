library(dplyr)
library(MASS)
library(mnormt)
library(foreach)
library(doRNG)
library(metafor)
s = c(15,50)
## true parameters setting ----------
 ## #the number of population studies
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


##
## calculate alpha in P(alpha+beta*t), beta=2 by Hu et al (2024)
##

beta=2
fa = function(alpha) mean(pnorm(alpha + beta*ti)) -0.7

rtimes = 1000
alpha.1bn = NULL
alpha.2hn = NULL
alpha.2bn = NULL
S=50

ncores = min(120, parallel::detectCores()-1)


cl = parallel::makeCluster(ncores, "SOCK")
doSNOW::registerDoSNOW(cl)
set.seed(2025)
for(i in 1:nrow(set)){#
  al.1bn = foreach(r=1:rtimes, .combine = rbind,.packages=c("mnormt","dplyr","metafor"))  %dorng%  {
    
    set.gr = set[i,]
    pmax = set.gr$pmax
    pmin = set.gr$pmin
    
    ## create 1BN datasets
    plist = gendata.1bn(
      s=S, 
      n_min = round(set.gr$nmin/2), n_max=round(set.gr$nmax/2),
      theta=set.gr$t.theta,
      tau=set.gr$t.tau,
      rho=set.gr$t.rho,
      Pnmax = pmax, Pnmin = pmin)
    
    
    ## population BN data
    pdata.1bn = plist$p.dt
    pdata.1bn2 = escalc(measure="PLO", xi=y, ni=n, data=pdata.1bn)
    ti = with(pdata.1bn2, yi/sqrt(vi))
    
    uniroot(fa, c(-5,5), extendInt = "yes")$root
  }
  
  alpha.1bn =c(alpha.1bn, mean(al.1bn, na.rm = TRUE))

  
  al.2hn = foreach(r=1:rtimes, .combine = rbind,.packages=c("mnormt","dplyr","metafor"))  %dorng%  {
    
    set.gr = set[i,]
    pmax = set.gr$pmax
    pmin = set.gr$pmin
    
      ## create HN datasets
    plist = gendata.2hn(
      s=S, 
      n_min = set.gr$nmin, n_max=set.gr$nmax,
      gr=set.gr$grp.r,
      theta=set.gr$t.theta,
      tau=set.gr$t.tau,
      rho=set.gr$t.rho,
      y_min=set.gr$ymin,y_max=set.gr$ymax,
      Pnmax = pmax, Pnmin = pmin)
    
    
    ## population 2HN data
    pdata.2hn = plist$p.dt
    pdata.2hn2 = escalc(measure="OR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0, data=pdata.2hn)
    ti = with(pdata.2hn2, yi/sqrt(vi))
    
    uniroot(fa, c(-5,5), extendInt = "yes")$root
  }
  
  alpha.2hn =c(alpha.2hn, mean(al.2hn, na.rm = TRUE))
  
  al.2bn = foreach(r=1:rtimes, .combine = rbind,.packages=c("mnormt","dplyr","metafor"))  %dorng%  {
    
    set.gr = set[i,]
    pmax = set.gr$pmax
    pmin = set.gr$pmin
    
    ## create 2BN datasets
    plist = gendata.2bn(
      s=S, 
      n_min = set.gr$nmin, n_max=set.gr$nmax,
      p0=set$p0,
      gr=set.gr$grp.r,
      theta=set.gr$t.theta,
      tau=set.gr$t.tau,
      rho=set.gr$t.rho,
      Pnmax = pmax, Pnmin = pmin)
    
    
    ## population BBN data
    pdata.2bn = plist$p.dt
    pdata.2bn2 = escalc(measure="OR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0, data=pdata.2bn)
    ti = with(pdata.2bn2, yi/sqrt(vi))
    
    uniroot(fa, c(-5,5), extendInt = "yes")$root
  }
  alpha.2bn =c(alpha.2bn, mean(al.2bn, na.rm = TRUE))
}
parallel::stopCluster(cl)

set$alpha.2hn = alpha.2hn
set$alpha.2bn = alpha.2bn
set$alpha.1bn = alpha.1bn
set$beta = rep(beta, nrow(set))

save(set, file = "set.RData")
