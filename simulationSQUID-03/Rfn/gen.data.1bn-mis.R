## Generate meta-analysis according to the BN-GLMM
## N are sampled from uniform distributions

gendata.1bn.hedges = function(
    s, 
    n_min, n_max,
    theta,tau,rho,
    Pnmax, Pnmin,
    cutoff=c(0.01,0.05,0.1,0.9), ## <0.05, <0.1, others
    wi=c(0.99,0.99,0.2,0.7,0.9)){
  
  ni = runif( s, min = n_min, max = n_max )%>%round()
  ni = ifelse(ni<10,10,ni)
  
  # generate deltai and thetai
  sigma=matrix(c(tau^2,rho*tau,rho*tau,1),2,2)
  m = MASS::mvrnorm(s,c(theta,0),sigma)
  thetai=m[,1]
  deltai=m[,2]
  
  # generate yi
  yi = sapply(1:s, function(i) rbinom(1, ni[i], plogis(thetai[i])) )
  
  p.dt = data.frame(y=yi,n=ni)
  
  p.dt2 = escalc(measure="PLO", xi=y, ni=n, data=p.dt) %>% 
    mutate(ti=yi/sqrt(vi)) %>%
    mutate(pvalue = 2*pnorm(abs(ti), lower.tail = FALSE)) %>%
    mutate(wi=ifelse(pvalue<cutoff[1], wi[1], 
                     ifelse(pvalue<cutoff[2], wi[2], 
                            ifelse(pvalue<cutoff[3], wi[3],
                                   ifelse(pvalue<cutoff[4], wi[4], wi[5])))))
  
  z=rbinom(nrow(p.dt2),1, p.dt2$wi)
  s.dt=p.dt2[z>0,]
  
  
  res.list = list(
    p.dt=p.dt2,
    s.dt=s.dt,
    Miss=c(M.o=sum(z==0)), 
    p=c(p.e=mean(p.dt2$wi), p.o=mean(z>0))
  )
  return(res.list)
  
}