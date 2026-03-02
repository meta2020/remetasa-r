## Generate meta-analysis according to the BN-GLMM
## N are sampled from uniform distributions

gendata.1bn.t = function(
    s, 
    n_min, n_max,
    theta,tau,rho,
    Pnmax, Pnmin,
    alpha, beta){
  
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
    mutate(pi=pnorm(alpha+beta*ti))
  
  z=sapply(1:nrow(p.dt2),function(i) rbinom(1,1,p.dt2$pi[i]))
  s.dt=p.dt2[z>0,]
  
  
  res.list = list(
    p.dt=p.dt2,
    s.dt=s.dt,
    Miss=c(M.o=sum(z==0)), 
    p=c(p.e=mean(p.dt2$pi), p.o=mean(z>0))
  )
  return(res.list)
  
}