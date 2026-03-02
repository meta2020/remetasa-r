## Generate meta-analysis according to the BN-GLMM
## N are sampled from uniform distributions

gendata.1bn = function(
    s, 
    n_min, n_max,
    gr,
    theta,tau,rho,
    Pnmax, Pnmin){
  
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
  ## selective process
  n_min=min(ni) 
  n_max=max(ni)
  
  a1=(qnorm(Pnmax)-qnorm(Pnmin))/(sqrt(n_max)-sqrt(n_min))
  a0=qnorm(Pnmax)-a1*sqrt(n_max)
  
  zi=a0+a1*sqrt(ni)+deltai
  z=1*(zi>0)
  
  pz=pnorm(a0+a1*sqrt(ni))
  pz.s = pz[z>0]
  M= sum((1-pz.s)/pz.s)%>%round()
  p.dt$z=z
  p.dt$pz=pz
  
  s.dt=p.dt[z>0,]
  
  res.list = list(
    p.dt=p.dt,
    s.dt=s.dt,
    a=c(a1=a1,a0=a0),
    Miss=c(M.e=M, M.o=sum(z==0)), 
    p=c(p.e=mean(pz), p.o=mean(z>0))
  )
  return(res.list)
  
}