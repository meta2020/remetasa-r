## Generate meta-analysis according to the HN-GLMM 
## N are sampled from uniform distributions


gendata.2hn = function(
    s, 
    n_min, n_max,
    gr,
    theta,tau,rho,
    y_min,y_max,
    Pnmax, Pnmin){
  
  n = runif( s, min = n_min, max = n_max )
  n = ifelse(n<20,20,n)
  
  n0i = (n*(1 / (1 + gr))) %>% round()
  n1i = (n-n0i)%>% round() # #treatment subjects
  ni=n0i+n1i
  
  # generate deltai and thetai
  sigma=matrix(c(tau^2,rho*tau,rho*tau,1),2,2)
  m = MASS::mvrnorm(s,c(theta,0),sigma)
  thetai=m[,1]
  deltai=m[,2]
  
  # generate yi
  yi  = runif(s, min =y_min, max = y_max) %>% round()
  yi[which(yi>ni)]=ni[which(yi>ni)]

  # y1i = sapply(1:s, function(i) MCMCpack::rnoncenhypergeom(1, n1i[i], n0i[i], yi[i], exp(thetai[i])))
  y1i = sapply(1:s, function(i) BiasedUrn::rFNCHypergeo(nran=1, m1=n1i[i], m2=n0i[i], n=yi[i], odds=exp(thetai[i])))
  y0i = yi- y1i 

  p.dt = data.frame(y1=y1i,y0=y0i,n1=n1i,n0=n0i,n=ni)
  
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