## Generate meta-analysis according to the HN-GLMM 
## N are sampled from uniform distributions


gendata.2hn.hedges = function(
    s, 
    n_min, n_max,
    gr,
    theta,tau,rho,
    y_min,y_max,
    Pnmax, Pnmin,
    cutoff=c(0.01,0.05,0.1), ## <0.05, <0.1, others
    wi=c(1,1,0.5,0.3)){
  
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

  p.dt2 = escalc(measure="OR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0, data = p.dt) %>% 
      mutate(ti=yi/sqrt(vi)) %>%
      mutate(pvalue = 2*pnorm(abs(ti), lower.tail = FALSE)) %>%
      mutate(wi=ifelse(pvalue<cutoff[1], wi[1], ifelse(pvalue<cutoff[2], wi[2], ifelse(pvalue<cutoff[3], wi[3], wi[4]))))
    
  z=rbinom(nrow(p.dt),1, p.dt2$wi)
  s.dt=p.dt[z>0,]
  p.dt$pvalue=p.dt2$pvalue
  
  res.list = list(
    p.dt=p.dt,
    s.dt=s.dt,
    Miss=c(M.o=sum(z==0)), 
    p=c(p.e=mean(p.dt2$wi), p.o=mean(z>0))
  )
  return(res.list)
  
}