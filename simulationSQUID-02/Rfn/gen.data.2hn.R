## Generate meta-analysis according to the HN-GLMM 
## N are sampled from uniform distributions
##
## Generate population data and selected data (correct selection)
##

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

##
## when selection function is misspecified
##
gendata.2hn.hedges = function(
    s, 
    n_min, n_max,
    gr,
    theta,tau,rho,
    y_min,y_max,
    cutoff, ## <0.05, <0.1, others
    wi){
  
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
      mutate(wi=ifelse(pvalue<cutoff[1], wi[1], 
                      ifelse(pvalue<cutoff[2], wi[2], 
                             ifelse(pvalue<cutoff[3], wi[3],
                                    ifelse(pvalue<cutoff[4], wi[4], wi[5])))))
    
  z=rbinom(nrow(p.dt),1, p.dt2$wi)
  s.dt=p.dt2[z>0,]
  
  res.list = list(
    p.dt=p.dt2,
    s.dt=s.dt,
    Miss=c(M.o=sum(z==0)), 
    p=c(p.e=mean(p.dt2$wi), p.o=mean(z>0))
  )
  return(res.list)
  
}


gendata.2hn.se = function(
    s, 
    n_min, n_max,
    gr,
    theta,tau,rho,
    y_min,y_max,
    Psemax, Psemin){
  
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
  p.dt2 = escalc(measure="OR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0, data=p.dt)%>% 
      mutate(preci=1/sqrt(vi))
  preci = p.dt2$preci

    ## selective process
    semin=min(preci) 
    semax=max(preci)
    
    a1=(qnorm(Psemax)-qnorm(Psemin))/((semax)-(semin))
    a0=qnorm(Psemax)-a1*(semax)
    
    zi=a0+a1*preci+deltai
    wi=pnorm(a0+a1*preci)

    wi.s = wi[zi>0]
    M= sum((1-wi.s)/wi.s)%>%round()
    p.dt$zi=zi
    p.dt$wi=wi
    
    s.dt=p.dt[zi>0,]
    
    res.list = list(
      p.dt=p.dt,
      s.dt=s.dt,
      a=c(a1=a1,a0=a0),
      Miss=c(M.e=M, M.o=sum(zi<=0)), 
      p=c(p.e=mean(wi), p.o=mean(zi>0))
    )
    return(res.list)   
  
}


gendata.2hn.t = function(
    s, 
    n_min, n_max,
    gr,
    theta,tau,rho,
    y_min,y_max,
    alpha, beta){
  
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
  p.dt2 = escalc(measure="OR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0, data=p.dt)%>% 
      mutate(ti = yi/sqrt(vi)) %>%
      mutate(zi = alpha+beta*ti+deltai)%>%
      mutate(wi = pnorm(alpha+beta*ti))
  
  s.dt=p.dt2[p.dt2$zi>0,]
  
  res.list = list(
    p.dt=p.dt2,
    s.dt=s.dt
  )
    return(res.list)   
  
}