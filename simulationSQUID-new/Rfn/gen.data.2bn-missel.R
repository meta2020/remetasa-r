## Gnerate meta-analysis by group according to the BN model
##
## From Ao
##
gendata.2bn.hedges = function(
    s, 
    n_min, n_max, p0,
    gr,
    theta,tau,rho,
    Pnmax, Pnmin){
  
    n = runif( s, min = n_min, max = n_max )
    n = ifelse(n<20,20,n)
    
    n0i = (n*(1 / (1 + gr))) %>% round()
    n1i = (n-n0i)%>% round() # #treatment subjects
    ni=n0i+n1i
    
    ## thetai=log(ORi) and deltai 
    sigma=matrix(c(tau^2,rho*tau,rho*tau,1),2,2)
    m = MASS::mvrnorm(s,c(theta,0),sigma)
    thetai=m[,1]
    deltai=m[,2]
  
    ORi = exp(thetai) # empirical odds ratio
    # pi0 = runif(s,0.2,0.9)
    pi0 = plogis(rnorm(s, qlogis(p0), tau / 2)) # prob of control in study i
    m = ORi * pi0 / ( 1 - pi0 )
    pi1 = m / ( 1 + m )   # prob of event in study i

    all.dat = vapply(1:s,function(i){
      yi1 = sum(rbinom(n1i[i],1,pi1[i]))
      ni1 = n1i[i]
      yi0 = sum(rbinom(n0i[i],1,pi0[i]))
      ni0 = n0i[i]
      c(yi1,ni1,yi0,ni0)
      },c(y1=0,n1=0,y0=0,n0=0))%>%t()%>%data.frame()
    
    p.dt = data.frame(y1=all.dat$y1,y0=all.dat$y0,
                      n1=all.dat$n1,n0=all.dat$n0,n=ni)
    
    ## selective process
   
    
    
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

