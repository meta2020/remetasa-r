## Generate meta-analysis according to the HN-GLMM 
## N are sampled from uniform distributions


gendata.2hn.hedges = function(
    s, 
    n_min, n_max,
    gr,
    theta,tau,rho,
    y_min,y_max,
    Pnmax, Pnmin,
    cutoff=c(0.05,0.1), ## <0.05, <0.1, others
    wi=c(0.5,0.3)){
  
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

  p.dt = data.frame(y1=y1i,y0=y0i,n1=n1i,n0=n0i,n=ni)%>%
    mutate( cx = ( ( y1 == 0 ) | (y1 ==n1 ) | (y0 == 0 ) | (y0 == n0) )*0.5 ) %>%## selective process
    mutate( y = log( ( y1 + cx )*( n0-y0 + cx )/( n1-y1 + cx )/( y0 + cx ) ), 
            v = 1/( y1 + cx ) + 1/( n0-y0 + cx ) + 
              1/( n1-y1 + cx ) + 1/( y0 + cx )) %>%
    mutate( t = y/sqrt(v) ) %>% 
    mutate(pvalue = pnorm(abs(t), lower.tail = FALSE))%>%
    mutate(wi=ifelse(pvalue<cutoff[1], 1, ifelse(pvalue<cutoff[2], wi[1], wi[2])))
  
  z=rbinom(nrow(p.dt),1, p.dt$wi)
  
  s.dt=p.dt[z>0,]
  
  res.list = list(
    p.dt=p.dt,
    s.dt=s.dt,
    Miss=c(M.o=sum(z==0)), 
    p=c(p.e=mean(p.dt$wi), p.o=mean(z>0))
  )
  return(res.list)
  
}