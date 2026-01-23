## Generate meta-analysis according to the BN-GLMM
## N are sampled from uniform distributions

gendata.1bn.hedges = function(
    s, 
    n_min, n_max,
    theta,tau,rho,
    Pnmax, Pnmin,
    cutoff=c(0.05,0.1), ## <0.05, <0.1, others
    wi=c(0.5,0.3)){
  
  ni = runif( s, min = n_min, max = n_max )%>%round()
  ni = ifelse(ni<10,10,ni)
  
  # generate deltai and thetai
  sigma=matrix(c(tau^2,rho*tau,rho*tau,1),2,2)
  m = MASS::mvrnorm(s,c(theta,0),sigma)
  thetai=m[,1]
  deltai=m[,2]
  
  # generate yi
  yi = sapply(1:s, function(i) rbinom(1, ni[i], plogis(thetai[i])) )
  
  p.dt = data.frame(y=yi,n=ni)%>%
    mutate( cx = ( ( y == 0 ) | (y ==n ) )*0.5 ) %>%## selective process
    mutate( prop = plogis( ( y + cx )/( n-y + cx )), 
            v = 1/( y + cx ) + 1/( n-y + cx )) %>%
    mutate( t = prop/sqrt(v) ) %>% 
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