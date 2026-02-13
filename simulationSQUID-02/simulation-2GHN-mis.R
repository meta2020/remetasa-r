##
## Compare all methods based on HN model data generating process (Add Hu et al.)
##
log.name = paste0("log-mis1-2GHN",as.numeric(Sys.time()),".txt")
sink(log.name)
msg_file = file(log.name, open="at")
sink(msg_file, type = "message")


rm(list=ls())

## load functions and settings
file.sources = list.files("Rfn/")
sapply(paste0("Rfn/", file.sources), source)

## load scenarios
load("scenarios/set.RData")


## settings for simulations
rtimes=1000

set.eps = 1e-4
set.int.lmt = 10
set.cub.tol = 1e-4

set.cutoff=c(0.01,0.05,0.1)
set.wi=c(0.99, 0.99, 0.4, 0.3)

set.tau.bound = 2
set.mu.bound = abs(-2)*3

## SIMULATION  ----------------------------------
ncores = min(120, parallel::detectCores())
cl = parallel::makeCluster(ncores, "SOCK")
doSNOW::registerDoSNOW(cl)

message(paste0("Start",Sys.time()))


set.seed(2025)
for(S in s[1]){
for(i in 1:nrow(set)){ 

##-- Simulation 1: HN model based -------


  DATA = foreach(r=1:rtimes, .combine = rbind,.packages=c("mnormt","dplyr","metafor"), .errorhandling = "remove")  %dorng%  {

    set.gr = set[i,]
    # pmax = set.gr$pmax
    # pmin = set.gr$pmin
    
    
    
    ## create datasets
    plist = gendata.2hn.hedges(
      s=S, 
      n_min = set.gr$nmin, n_max=set.gr$nmax,
      gr=set.gr$grp.r,
      theta=set.gr$t.theta,
      tau=set.gr$t.tau,
      rho=set.gr$t.rho,
      y_min=set.gr$ymin, y_max=set.gr$ymax,
      Pnmax = pmax, Pnmin = pmin,
      cutoff=set.cutoff, ## <0.05, <0.1, others
      wi=set.wi) ## weight=1,0.5,0.3

    
    ## selected data and population data
    sdata = plist$s.dt
    pdata = plist$p.dt

    ## data with n and lnOR
    data = lapply(list(pdata, sdata),
                  function(data) escalc(measure="OR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0, data=data)
    )
    lpdata = data[[1]]
    lsdata = data[[2]]

    
    ## proportion of rare events
    p.small = lpdata %>%summarise(prop = mean(y1 <=3 | y0 <=3))%>%c()
    s.small = lsdata %>%summarise(prop = mean(y1 <=3 | y0 <=3))%>%c()

    ## population and selected models 
    
    ## set parset list
    parset.nn = list(
      mu.bound = set.mu.bound, 
      tau.bound = set.tau.bound,
      eps = set.eps,
      init.vals = c(set.gr$t.theta+round(runif(1,-0.2,0.2),2),
                    set.gr$t.tau  +round(runif(1,-0.2,0.2),2))
    )
    
    ## estimation without/with PB: NN, HN-GLMM, BN-GLMM on pdata and sdata
    fit.nn = lapply(
      list(lpdata,lsdata), 
      function(data) with(data, NN_LMM(yi, vi, parset=parset.nn)))
    
    pnn = c(fit.nn[[1]]$mu, fit.nn[[1]]$tau, rho=NA, rho.se=NA,
             cv = ifelse(is.null(fit.nn[[1]]$opt$convergence), NA, fit.nn[[1]]$opt$convergence))
    snn = c(fit.nn[[2]]$mu, fit.nn[[2]]$tau, rep(NA,2),
             cv = ifelse(is.null(fit.nn[[2]]$opt$convergence), NA, fit.nn[[2]]$opt$convergence))
    
    ## initial values for GLMM
    parset.glmm = list(
      mu.bound = set.mu.bound,
      tau.bound = set.tau.bound,
      eps = set.eps,
      integ.limit = set.int.lmt, 
      cub.tol = set.cub.tol,
      init.vals = c(set.gr$t.theta+round(runif(1,-0.2,0.2),2),
                    set.gr$t.tau  +round(runif(1,-0.2,0.2),2))
      )
    
    fit.hn = lapply(
      list(lpdata, lsdata), 
      function(data) with(data, HN_GLMM(y0, y1, n0, n1, parset = parset.glmm)))
    phn = c(fit.hn[[1]]$mu, fit.hn[[1]]$tau, rep(NA,2),
             cv = ifelse(is.null(fit.hn[[1]]$opt$convergence), NA, fit.hn[[1]]$opt$convergence))
    shn = c(fit.hn[[2]]$mu, fit.hn[[2]]$tau, rep(NA,2),
             cv = ifelse(is.null(fit.hn[[2]]$opt$convergence), NA, fit.hn[[2]]$opt$convergence))
    
    fit.bn = lapply(
      list(lpdata, lsdata), 
      function(data) with(data, BN_GLMM(y0, y1, n0, n1, parset = parset.glmm)))
    pbn = c(fit.bn[[1]]$mu, fit.bn[[1]]$tau, rep(NA,2),
             cv = ifelse(is.null(fit.bn[[1]]$opt$convergence), NA, fit.bn[[1]]$opt$convergence))
    sbn = c(fit.bn[[2]]$mu, fit.bn[[2]]$tau, rep(NA,2),
             cv = ifelse(is.null(fit.bn[[2]]$opt$convergence), NA, fit.bn[[2]]$opt$convergence))
   
    
    ## adjusted models 
    ## Copas methods initial values
    parset.copas = list(
    mu.bound = set.mu.bound,
    tau.bound = set.tau.bound,
    estimate.rho = TRUE, 
    eps = set.eps,
    init.vals = c(set.gr$t.theta+round(runif(1,-0.2,0.2),2),
                  set.gr$t.tau  +round(runif(1,-0.2,0.2),2), 
                  set.gr$t.rho  +round(runif(1,-0.2,0.1),2)) ## initials for mu tau and rho
    )
    ## proposed method initial values
    parset.new = list(
      mu.bound = set.mu.bound,
      tau.bound = set.tau.bound,
      estimate.rho = TRUE, 
      eps = set.eps,
      integ.limit = set.int.lmt, 
      cub.tol = set.cub.tol,
      init.vals = c(set.gr$t.theta+round(runif(1,-0.2,0.2),2),
                    set.gr$t.tau  +round(runif(1,-0.2,0.2),2), 
                    set.gr$t.rho  +round(runif(1,-0.2,0.1),2)) ## initials for mu tau and rho
    )
    
    ## HTJ method initial values
    parset.htj = list(
      mu.bound = set.mu.bound,
      tau.bound = set.tau.bound,
      beta.bound = 5,
      alpha.bound = 10,
      eps = set.eps,
      integ.limit = set.int.lmt, 
      cub.tol = set.cub.tol,
      init.vals = c(set.gr$t.theta+round(runif(1,-0.2,0.2),2),
                    set.gr$t.tau  +round(runif(1,-0.2,0.2),2), 
                    round(runif(1,0.5,1),2)) ## initial value for mu tau and beta
    )
    
    nmin = min(lsdata$n)
    nmax = max(lsdata$n)
    pmax = lsdata$wi[which.max(lsdata$n)]
    pmin = lsdata$wi[which.min(lsdata$n)]
    p = nrow(sdata)/nrow(pdata)
    
    if(p==1 || nrow(sdata)==0) next else {
      
      func_list = list(
        function(x) with(x, 
                         COPAS1999(yi, n, Pnmax = pmax, Pnmin = pmin, parset=parset.copas)),
        function(x) with(x, 
                         COPAS2000(yi, vi, Psemax = pmax, Psemin = pmin, parset=parset.copas)),
        function(x) with(x,
                         COPAS_HNGLMM(y0, y1, n0, n1, Pnmax = pmax, Pnmin = pmin,
                                      n_min = nmin, n_max = nmax, parset=parset.new)),
        function(x) with(x,
                         COPAS_BNGLMM(y0, y1, n0, n1, Pnmax = pmax, Pnmin = pmin,
                                      n_min = nmin, n_max = nmax, parset=parset.new)),
        function(x) with(x,
                         HTJ_HNGLMM(y0, y1, n0, n1, p=p, parset=parset.htj)),
        function(x) with(x,
                         HTJ_BNGLMM(y0, y1, n0, n1, p=p, parset=parset.htj))
      )
      
      adj.list = suppressWarnings(lapply(func_list, function(f) f(lsdata)))
      
      copas99 = c(adj.list[[1]]$mu, adj.list[[1]]$tau,adj.list[[1]]$rho,
                  cv = ifelse(is.null(adj.list[[1]]$opt$convergence), NA, adj.list[[1]]$opt$convergence))
      copas20 = c(adj.list[[2]]$mu, adj.list[[2]]$tau,adj.list[[2]]$rho,
                  cv = ifelse(is.null(adj.list[[2]]$opt$convergence), NA, adj.list[[2]]$opt$convergence))
      adjhn   = c(adj.list[[3]]$mu, adj.list[[3]]$tau,adj.list[[3]]$rho,
                  cv = ifelse(is.null(adj.list[[3]]$opt$convergence), NA, adj.list[[3]]$opt$convergence))
      adjbn   = c(adj.list[[4]]$mu, adj.list[[4]]$tau,adj.list[[4]]$rho,
                  cv = ifelse(is.null(adj.list[[4]]$opt$convergence), NA, adj.list[[4]]$opt$convergence))
      htjhn   = c(adj.list[[5]]$mu, adj.list[[5]]$tau, NA, NA,
                  cv = ifelse(is.null(adj.list[[5]]$opt$convergence), NA, adj.list[[5]]$opt$convergence))
      htjbn   = c(adj.list[[6]]$mu, adj.list[[6]]$tau, NA, NA,
                  cv = ifelse(is.null(adj.list[[6]]$opt$convergence), NA, adj.list[[6]]$opt$convergence))
      
      
      res.est = rbind(pnn,phn,pbn,
                      snn,shn,sbn,
                      copas99,copas20,
                      adjhn,adjbn,
                      htjhn, htjbn
                      )
      
    }
    
    res = cbind(res.est,
                p.prop=c(p.small, rep(NA,11)),
                s.prop=c(s.small, rep(NA,11)),
                n.pub=c(nrow(sdata), rep(NA,11)),
                p.pub=c(p, rep(NA,11)))

    res

  }
  save(DATA,file = paste0("res-2GHN-mis1/data-set-",i,"-S",S,".RData"))
  message(paste0("Finish-2GHN-mis1-",S,"-",i,":",Sys.time()))

}}


parallel::stopCluster(cl)
message(paste0("Finish",Sys.time()))

sink(type = "message")
# close(msg_file)
sink()
