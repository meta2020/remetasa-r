library(dplyr)
library(kableExtra)

load("scenarios/set.RData")

## create a summary table for 2-sample HN/BN----

sum.event1 = function(S){

  nset=nrow(set)
  
  s.prop.sum=NULL
  n.prop.sum=NULL
  for(i in 1:nset){
    load(paste0("res-2GHN-new/data-set-",i,"-S",S,".RData"))
    rtimes = (dim(DATA)/12)[1]
    DATA0 = DATA %>% t()%>% as.numeric() %>% 
      array(., dim = c(12, 12, rtimes),
            dimnames = list(colnames(DATA),rownames(DATA)[1:12],c(1:rtimes)))
    
    s.prop = do.call(rbind, lapply(1:rtimes, function(i) DATA0[9,,i]))
    n.prop = do.call(rbind, lapply(1:rtimes, function(i) DATA0[10,,i]))
    
    s.prop.sum = rbind(s.prop.sum, colMeans(s.prop, na.rm = T))
    n.prop.sum = rbind(n.prop.sum, colMeans(n.prop, na.rm = T))
  }
  s.prop.sum2=NULL
  n.prop.sum2=NULL
  for(i in 1:nset){
    load(paste0("res-2GBN-new/data-set-",i,"-S",S,".RData"))
    rtimes = (dim(DATA)/12)[1]
    DATA0 = DATA %>% t()%>% as.numeric() %>% 
      array(., dim = c(12, 12, rtimes),
            dimnames = list(colnames(DATA),rownames(DATA)[1:12],c(1:rtimes)))
    
    s.prop = do.call(rbind, lapply(1:rtimes, function(i) DATA0[9,,i]))
    n.prop = do.call(rbind, lapply(1:rtimes, function(i) DATA0[10,,i]))
    
    s.prop.sum2 = rbind(s.prop.sum2, colMeans(s.prop, na.rm = T))
    n.prop.sum2 = rbind(n.prop.sum2, colMeans(n.prop, na.rm = T))
  }
 
  
  pt=sprintf("U[%d, %d]",set$nmin,set$nmax)
  pt[-seq(1,18,6)]=""
  grp=sprintf("%d:1",set$grp)
  grp[-seq(1,18,3)]=""
  
RE.prop=cbind.data.frame(
  s.prop.sum[,1]*100,n.prop.sum[,1]*100,
  s.prop.sum2[,1]*100,n.prop.sum2[,1]*100)

RE.prop.all = cbind.data.frame(S=c(S, rep("",17)),
                               pt, grp, tau2=(set$t.tau)^2,
                               round(RE.prop,1))
return(RE.prop.all) 
}  
  



RE=rbind.data.frame(sum.event1(15), sum.event1(50))

RE%>%kbl(., 
             format = "latex",
             longtable = F, 
             booktabs = T, 
             col.names = c("S","Patients","T:C","$\\tau^2$",
                           "HN$_P$", 
                           "HN$_O$", 
                           "2SBN$_P$", 
                           "2SBN$_O$"),
             align = "r",
             linesep = c('', '','\\addlinespace'),
             escape = FALSE,
             label = "prop1",
             caption = "Proportions of studies with rare events under various scenarios of HN and 2-sample BN model based data generating processes.") %>% 
  footnote(general = "
           Results of HN$_P$ and 2SBN$_P$ are based on the population studies;
           results of HN$_O$ and 2SBN$_O$ are based on the published studies.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")



## create a summary table for 1-sample BN----

sum.event2 = function(S){
  
  nset=6
  s.prop.sum3=NULL
  n.prop.sum3=NULL
  for(i in c(1:3,7:9)){
    load(paste0("res-1GBN-new/data-set-",i,"-S",S,".RData"))
    rtimes = (dim(DATA)/8)[1]
    DATA0 = DATA %>% t()%>% as.numeric() %>% 
      array(., dim = c(12, 8, rtimes),
            dimnames = list(colnames(DATA),rownames(DATA)[1:8],c(1:rtimes)))
    
    s.prop = do.call(rbind, lapply(1:rtimes, function(i) DATA0[9,,i]))
    n.prop = do.call(rbind, lapply(1:rtimes, function(i) DATA0[10,,i]))
    
    s.prop.sum3 = rbind(s.prop.sum3, colMeans(s.prop, na.rm = T))
    n.prop.sum3 = rbind(n.prop.sum3, colMeans(n.prop, na.rm = T))
  }
  
  
  pt=sprintf("U[%d, %d]",round(set$nmin/2),round(set$nmax/2))
  pt[-seq(1,18,6)]=""
  
  RE.prop=cbind.data.frame(
    s.prop.sum3[,1]*100,n.prop.sum3[,1]*100)
  
  RE.prop.all = cbind.data.frame(S=c(S, rep("",5)),
                                 pt[c(1:3,7:9)], tau2=(set$t.tau[c(1:3,7:9)])^2,
                                 round(RE.prop,1))
  return(RE.prop.all) 
}  




RE=rbind.data.frame(sum.event2(15), sum.event2(50))

RE%>%kbl(., 
         format = "latex",
         longtable = F, 
         booktabs = T, 
         col.names = c("S","Patients","$\\tau^2$",
                       "1SBN$_P$", 
                       "1SBN$_O$"),
         position ="!htbp",
         align = "r",
         linesep = c('', '','\\addlinespace'),
         escape = FALSE,
         label = "prop2",
         caption = "Proportions of studies with rare events under various scenarios of 1-sample BN model based data generating processes.") %>% 
  footnote(general = "
           1SBN$_P$ are based on the population studies;
           1SBN$_O$ are based on the published studies.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")

