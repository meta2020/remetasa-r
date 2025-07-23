library(dplyr)
library(kableExtra)

source("Rfn/set.R")

## create a summary table

sum.tab3 = function(S){
  
  rtimes=1000  

  nset=nrow(set)
  
  s.prop.sum=NULL
  n.prop.sum=NULL
  for(i in 1:nset){
    load(paste0("res-HN1/data-set-",i,"-S",S,".RData"))
    DATA0 = DATA %>% t()%>% as.numeric() %>% 
      array(., dim = c(11, 10, rtimes),
            dimnames = list(colnames(DATA),rownames(DATA)[1:10],c(1:rtimes)))
    
    s.prop = do.call(rbind, lapply(1:rtimes, function(i) DATA0[9,,i]))
    n.prop = do.call(rbind, lapply(1:rtimes, function(i) DATA0[10,,i]))
    
    s.prop.sum = rbind(s.prop.sum, colMeans(s.prop, na.rm = T))
    n.prop.sum = rbind(n.prop.sum, colMeans(n.prop, na.rm = T))
  }
  s.prop.sum2=NULL
  n.prop.sum2=NULL
  for(i in 1:nset){
    load(paste0("res-2GBN/data-set-",i,"-S",S,".RData"))
    DATA0 = DATA %>% t()%>% as.numeric() %>% 
      array(., dim = c(11, 10, rtimes),
            dimnames = list(colnames(DATA),rownames(DATA)[1:10],c(1:rtimes)))
    
    s.prop = do.call(rbind, lapply(1:rtimes, function(i) DATA0[9,,i]))
    n.prop = do.call(rbind, lapply(1:rtimes, function(i) DATA0[10,,i]))
    
    s.prop.sum2 = rbind(s.prop.sum2, colMeans(s.prop, na.rm = T))
    n.prop.sum2 = rbind(n.prop.sum2, colMeans(n.prop, na.rm = T))
  }
  s.prop.sum3=NULL
  n.prop.sum3=NULL
  for(i in 1:nset){
    load(paste0("res-BNprop/data-set-",i,"-S",S,".RData"))
    DATA0 = DATA %>% t()%>% as.numeric() %>% 
      array(., dim = c(11, 7, rtimes),
            dimnames = list(colnames(DATA),rownames(DATA)[1:7],c(1:rtimes)))
    
    s.prop = do.call(rbind, lapply(1:rtimes, function(i) DATA0[9,,i]))
    n.prop = do.call(rbind, lapply(1:rtimes, function(i) DATA0[10,,i]))
    
    s.prop.sum3 = rbind(s.prop.sum3, colMeans(s.prop, na.rm = T))
    n.prop.sum3 = rbind(n.prop.sum3, colMeans(n.prop, na.rm = T))
  }
  
  
  pt=sprintf("U[%d, %d]",set$nmin,set$nmax)
  pt[-seq(1,18,6)]=""
  grp=sprintf("%d:1",set$grp)
  grp[-seq(1,18,3)]=""
  
RE.prop=cbind.data.frame(
  s.prop.sum[,1]*100,s.prop.sum2[,1]*100,s.prop.sum3[,1]*100,
  n.prop.sum[,1]*100,n.prop.sum2[,1]*100,n.prop.sum3[,1]*100)

RE.prop.all = cbind.data.frame(S=c(15, rep("",17),50, rep("",17)),
                               pt, grp, tau2=(set$t.tau)^2,
                               round(RE.prop,1))
return(RE.prop.all) 
}  
  



RE=rbind.data.frame(sum.tab3(15), sum.tab3(50))



RE%>%kbl(., 
             format = "html",
             longtable = F, 
             booktabs = T, 
             col.names = c("S","Patients","T:C","$\\tau^2$",
                           "HN$_P$", "Real-Word$_P$", "SGBN$_P$", 
                           "HN$_O$", "Real-Word$_O$", "SGBN$_O$"),
             # digits = 1,
             align = "r",
             linesep = c('', '','\\addlinespace'),
             escape = FALSE,
             label = "set1",
             caption = "RE") %>% 
  footnote(general = "
           NN$_P$, HN$_P$, and BN$_P$ are the estimates based on the population studies;
           NN$_O$, HN$_O$, and BN$_O$ are the estimates based on the published studies.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")

