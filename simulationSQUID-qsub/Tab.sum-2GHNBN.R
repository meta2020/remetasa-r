library(dplyr)
library(kableExtra)

load("scenarios/set.RData")

## create a summary table

sum.tab = function(S){
  
# rtimes=1000  
t.mu=set$t.theta[1]
CP.sum=NULL
theta.sum=NULL
tau.sum=NULL
n.sum=NULL
cv.sum=NULL

nset=nrow(set)
for(i in 1:nset){
  load(paste0("res-2GHN-new/data-set-",i,"-S",S,".RData"))
  rtimes = (dim(DATA)/12)[1]
  DATA0 = DATA %>% t()%>% as.numeric() %>% 
    array(., dim = c(12, 12, rtimes),
          dimnames = list(colnames(DATA),rownames(DATA)[1:12],c(1:rtimes)))
  
  ## remove nonconverged values
  for(j in 1:rtimes){
    DATA0[1:7,,j][,(DATA0[,,j][8,]!=0)]=NA
  }
  
  mu    = do.call(rbind, lapply(1:rtimes, function(i) DATA0[1,,i]))
  mu.se = do.call(rbind, lapply(1:rtimes, function(i) DATA0[2,,i]))
  tau   = do.call(rbind, lapply(1:rtimes, function(i) DATA0[5,,i]))
  N   = do.call(rbind, lapply(1:rtimes, function(i) DATA0[11,,i]))
  CV    = do.call(rbind, lapply(1:rtimes, function(i) DATA0[8,,i]))
  
  CI.lb = mu-1.96*mu.se
  CI.ub = mu+1.96*mu.se
  CP=(CI.lb<t.mu)&(CI.ub>t.mu)
  CP.sum=rbind(CP.sum, colMeans(CP, na.rm = T))
  
  theta.sum = rbind(theta.sum, colMeans(mu, na.rm = T))
  tau.sum   = rbind(tau.sum, colMeans(tau, na.rm = T))
  n.sum     = rbind(n.sum, colMeans(N, na.rm = T))
  cv.sum    = rbind(cv.sum, colMeans(CV, na.rm = T))
}
## 1:1 set with biased and adjusted
PM=theta.sum[,1:3]-t.mu
PM.sp=sprintf("%.1f", PM*100)%>%matrix(,nrow=nset)%>%as.data.frame()
colnames(PM.sp)=colnames(PM)

BIAS=theta.sum[,(4:6)]-t.mu
BIAS.sp=sprintf("%.1f", BIAS*100)%>%matrix(,nrow=nset)%>%as.data.frame()
colnames(BIAS.sp)=colnames(BIAS)


SA=theta.sum[,-(1:6)]-t.mu
SA.CP=100*CP.sum[,-(1:6)]
SA.sp=sprintf("%.1f (%.1f)", SA*100, SA.CP)%>%matrix(,nrow=nset)%>%as.data.frame()
colnames(SA.sp)=colnames(SA)


pt=sprintf("U[%d, %d]",set$nmin,set$nmax)
pt[-seq(1,18,6)]=""
grp=sprintf("%d:1",set$grp)
grp[-seq(1,18,3)]=""

DF = cbind.data.frame(S=c(S, rep("",17)),
                      pt, grp, tau2=(set$t.tau)^2, N=round(n.sum[,1],1),
                      POP=PM.sp,PB=BIAS.sp,SA=SA.sp)

DF.cv = cbind.data.frame(S=c(S, rep("",17)),
                         pt, grp, tau2=(set$t.tau)^2, N=round(n.sum[,1],1), 
                         100-cv.sum*100)
rownames(DF.cv)=NULL

tPM=tau.sum[,1:3]
tPM.sp=sprintf("%.2f", tPM)%>%matrix(,nrow=nset)%>%as.data.frame()
colnames(tPM.sp)=colnames(tPM)

tBIAS=tau.sum[,(4:6)]
tBIAS.sp=sprintf("%.2f", tBIAS)%>%matrix(,nrow=nset)%>%as.data.frame()
colnames(tBIAS.sp)=colnames(tBIAS)


tSA=tau.sum[,-(1:6)]
tSA.sp=sprintf("%.2f", tSA)%>%matrix(,nrow=nset)%>%as.data.frame()
colnames(tSA.sp)=colnames(tSA)


tDF = cbind.data.frame(S=c(S, rep("",17)),
                      pt, grp, tau2=(set$t.tau)^2,N=round(n.sum[,1],1), 
                      POP=tPM.sp,PB=tBIAS.sp,SA=tSA.sp)

DF.list = list(DF, DF.cv, tDF)

return(DF.list)
}


##S=50----------

DF.all=rbind.data.frame(sum.tab(15)[[1]], sum.tab(50)[[1]])
DF.cv.all=rbind.data.frame(sum.tab(15)[[2]], sum.tab(50)[[2]])
tDF.all=rbind.data.frame(sum.tab(15)[[3]], sum.tab(50)[[3]])



tDF.all%>%kbl(., 
         format = "html",
         longtable = F, 
         booktabs = T, 
         col.names = c("$S$","Patients","T:C","$\\tau^2$","$N$",
                       "NN$_P$", "HN$_P$", "CBN$_P$", 
                       "NN$_O$", "HN$_O$", "CBN$_O$", 
                       "CN (CP)", "CS (CP)", 
                       "HN$^\\text{Prop}$ (CP)","CBN$^\\text{Prop}$ (CP)",
                       "HN$^\\text{HTJ}$ (CP)","CBN$^\\text{HTJ}$ (CP)"),
         align = "r",
         linesep = c('', '','\\addlinespace'),
         escape = FALSE,
         label = "set1",
         position ="!htbp",
         caption = "Simulation1") %>% 
  footnote(general = "
           Results of NN$_P$, HN$_P$, and CBN$_P$ are based on the population studies;
           results of NN$_O$, HN$_O$, and CBN$_O$ are based on the published studies;
           CN and CS indicate the Copas-N and Copas-Shi methods;
           HN$^\\text{Prop}$ and CBN$^\\text{Prop}$ indicate the proposed HN or CBN model based sensitivity analysis methods;
           CP indicates the coverage probability.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")

