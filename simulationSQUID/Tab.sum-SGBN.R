library(dplyr)
library(kableExtra)

## --------------------------------

source("Rfn/set.R")

## create a summary table

sum.tab2 = function(S){
  
rtimes=1000
t.mu=set$t.theta[1]
CP.sum=NULL
theta.sum=NULL
tau.sum=NULL
n.sum=NULL
cv.sum=NULL
nset=nrow(set)/2

for(i in c(1:3,7:9,13:15)){
  # i=2
  load(paste0("res-SGBN/data-set-",i,"-S",S,".RData"))
  DATA0 = DATA %>% t()%>% as.numeric() %>% 
    array(., dim = c(11, 7, rtimes),
          dimnames = list(colnames(DATA),rownames(DATA)[1:7],c(1:rtimes)))
  
  ## remove nonconverged values
  for(j in 1:rtimes){
    DATA0[1:7,,j][,(DATA0[,,j][8,]!=0)]=NA
  }
  # x=do.call(rbind, lapply(1:rtimes, function(i) DATA0[1,,i]))
  # boxplot(x[,])
  # abline(h=-2)
  
  mu    = do.call(rbind, lapply(1:rtimes, function(i) DATA0[1,,i]))
  mu.se = do.call(rbind, lapply(1:rtimes, function(i) DATA0[2,,i]))
  tau   = do.call(rbind, lapply(1:rtimes, function(i) DATA0[5,,i]))
  N     = do.call(rbind, lapply(1:rtimes, function(i) DATA0[11,,i]))
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
PM=theta.sum[,1:2]-t.mu
PM.sp=sprintf("%.1f", PM*100)%>%matrix(,nrow=nset)%>%as.data.frame()
colnames(PM.sp)=colnames(PM)

BIAS=theta.sum[,(3:4)]-t.mu
BIAS.sp=sprintf("%.1f", BIAS*100)%>%matrix(,nrow=nset)%>%as.data.frame()
colnames(BIAS.sp)=colnames(BIAS)

SA=theta.sum[,-(1:4)]-t.mu
SA.CP=100*CP.sum[,-(1:4)]
SA.sp=sprintf("%.1f (%.1f)", SA*100, SA.CP)%>%matrix(,nrow=nset)%>%as.data.frame()
colnames(SA.sp)=colnames(SA)

pt=sprintf("U[%d, %d]",set$nmin,set$nmax)
pt[-seq(1,18,6)]=""

DF = cbind.data.frame(S=c(S, rep("",8)),
                      pt=pt[c(1:3,7:9,13:15)], 
                      tau2=(set$t.tau[c(1:3,7:9,13:15)])^2, 
                      N=round(n.sum[,1],1),
                      POP=PM.sp,PB=BIAS.sp,SA=SA.sp)

DF.cv = cbind.data.frame(S=c(S, rep("",8)),
                         pt=pt[c(1:3,7:9,13:15)], 
                         tau2=(set$t.tau[c(1:3,7:9,13:15)])^2, 
                         N=round(n.sum[,1],1),
                         100-cv.sum*100)

rownames(DF.cv)=NULL

DF.list = list(DF, DF.cv)

return(DF.list)
}


DF.all=rbind.data.frame(sum.tab2(15)[[1]], sum.tab2(50)[[1]])
DF.cv.all=rbind.data.frame(sum.tab2(15)[[2]], sum.tab2(50)[[2]])


DF.cv.all%>%kbl(., 
         format = "latex",
         longtable = F, 
         booktabs = T, 
         col.names = c("S","Patients","$\\tau^2$","N",
                       "NN$_P$", "BN$_P$", 
                       "NN$_O$", "BN$_O$", 
                       "C-N (CP)", "C-H (CP)", 
                       "BN$_{Prop}$ (CP)"),
         # digits = 1,
         align = "r",
         linesep = c('', '','\\addlinespace'),
         escape = FALSE,
         label = "set1",
         position ="!htbp",
         caption = "Simulation: BN") %>% 
  # add_header_above(c("Patients","$\\tau^2$",
  #                    "NN$_O$", "HN$_O$", "BN$_O$", 
  #                    "CopasN (CP)", "CopasH (CP)", 
  #                    "HN$_{SA}$ (CP)","BN$_{SA}$ (CP)"), escape = FALSE) %>% 
  # kable_styling(font_size = 9) %>%
  footnote(general = "
           Patients indicates the median of patients;
           NN$_P$ and BN$_P$ are the estimates based on the population studies;
           NN$_O$ and BN$_O$ are the estimates based on the published studies;
           CopasN and CopasS are the Copas-N and Copas-Shi methods;
           BN$_{prop}$ are the proposed sensitivity analysis methods;
           CP indicates the coverage probability.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")

