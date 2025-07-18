library(dplyr)
library(kableExtra)

## --------------------------------

source("Rfn/set.R")
rtimes=1000
## create a summary table
t.mu=set$t.theta[1]
S=15
CP.sum=NULL
bias.sum=NULL
sum.cv.all=NULL
nset=nrow(set)
for(i in 1:nrow(set)){
  # i=2
  load(paste0("res-all-BN2/data-set-",i,"-S",S,".RData"))
  DATA0 = DATA %>% t()%>% as.numeric() %>% 
    array(., dim = c(11, 7, rtimes),
          dimnames = list(colnames(DATA),rownames(DATA)[1:7],c(1:rtimes)))
  
  ## remove nonconverged values
  for(j in 1:rtimes){
    DATA0[1:7,,j][,(DATA0[,,j][8,]!=0)]=NA
  }
  x=do.call(rbind, lapply(1:rtimes, function(i) DATA0[1,,i]))
  boxplot(x[,])
  abline(h=-2)
  
  mu = do.call(rbind, lapply(1:rtimes, function(i) DATA0[1,,i]))
  mu.se = do.call(rbind, lapply(1:rtimes, function(i) DATA0[2,,i]))
  CI.lb = mu-1.96*mu.se
  CI.ub = mu+1.96*mu.se
  CP=(CI.lb<t.mu)&(CI.ub>t.mu)
  CP.sum=rbind(CP.sum, colMeans(CP, na.rm = T))
  bias.sum=rbind(bias.sum, colMeans(mu, na.rm = T))
  mean_matrix = apply(DATA0, c(1, 2), function(x) mean(x, na.rm=T))[8,]
  sum.cv.all=rbind(sum.cv.all, mean_matrix)
}
## 1:1 set with biased and adjusted
PM=bias.sum[,1:2]-t.mu
PM.sp=sprintf("%.3f", PM)%>%matrix(,nrow=nset)%>%as.data.frame()
colnames(PM.sp)=colnames(PM)

BIAS=bias.sum[,(3:4)]-t.mu
BIAS.sp=sprintf("%.3f", BIAS)%>%matrix(,nrow=nset)%>%as.data.frame()
colnames(BIAS.sp)=colnames(BIAS)

SA=bias.sum[,-(1:4)]-t.mu
SA.CP=100*CP.sum[,-(1:6)]
SA.sp=sprintf("%.3f (%.2f)", SA, SA.CP)%>%matrix(,nrow=nset)%>%as.data.frame()
colnames(SA.sp)=colnames(SA)


DF = cbind.data.frame(n.med=set$n.median, grp=set$grp.r, tau2=(set$t.tau)^2,
                      POP=PM.sp,PB=BIAS.sp,SA=SA.sp)

DF.cv = cbind.data.frame(n.med=set$n.median, grp=set$grp.r, tau2=(set$t.tau)^2,
                      sum.cv.all*100)
rownames(DF.cv)=NULL


DF%>%kbl(., 
         format = "html",
         longtable = F, 
         booktabs = T, 
         col.names = c("Patients","T:C","$\\tau^2$",
                       "NN$_P$", "BN$_P$", 
                       "NN$_O$", "BN$_O$", 
                       "CopasN (CP)", "CopasH (CP)", 
                       "BN$_{Prop}$ (CP)"),
         # digits = 1,
         align = "r",
         linesep = c('', '','\\addlinespace'),
         escape = FALSE,
         label = "set1",
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

