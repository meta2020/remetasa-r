library(dplyr)
library(kableExtra)

## --------------------------------

source("Rfn/set.R")

## create a summary table

sum.tab2 = function(S){
  
rtimes=1000
tau.sum=NULL

nset=nrow(set)/2
for(i in c(1:3,7:9,13:15)){
  # i=2
  load(paste0("res-BNprop/data-set-",i,"-S",S,".RData"))
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
  tau   = do.call(rbind, lapply(1:rtimes, function(i) DATA0[5,,i]))
  tau.sum   = rbind(tau.sum, colMeans(tau, na.rm = T))
}
## 1:1 set with biased and adjusted
PM=tau.sum[,1:2]
PM.sp=sprintf("%.2f", PM)%>%matrix(,nrow=nset)%>%as.data.frame()
colnames(PM.sp)=colnames(PM)

BIAS=tau.sum[,(3:4)]
BIAS.sp=sprintf("%.2f", BIAS)%>%matrix(,nrow=nset)%>%as.data.frame()
colnames(BIAS.sp)=colnames(BIAS)

SA=tau.sum[,-(1:4)]
SA.sp=sprintf("%.2f", SA)%>%matrix(,nrow=nset)%>%as.data.frame()
colnames(SA.sp)=colnames(SA)

pt=sprintf("U[%d, %d]",set$nmin,set$nmax)
pt[-seq(1,18,6)]=""
grp=sprintf("%d:1",set$grp)
grp[-seq(1,18,3)]=""

DF = cbind.data.frame(S=c(S, rep("",8)),
                      pt=pt[c(1:3,7:9,13:15)], 
                      tau2=(set$t.tau[c(1:3,7:9,13:15)])^2, 
                      N=round(n.sum[,1],1),
                      POP=PM.sp,PB=BIAS.sp,SA=SA.sp)

DF.list = list(DF)

return(DF.list)
}


DF.all=rbind.data.frame(sum.tab2(15)[[1]], sum.tab2(50)[[1]])


DF.all%>%kbl(., 
         format = "latex",
         longtable = F, 
         booktabs = T, 
         col.names = c("S","Patients","T:C","$\\tau^2$",
                       "NN$_P$", "BN$_P$", 
                       "NN$_O$", "BN$_O$", 
                       "CN", "CH", 
                       "BN$_{Prop}$"),
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
           CN and CS are the Copas-N and Copas-Shi methods;
           BN$_{prop}$ are the proposed sensitivity analysis methods;
           CP indicates the coverage probability.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")

