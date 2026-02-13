library(dplyr)
library(kableExtra)

load("scenarios/set.RData")

set$t.tau=round(set$t.tau^2,1)
set$grp.r=sprintf("%d:1", set$grp.r)
set$pmax=sprintf("(%.2f, %.2f)", set$pmax, set$pmin)
set$ymin=sprintf("[%d, %d]", set$ymin, set$ymax)
set$nmin=sprintf("[%d, %d]", set$nmin, set$nmax)
set=set[,c("t.theta","t.rho","ymin","pmax","nmin","p0","grp.r","t.tau")]
set[-1,1:4]=""
set[-seq(1,18,6),5:6]=""
set[-seq(1,18,3),7]=""


set%>%kbl(., 
         format = "html",
         longtable = F, 
         booktabs = T, 
         col.names = c("$\\theta$","$\\rho$","Total event ($y_i$)","$(P_{max}, P_{min})$",
                       "Total subjects ($n_i$)", "$p_0$","T:C","$\\tau^2$" 
                       ),
         position ="!htbp",
         align = "r",
         linesep = c('', '','\\addlinespace'),
         escape = FALSE,
         label = "set",
         caption = "Scenarios for simulating meta-analysis of odds ratios") %>% 
  footnote(general = "
           T:C indicates Treatment:Control.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")


load("scenarios/set.RData")
set$t.tau=round(set$t.tau^2,1)
set$grp.r=sprintf("%d:1", set$grp.r)
set$pmax=sprintf("(%.2f, %.2f)", set$pmax, set$pmin)
set$nmin=sprintf("[%d, %d]", round(set$nmin/2), round(set$nmax/2))
set=set[,c("t.theta","t.rho","pmax","nmin","t.tau")]
set1=set[c(1:3,7:9),]
set1[-1,1:3]=""
set1[-c(1,4),4]=""
rownames(set1)=NULL
set1%>%kbl(., 
          format = "latex",
          longtable = F, 
          booktabs = T, 
          col.names = c("$\\theta$","$\\rho$","$(P_{max}, P_{min})$",
                        "Total subjects ($n_i$)", "$\\tau^2$" 
          ),
          position ="!htbp",
          align = "r",
          linesep = c('', '','\\addlinespace'),
          escape = FALSE,
          label = "set",
          caption = "Scenarios for simulating meta-analysis of proportions")
