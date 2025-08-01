library(dplyr)
library(kableExtra)

source("Rfn/set.R")

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
         format = "latex",
         longtable = F, 
         booktabs = T, 
         col.names = c("$\\theta$","$\\rho$","Total event ($y_i$)","$(P_{max}, P_{min})$",
                       "Total subjects ($n_i$)", "$p_0$","T:C","$\\tau^2$" 
                       ),
         # digits = 1,
         align = "r",
         linesep = c('', '','\\addlinespace'),
         escape = FALSE,
         label = "set1",
         caption = "tab:set") %>% 
  footnote(general = "
           T:C indicates Treatment:Control.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")

