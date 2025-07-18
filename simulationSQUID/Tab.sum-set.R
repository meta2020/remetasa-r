library(dplyr)
library(kableExtra)

source("Rfn/set.R")

set$t.tau=round(set$t.tau^2,1)

set[,-4]%>%kbl(., 
         format = "html",
         longtable = F, 
         booktabs = T, 
         col.names = c("$\\theta$","$\\tau^2$","$\\rho$",
                       "T:C", "HN$_P$", "BN$_P$", 
                       "NN$_O$", "HN$_O$", "BN$_O$", 
                       "C-N (CP)", "C-S (CP)", 
                       "HN$_{Prop}$ (CP)","BN$_{Prop}$ (CP)"),
         # digits = 1,
         align = "r",
         linesep = c('', '','\\addlinespace'),
         escape = FALSE,
         label = "set1",
         caption = "Simulation: HN") %>% 
  # add_header_above(c("Patients","$\\tau^2$",
  #                    "NN$_O$", "HN$_O$", "BN$_O$", 
  #                    "CopasN (CP)", "CopasH (CP)", 
  #                    "HN$_{SA}$ (CP)","BN$_{SA}$ (CP)"), escape = FALSE) %>% 
  # kable_styling(font_size = 9) %>%
  footnote(general = "
           Patients indicates the median of patients;
           NN$_P$, HN$_P$, and BN$_P$ are the estimates based on the population studies;
           NN$_O$, HN$_O$, and BN$_O$ are the estimates based on the published studies;
           C-N and C-S are the Copas-N and Copas-Shi methods;
           HN$_{prop}$ and BN$_{prop}$ are the proposed sensitivity analysis methods;
           CP indicates the coverage probability.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")

