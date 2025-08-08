library(magrittr)
library(kableExtra)

read.csv("niel-weise21.csv")%>%
  kbl(.,
      format = "latex",
      longtable = F, 
      booktabs = T, 
      digits = 3,
      align = "r",
      linesep = c(rep("",9), "\\addlinespace"),
      escape = FALSE,
      caption = "Summary of the estimations of different sensitivity analysis methods",
      label = "tab1"))
