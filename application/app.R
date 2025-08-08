library(magrittr)
library(kableExtra)

read.csv("egger2001.csv")%>%
  kbl(.,
      format = "latex",
      longtable = F, 
      booktabs = T, 
      align = "r",
      escape = FALSE,
      caption = "Data of Example 1",
      label = "app1")
