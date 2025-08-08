library(magrittr)
library(kableExtra)

read.csv("pritz1997.csv")[,-2]%>%
  kbl(.,
      format = "latex",
      longtable = F, 
      booktabs = T, 
      align = "r",
      escape = FALSE,
      caption = "Data of Example 1",
      label = "app1")
