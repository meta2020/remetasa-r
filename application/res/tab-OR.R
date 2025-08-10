rm(list = ls())
sink("output.txt") 
library(magrittr)

## TABLE
load("app1-nopb.RData")
load("app1-OR.RData")
load("app1-t-OR.RData")


tab1 = data.frame(
  p = sprintf("(%.2f, %.3f)", tab1_all$pnmin, tab1_all$pnmax),
  Mp = tab1_all$M.p,
  HN = sprintf("%.3f (%.3f, %.3f)", tab1_all$HN.mu, tab1_all$HN.mu.lb, tab1_all$HN.mu.ub),
  HN.tau = sprintf("%.3f", lnOR_COPAS_HNGLMM[,5]),
  HN.rho = sprintf("%.3f", lnOR_COPAS_HNGLMM[,7]),
  BN = sprintf("%.3f (%.3f, %.3f)", tab1_all$BN.mu, tab1_all$BN.mu.lb, tab1_all$BN.mu.ub),
  BN.tau = sprintf("%.3f", lnOR_COPAS_BNGLMM[,5]),
  BN.rho = sprintf("%.3f", lnOR_COPAS_BNGLMM[,7])
)
colnames(tab1) = c("$(\\Pmin, \\Pmax)$", "$M$", 
                   "$\\theta$ (95\\% CI)", "$\\tau$", "$\\rho$", 
                   "$\\theta$ (95\\% CI)", "$\\tau$", "$\\rho$")
kbl(tab1, 
    format = "latex",
    longtable = F, 
    booktabs = T, 
    digits = 3,
    align = "r",
    linesep = c(rep("",9), "\\addlinespace"),
    escape = FALSE,
    position ="!htbp",
    caption = "Example 1: summary of the estimations of different sensitivity analysis methods",
    label = "tab1-1")%>% 
  add_header_above(c(" ","", 
                     "The proposed HN model based method" = 3, 
                     "The proposed CBN model based method" = 3
  ))%>%
  footnote(general = "$M$ indicates the number of potentially unpublished studies; 
           CI indicates the confidence interval.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")



tab2 = data.frame(
  p = sprintf("(%.2f, %.3f)", tab1_all$pnmin, tab1_all$pnmax),
  M1 = tab1_all$M.c1,
  CN1 = sprintf("%.3f (%.3f, %.3f)", lnOR_COPAS1999_1$mu, lnOR_COPAS1999_1$mu.lb, lnOR_COPAS1999_1$mu.ub),
  CN1.tau = sprintf("%.3f", lnOR_COPAS1999_1$tau),
  CN1.rho = sprintf("%.3f", lnOR_COPAS1999_1$rho),
  M2 = tab1_all$M.c2,
  CN2 = sprintf("%.3f (%.3f, %.3f)", lnOR_COPAS1999_2$mu, lnOR_COPAS1999_2$mu.lb, lnOR_COPAS1999_2$mu.ub),
  CN2.tau = sprintf("%.3f", lnOR_COPAS1999_2$tau),
  CN2.rho = sprintf("%.3f", lnOR_COPAS1999_2$rho)
)
colnames(tab2) = c("$(\\Pmin, \\Pmax)$", 
                   "$M_1$", 
                   "$\\theta$ (95\\% CI)", "$\\tau$", "$\\rho$", 
                   "$M_2$", 
                   "$\\theta$ (95\\% CI)", "$\\tau$", "$\\rho$")
kbl(tab2, 
    format = "latex",
    longtable = F, 
    booktabs = T, 
    digits = 3,
    align = "r",
    linesep = c(rep("",9), "\\addlinespace"),
    escape = FALSE,
    position ="!htbp",
    caption = "Example 1: summary of the estimations of different sensitivity analysis methods",
    label = "tab1-2")%>% 
  add_header_above(c(" ",
                     "The Copas-N method (only0)" = 4,
                     "The Copas-N method (all)" = 4))%>%
  footnote(general = "$M_1$ and $M_2$ indicate the number of potentially unpublished studies; 
  only 0 indicates continuity correction for only studies with 0 cells;
  all indicates continuity correction for all the studies; 
           CI indicates the confidence interval.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")



tab3 = data.frame(
  p = sprintf("(%.2f, %.3f)", tab1_all$pnmin, tab1_all$pnmax),
  M1 = tab1_all$M.c1,
  CN1 = sprintf("%.3f (%.3f, %.3f)", lnOR_COPAS2000_1$mu, lnOR_COPAS2000_1$mu.lb, lnOR_COPAS2000_1$mu.ub),
  CN1.tau = sprintf("%.3f", lnOR_COPAS2000_1$tau),
  CN1.rho = sprintf("%.3f", lnOR_COPAS2000_1$rho),
  M2 = tab1_all$M.c2,
  CN2 = sprintf("%.3f (%.3f, %.3f)", lnOR_COPAS2000_2$mu, lnOR_COPAS2000_2$mu.lb, lnOR_COPAS2000_2$mu.ub),
  CN2.tau = sprintf("%.3f", lnOR_COPAS2000_2$tau),
  CN2.rho = sprintf("%.3f", lnOR_COPAS2000_2$rho)
)
colnames(tab3) = c("$(\\Pmin, \\Pmax)$", 
                   "$M_1$", 
                   "$\\theta$ (95\\% CI)", "$\\tau$", "$\\rho$", 
                   "$M_2$", 
                   "$\\theta$ (95\\% CI)", "$\\tau$", "$\\rho$")
kbl(tab3, 
    format = "latex",
    longtable = F, 
    booktabs = T, 
    digits = 3,
    align = "r",
    linesep = c(rep("",9), "\\addlinespace"),
    escape = FALSE,
    position ="!htbp",
    caption = "Example 1: summary of the estimations of different sensitivity analysis methods",
    label = "tab1-3")%>% 
  add_header_above(c(" ",
                     "The Copas-Shi method (only0)" = 4,
                     "The Copas-Shi method (all)" = 4))%>%
  footnote(general = "$M_1$ and $M_2$ indicate the number of potentially unpublished studies; 
  only 0 indicates continuity correction for only studies with 0 cells;
  all indicates continuity correction for all the studies; 
           CI indicates the confidence interval.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")




tab4 = data.frame(
  p = sprintf("%.2f", res.t.all$p),
  M1 = res.t.all$M.t,
  CN1 = sprintf("%.3f (%.3f, %.3f)", res.t.all$hn.mean.1, res.t.all$hn.lower.1, res.t.all$hn.upper.1),
  CN1.tau = sprintf("%.3f", res.t.all$hn.mean.2),
  M2 = res.t.only0$M.t,
  CN2 = sprintf("%.3f (%.3f, %.3f)", res.t.only0$hn.mean.1, res.t.only0$hn.lower.1, res.t.only0$hn.upper.1),
  CN2.tau = sprintf("%.3f", res.t.only0$hn.mean.2)
)
colnames(tab4) = c("$p$", 
                   "$M_1$", 
                   "$\\theta$ (95\\% CI)", "$\\tau$",  
                   "$M_2$", 
                   "$\\theta$ (95\\% CI)", "$\\tau$")
kbl(tab4, 
    format = "latex",
    longtable = F, 
    booktabs = T, 
    digits = 3,
    align = "r",
    position ="!htbp",
    linesep = c(rep("",9), "\\addlinespace"),
    escape = FALSE,
    caption = "Example 1: summary of the estimations of different sensitivity analysis methods",
    label = "tab1-4")%>% 
  add_header_above(c(" ","",
                     "The t-statistic and HN model based method (only0)" = 2,
                     "The t-statistic and HN model based method (all)" = 2))%>%
  footnote(general = "$M$ indicates the number of potentially unpublished studies; 
  only 0 indicates continuity correction for only studies with 0 cells;
  all indicates continuity correction for all the studies; 
           CI indicates the confidence interval.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")
sink()