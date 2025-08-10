rm(list = ls())
sink("output.txt") 
library(magrittr)

## TABLE
load("app4-nopb.RData")
load("app4-prop.RData")
load("app4-t-prop.RData")


tab1 = data.frame(
  p = sprintf("(%.2f, %.3f)", tab1_all$pnmin, tab1_all$pnmax),
  Mp = tab1_all$M.p,
  HN = sprintf("%.3f (%.3f, %.3f)", tab1_all$BN.mu, tab1_all$BN.mu.lb, tab1_all$BN.mu.ub),
  HN.tau = sprintf("%.3f", lgtP_COPAS_BNGLMM$tau),
  HN.rho = sprintf("%.3f", lgtP_COPAS_BNGLMM$rho)
)
colnames(tab1) = c("$(\\Pmin, \\Pmax)$", "$M$", 
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
                     "The proposed 1SBN model based method" = 3
  ))%>%
  footnote(general = "$M$ indicates the number of potentially unpublished studies; 
           CI indicates the confidence interval.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")



tab2 = data.frame(
  p = sprintf("(%.2f, %.3f)", tab1_all$pnmin, tab1_all$pnmax),
  M1 = tab1_all$M.c1,
  CN1 = sprintf("%.3f (%.3f, %.3f)", lgtP_COPAS1999_1$mu, lgtP_COPAS1999_1$mu.lb, lgtP_COPAS1999_1$mu.ub),
  CN1.tau = sprintf("%.3f", lgtP_COPAS1999_1$tau),
  CN1.rho = sprintf("%.3f", lgtP_COPAS1999_1$rho),
  M2 = tab1_all$M.c2,
  CN2 = sprintf("%.3f (%.3f, %.3f)", lgtP_COPAS1999_2$mu, lgtP_COPAS1999_2$mu.lb, lgtP_COPAS1999_2$mu.ub),
  CN2.tau = sprintf("%.3f", lgtP_COPAS1999_2$tau),
  CN2.rho = sprintf("%.3f", lgtP_COPAS1999_2$rho)
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
  CN1 = sprintf("%.3f (%.3f, %.3f)", lgtP_COPAS2000_1$mu, lgtP_COPAS2000_1$mu.lb, lgtP_COPAS2000_1$mu.ub),
  CN1.tau = sprintf("%.3f", lgtP_COPAS2000_1$tau),
  CN1.rho = sprintf("%.3f", lgtP_COPAS2000_1$rho),
  M2 = tab1_all$M.c2,
  CN2 = sprintf("%.3f (%.3f, %.3f)", lgtP_COPAS2000_2$mu, lgtP_COPAS2000_2$mu.lb, lgtP_COPAS2000_2$mu.ub),
  CN2.tau = sprintf("%.3f", lgtP_COPAS2000_2$tau),
  CN2.rho = sprintf("%.3f", lgtP_COPAS2000_2$rho)
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
  CN1 = sprintf("%.3f (%.3f, %.3f)", res.t.all$bn.mean.1, res.t.all$bn.lower.1, res.t.all$bn.upper.1),
  CN1.tau = sprintf("%.3f", res.t.all$bn.mean.2),
  CN2 = sprintf("%.3f (%.3f, %.3f)", res.t.only0$bn.mean.1, res.t.only0$bn.lower.1, res.t.only0$bn.upper.1),
  CN2.tau = sprintf("%.3f", res.t.only0$bn.mean.2)
)
colnames(tab4) = c("$p$", 
                   "$M$", 
                   "$\\theta$ (95\\% CI)", "$\\tau$",  
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