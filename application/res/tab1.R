
## TABLE
load("example-bias1.RData")

taba = data.frame(
  p = sprintf("(%.2f, %.3f)", tab1_all$pnmin, tab1_all$pnmax),
  Mp = tab1_all$M.p,
  HN = sprintf("%.3f (%.3f, %.3f)", tab1_all$HN.mu, tab1_all$HN.mu.lb, tab1_all$HN.mu.ub),
  HN.tau = sprintf("%.3f", lnOR_COPAS_HNGLMM[5,]),
  HN.rho = sprintf("%.3f", lnOR_COPAS_HNGLMM[7,]),
  BN = sprintf("%.3f (%.3f, %.3f)", tab1_all$BN.mu, tab1_all$BN.mu.lb, tab1_all$BN.mu.ub),
  BN.tau = sprintf("%.3f", lnOR_COPAS_BNGLMM[5,]),
  BN.rho = sprintf("%.3f", lnOR_COPAS_BNGLMM[7,]),
  CN = sprintf("%.3f (%.3f, %.3f)", tab1_all$CN.mu, tab1_all$CN.mu.lb, tab1_all$CN.mu.ub),
  CN.tau = sprintf("%.3f", lnOR_COPAS1999[5,]),
  CN.rho = sprintf("%.3f", lnOR_COPAS1999[7,]),
  Mc = tab1_all$M.c,
  CH = sprintf("%.3f (%.3f, %.3f)", tab1_all$CH.mu, tab1_all$CH.mu.lb, tab1_all$CH.mu.ub),
  CH.tau = sprintf("%.3f", lnOR_COPAS2000[5,]),
  CH.rho = sprintf("%.3f", lnOR_COPAS2000[7,])
)

colnames(taba) = c("$(\\Pmin, \\Pmax)$", "$M$", 
                   "$\\theta$ (95\\% CI)", "$\\tau$", "$\\rho$", 
                   "$\\theta$ (95\\% CI)", "$\\tau$", "$\\rho$", 
                   "$\\theta$ (95\\% CI)", "$\\tau$", "$\\rho$", 
                   "$M$", 
                   "$\\theta$ (95\\% CI)", "$\\tau$", "$\\rho$")
taba1 = taba[,1:8]
colnames(taba1) = c("$(\\Pmin, \\Pmax)$", "$M1$", 
                    "$\\theta$ (95\\% CI)", "$\\tau$ (SE)", "$\\rho$ (SE)", 
                    "$\\theta$ (95\\% CI)", "$\\tau$ (SE)", "$\\rho$ (SE)")
taba2 = taba[,c(1,2,9:15)]
colnames(taba2) = c("$(\\Pmin, \\Pmax)$", "$M1$", 
                    "$\\theta$ (95\\% CI)", "$\\tau$ (SE)", "$\\rho$ (SE)", 
                    "$M2$",
                    "$\\theta$ (95\\% CI)", "$\\tau$ (SE)", "$\\rho$ (SE)")

kbl(taba1, 
    format = "latex",
    longtable = F, 
    booktabs = T, 
    digits = 3,
    align = "r",
    linesep = c(rep("",9), "\\addlinespace"),
    escape = FALSE,
    caption = "Summary of the estimations of different sensitivity analysis methods",
    label = "tab1")%>% 
  add_header_above(c(" ","", 
                     "The proposed HN method" = 3, 
                     "The proposed BN method" = 3
  ))%>%
  footnote(general = "$M1$ indicates the number of potentially unpublished studies; 
           CI indicates the confidence interval; 
           SE indicates the standard error.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")

kbl(taba2, 
    format = "latex",
    longtable = F, 
    booktabs = T, 
    digits = 3,
    align = "r",
    linesep = c(rep("",9), "\\addlinespace"),
    escape = FALSE,
    caption = "Summary of the estimations of different sensitivity analysis methods",
    label = "tab2")%>% 
  add_header_above(c(" ","", 
                     "The Copas-N method" = 3,"",
                     "The Copas-Shi method" = 3))%>%
  footnote(general = "$M1$ and $M2$ indicates the number of potentially unpublished studies from the Copas-N and Copas-Shi methods; 
           CI indicates the confidence interval; 
           SE indicates the standard error.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")
