## TABLE
load("example-bias2.RData")

taba = data.frame(
  p = sprintf("(%.2f, %.3f)", tab1_all$pnmin, tab1_all$pnmax),
  Mp = tab1_all$M.p,
  BN = sprintf("%.3f (%.3f, %.3f)", tab1_all$BN.mu, tab1_all$BN.mu.lb, tab1_all$BN.mu.ub),
  BN.tau = sprintf("%.3f", lgtP_COPAS_BNGLMM[5,]),
  BN.rho = sprintf("%.3f", lgtP_COPAS_BNGLMM[7,]),
  CN = sprintf("%.3f (%.3f, %.3f)", tab1_all$CN.mu, tab1_all$CN.mu.lb, tab1_all$CN.mu.ub),
  CN.tau = sprintf("%.3f", lgtP_COPAS1999[5,]),
  CN.rho = sprintf("%.3f", lgtP_COPAS1999[7,]),
  Mc = tab1_all$M.c,
  CH = sprintf("%.3f (%.3f, %.3f)", tab1_all$CH.mu, tab1_all$CH.mu.lb, tab1_all$CH.mu.ub),
  CH.tau = sprintf("%.3f", lgtP_COPAS2000[5,]),
  CH.rho = sprintf("%.3f", lgtP_COPAS2000[7,])
)
colnames(taba) = c("$(\\Pmin, \\Pmax)$", "$M$", 
                   "$\\theta$ (95\\% CI)", "$\\tau$", "$\\rho$", 
                   "$M$", "$\\theta$ (95\\% CI)", "$\\tau$", "$\\rho$")
taba1 = taba[,1:5]
colnames(taba1) = c("$(\\Pmin, \\Pmax)$", "$M1$", 
                    "$\\theta$ (95\\% CI)", "$\\tau$", "$\\rho$")
taba2 = taba[,c(1:2,6:12)]
colnames(taba2) = c("$(\\Pmin, \\Pmax)$", "$M1$", 
                    "$\\theta$ (95\\% CI)", "$\\tau$", "$\\rho$",
                    "$M2$", 
                    "$\\theta$ (95\\% CI)", "$\\tau$", "$\\rho$")
kbl(taba1, 
    format = "latex",
    longtable = F, 
    booktabs = T, 
    digits = 3,
    align = "r",
    linesep = c(rep("",9), "\\addlinespace"),
    escape = FALSE,
    caption = "Summary of the estimations of different sensitivity analysis methods",
    label = "tab3")%>% 
  add_header_above(c("", "","The proposed BN method" = 3))%>%
  footnote(general = "$M1$ indicates the number of potentially unpublished studies; 
           CI indicates the confidence interval; 
           SE indicates the standard error;
           NaN indicates an undefined value.", 
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
    label = "tab4")%>% 
  add_header_above(c("", "",
                     "The Copas-N method" = 3,"",
                     "Copas-Shi method" = 3))%>%
  footnote(general = "$M1$ and $M2$ indicates the number of potentially unpublished studies from the Copas-N and Copas-Shi methods; 
           CI indicates the confidence interval; 
           SE indicates the standard error;
           NaN indicates an undefined value.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")