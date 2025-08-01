##
## PRINT THE RESULTS
##
rm(list=ls())
library(ggplot2)
library(gridExtra)
library(kableExtra)


## LOAD R DATA
load("example-bias1.RData")

ptheme= theme(panel.background = element_rect(fill = "white", colour = "grey50"),
             panel.grid.major = element_line(colour = "grey87"),
             axis.title = element_text(size = 14),
             plot.title = element_text(hjust = 0, size=14),
             legend.key = element_rect (fill = "white"),
             legend.position = c(0.5, 0.1),
             legend.title = element_blank(),
             legend.text = element_text(size = 14), 
             legend.background = element_rect(fill = "white", color = "black"))
py= scale_y_continuous(limits = c(-3,0), name = "lnOR", n.breaks = 10)
px1=scale_x_continuous(limits = c(0,40), name= "Missing studies",n.breaks = 5)
px2=scale_x_continuous(limits = c(0,40), name= "Missing studies",n.breaks = 5)
pguide=guide_legend(override.aes = list(lty = 1, size = 1))

## Pnmax = p
p1 = ggplot(tab1_all, aes(x = M.p)) +
  geom_ribbon(aes(ymin = HN.mu.lb, ymax = HN.mu.ub), alpha = 0.1, fill = "#377eb8", na.rm = TRUE) + 
  geom_line(aes(y = HN.mu.lb, colour="The proposed HN model based method"), lty=2, linewidth=1) +
  geom_line(aes(y = HN.mu.ub, colour="The proposed HN model based method"), lty=2, linewidth=1) +
  geom_point(aes(y = HN.mu, colour="The proposed HN model based method"), size=3) +
  geom_line(aes(y = HN.mu, colour="The proposed HN model based method"), lty=1, linewidth=1) +
  # px1 + 
  py + ptheme+
  labs(title = "(A)")+  
  scale_colour_manual(breaks = "The proposed HN model based method", values = "#377eb8", guide = pguide)


p2 = ggplot(tab1_all, aes(x = M.p)) +
  geom_ribbon(aes(ymin = BN.mu.lb, ymax = BN.mu.ub), alpha = 0.1, fill = "#e41a1c", na.rm = TRUE) + 
  geom_line(aes(y = BN.mu.lb, colour="The proposed BN model based method"), lty=2, size=1) +
  geom_line(aes(y = BN.mu.ub, colour="The proposed BN model based method"), lty=2, size=1) +
  geom_point(aes(y = BN.mu, colour="The proposed BN model based method"), size=3) +
  geom_line(aes(y = BN.mu, colour="The proposed BN model based method"), lty=1, size=1) +
  px1 + 
  py + ptheme+
  labs(title = "(B)")+  
  scale_colour_manual(breaks = "The proposed BN model based method", values = "#e41a1c", guide = pguide)

p3 = ggplot(tab1_all, aes(x = M.p)) +
  geom_ribbon(aes(ymin = CN1.mu.lb, ymax = CN1.mu.ub), alpha = 0.1, fill = "#4daf4a", na.rm = TRUE) + 
  geom_line(aes(y = CN1.mu.lb, colour="The Copas-N method (only0)"),lty=2, size=1) +
  geom_line(aes(y = CN1.mu.ub, colour="The Copas-N method (only0)"),lty=2, size=1) +
  geom_point(aes(y = CN1.mu, colour="The Copas-N method (only0)"), size=3) +
  geom_line(aes(y = CN1.mu, colour="The Copas-N method (only0)"),lty=1, size=1) +
  px1 + py + ptheme +
  labs(title = "(C)")+  
  scale_colour_manual(breaks = c("The Copas-N method (only0)"), 
                      values = c("#4daf4a"), guide = pguide)

p4 = ggplot(tab1_all, aes(x = M.c1)) +
  geom_ribbon(aes(ymin = CH1.mu.lb, ymax = CH1.mu.ub), alpha = 0.1, fill = "#984ea3", na.rm = TRUE) + 
  geom_line(aes(y = CH1.mu.lb, colour="The Copas-Shi method (only0)"),lty=2, size=1) +
  geom_line(aes(y = CH1.mu.ub, colour="The Copas-Shi method (only0)"),lty=2, size=1) +
  geom_point(aes(y = CH1.mu, colour="The Copas-Shi method (only0)"), size=3) +
  geom_line(aes(y = CH1.mu, colour="The Copas-Shi method (only0)"),lty=1, size=1) +
  geom_ribbon(aes(ymin = CH2.mu.lb, ymax = CH2.mu.ub), alpha = 0.1, fill = "#beaed4", na.rm = TRUE) + 
  geom_line(aes(y = CH2.mu.lb, colour="The Copas-Shi method (all)"),lty=2, size=1) +
  geom_line(aes(y = CH2.mu.ub, colour="The Copas-Shi method (all)"),lty=2, size=1) +
  geom_point(aes(y = CH2.mu, colour="The Copas-Shi method (all)"), size=3) +
  geom_line(aes(y = CH2.mu, colour="The Copas-Shi method (all)"),lty=1, size=1) +
  px2 + py + ptheme +
  labs(title = "(D)")+  
  scale_colour_manual(breaks = c("The Copas-Shi method (only0)","The Copas-Shi method (all)"), 
                      values = c("#984ea3","#beaed4"), guide = pguide)



p=grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
ggsave(filename = "plot.eps", plot = p, device = cairo_ps, width = 12, height = 12) 


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
