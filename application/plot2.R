##
## PRINT THE RESULTS
##

library(ggplot2)
library(gridExtra)
library(kableExtra)

## LOAD R DATA
load("example-bias2.RData")
load("example2-t-all.RData")
load("example2-t-only0.RData")
res.t=cbind.data.frame(only0=res.all.t0, all=res.t.all)

ptheme= theme(panel.background = element_rect(fill = "white", colour = "grey50"),
              panel.grid.major = element_line(colour = "grey87"),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0, size=14),
              legend.key = element_rect (fill = "white"),
              legend.position = c(0.5, 0.1),
              legend.title = element_blank(),
              legend.text = element_text(size = 14), 
              legend.background = element_rect(fill = "white", color = "black"))
py= scale_y_continuous(limits = c(-7,-2), name = "lnOR", n.breaks = 10)
py2= scale_y_continuous(limits = c(-5,5), name = "lnOR", n.breaks = 10)
px1=scale_x_reverse(n.breaks = 10, name="P(publishing studies with smallest sample size)")
px2=scale_x_reverse(n.breaks = 10, name="P(publishing studies with largest SE)")
px3=scale_x_reverse(n.breaks = 10, name="P(publishing studies from population)")
pguide=guide_legend(override.aes = list(lty = 1, size = 1))

## Pnmax = p
p1 = ggplot(tab1_all, aes(x = pnmin)) +
  geom_ribbon(aes(ymin = BN.mu.lb, ymax = BN.mu.ub), alpha = 0.1, fill = "#d7191c", na.rm = TRUE) + 
  geom_line(aes(y = BN.mu.lb, colour="The proposed BN method"), lty=2, size=1) +
  geom_line(aes(y = BN.mu.ub, colour="The proposed BN method"), lty=2, size=1) +
  geom_point(aes(y = BN.mu, colour="The proposed BN method"), size=3) +
  geom_line(aes(y = BN.mu, colour="The proposed BN method"), lty=1, size=1) +
  px1 + py + ptheme+
  labs(title = "(A)")+  
  scale_colour_manual(breaks = "The proposed BN method", values = "#e41a1c", guide = pguide)

p2 = ggplot(tab1_all, aes(x = pnmin)) +
  geom_ribbon(aes(ymin = CN1.mu.lb, ymax = CN1.mu.ub), alpha = 0.1, fill = "#4daf4a", na.rm = TRUE) + 
  geom_line(aes(y = CN1.mu.lb, colour="The Copas-N method (only0)"), lty=2, size=1) +
  geom_line(aes(y = CN1.mu.ub, colour="The Copas-N method (only0)"), lty=2, size=1) +
  geom_point(aes(y = CN1.mu, colour="The Copas-N method (only0)"), size=3) +
  geom_line(aes(y = CN1.mu, colour="The Copas-N method (only0)"), lty=1, size=1) +
  geom_ribbon(aes(ymin = CN2.mu.lb, ymax = CN2.mu.ub), alpha = 0.1, fill = "#b2df8a", na.rm = TRUE) + 
  geom_line(aes(y = CN2.mu.lb, colour="The Copas-N method (all)"), lty=2, size=1) +
  geom_line(aes(y = CN2.mu.ub, colour="The Copas-N method (all)"), lty=2, size=1) +
  geom_point(aes(y = CN2.mu, colour="The Copas-N method (all)"), size=3) +
  geom_line(aes(y = CN2.mu, colour="The Copas-N method (all)"), lty=1, size=1) +
  px1 + py + ptheme+
  labs(title = "(B)")+  
  scale_colour_manual(breaks = c("The Copas-N method (only0)","The Copas-N method (all)"), 
                      values = c("#4daf4a","#b2df8a"), guide = pguide)

p3 = ggplot(tab1_all, aes(x = pnmin)) +
  geom_ribbon(aes(ymin = CH1.mu.lb, ymax = CH1.mu.ub), alpha = 0.1, fill = "#984ea3", na.rm = TRUE) + 
  geom_line(aes(y = CH1.mu.lb, colour="The Copas-Shi method (only0)"), lty=2, size=1) +
  geom_line(aes(y = CH1.mu.ub, colour="The Copas-Shi method (only0)"), lty=2, size=1) +
  geom_point(aes(y = CH1.mu, colour="The Copas-Shi method (only0)"), size=3) +
  geom_line(aes(y = CH1.mu, colour="The Copas-Shi method (only0)"), lty=1, size=1) +
  geom_ribbon(aes(ymin = CH2.mu.lb, ymax = CH2.mu.ub), alpha = 0.1, fill = "#beaed4", na.rm = TRUE) + 
  geom_line(aes(y = CH2.mu.lb, colour="The Copas-Shi method (all)"), lty=2, size=1) +
  geom_line(aes(y = CH2.mu.ub, colour="The Copas-Shi method (all)"), lty=2, size=1) +
  geom_point(aes(y = CH2.mu, colour="The Copas-Shi method (all)"), size=3) +
  geom_line(aes(y = CH2.mu, colour="The Copas-Shi method (all)"), lty=1, size=1) +
  px2 + py + ptheme+
  labs(title = "(C)")+  
  scale_colour_manual(breaks = c("The Copas-Shi method (only0)","The Copas-Shi method (all)"), 
                      values = c("#984ea3","#beaed4"), guide = pguide)

p4 = ggplot(res.t, aes(x = only0.p)) +
  geom_ribbon(aes(ymin = only0.bn.lower.1, ymax = only0.bn.upper.1), alpha = 0.1, fill = "#ff7f00", na.rm = TRUE) + 
  geom_line(aes(y = only0.bn.lower.1, colour="The t-statistic based method (only0)"), lty=2, size=1) +
  geom_line(aes(y = only0.bn.upper.1, colour="The t-statistic based method (only0)"), lty=2, size=1) +
  geom_point(aes(y = only0.bn.mean.1, colour="The t-statistic based method (only0)"), size=3) +
  geom_line(aes(y = only0.bn.mean.1, colour="The t-statistic based method (only0)"), lty=1, size=1) +
  geom_ribbon(aes(ymin = all.bn.lower.1, ymax = all.bn.upper.1), alpha = 0.1, fill = "#fdc086", na.rm = TRUE) + 
  geom_line(aes(y = all.bn.lower.1, colour="The t-statistic based method (all)"), lty=2, size=1) +
  geom_line(aes(y = all.bn.upper.1, colour="The t-statistic based method (all)"), lty=2, size=1) +
  geom_point(aes(y = all.bn.mean.1, colour="The t-statistic based method (all)"), size=3) +
  geom_line(aes(y = all.bn.mean.1, colour="The t-statistic based method (all)"), lty=1, size=1) +
  px3 + py2 + ptheme+
  labs(title = "(C)")+  
  scale_colour_manual(breaks = c("The t-statistic based method (only0)","The t-statistic based method (all)"), 
                      values = c("#ff7f00","#fdc086"), guide = pguide)

p=grid.arrange(p1, p2, p3, ncol = 2, nrow = 2)
ggsave(filename = "plot2.eps", plot = p, device = cairo_ps, width = 12, height = 12) 


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
