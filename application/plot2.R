##
## PRINT THE RESULTS
##

library(ggplot2)
library(gridExtra)
library(kableExtra)

## LOAD R DATA
load("example-bias2.RData")

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
px1=scale_x_reverse(n.breaks = 10, name="P(publishing studies with smallest sample size)")
px2=scale_x_reverse(n.breaks = 10, name="P(publishing studies with largest SE)")
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
  geom_ribbon(aes(ymin = CN.mu.lb, ymax = CN.mu.ub), alpha = 0.1, fill = "#4daf4a", na.rm = TRUE) + 
  geom_line(aes(y = CN.mu.lb, colour="The Copas-N method"), lty=2, size=1) +
  geom_line(aes(y = CN.mu.ub, colour="The Copas-N method"), lty=2, size=1) +
  geom_point(aes(y = CN.mu, colour="The Copas-N method"), size=3) +
  geom_line(aes(y = CN.mu, colour="The Copas-N method"), lty=1, size=1) +
  px1 + py + ptheme+
  labs(title = "(B)")+  
  scale_colour_manual(breaks = "The Copas-N method", values = "#4daf4a", guide = pguide)

p3 = ggplot(tab1_all, aes(x = pnmin)) +
  geom_ribbon(aes(ymin = CH.mu.lb, ymax = CH.mu.ub), alpha = 0.1, fill = "#984ea3", na.rm = TRUE) + 
  geom_line(aes(y = CH.mu.lb, colour="The Copas-Shi method"), lty=2, size=1) +
  geom_line(aes(y = CH.mu.ub, colour="The Copas-Shi method"), lty=2, size=1) +
  geom_point(aes(y = CH.mu, colour="The Copas-Shi method"), size=3) +
  geom_line(aes(y = CH.mu, colour="The Copas-Shi method"), lty=1, size=1) +
  px2 + py + ptheme+
  labs(title = "(C)")+  
  scale_colour_manual(breaks = "The Copas-Shi method", values = "#984ea3", guide = pguide)

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
