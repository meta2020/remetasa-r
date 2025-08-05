##
## PRINT THE RESULTS
##
rm(list=ls())
library(ggplot2)
library(gridExtra)
library(kableExtra)

## LOAD R DATA
load("app4-prop.RData")
load("app3-t-prop.RData")
res.t=cbind.data.frame(only0=res.t.only0, all=res.t.all)

ptheme= theme(panel.background = element_rect(fill = "white", colour = "grey50"),
              panel.grid.major = element_line(colour = "grey87"),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0, size=14),
              legend.key = element_rect (fill = "white"),
              legend.position = c(0.5, 0.1),
              legend.title = element_blank(),
              legend.text = element_text(size = 14), 
              legend.background = element_rect(fill = "white", color = "black"))
py= scale_y_continuous(limits = c(-8,0), name = "lnOR", n.breaks = 10)
py2= scale_y_continuous(limits = c(-3.5,0.5), name = "lnOR", n.breaks = 10)
px1=scale_x_reverse(n.breaks = 10, name="P(publishing studies with smallest sample size)")
px2=scale_x_reverse(n.breaks = 10, name="P(publishing studies with largest SE)")
px3=scale_x_reverse(n.breaks = 10, name="P(publishing studies from population)")
pguide=guide_legend(override.aes = list(lty = 1, size = 1))

## Pnmax = p
p1 = ggplot(tab1_all, aes(x = pnmin)) +
  geom_ribbon(aes(ymin = BN.mu.lb, ymax = BN.mu.ub), alpha = 0.1, fill = "#d7191c", na.rm = TRUE) + 
  geom_line(aes(y = BN.mu.lb, colour="The proposed BN model based method"), lty=2, size=1) +
  geom_line(aes(y = BN.mu.ub, colour="The proposed BN model based method"), lty=2, size=1) +
  geom_point(aes(y = BN.mu, colour="The proposed BN model based method"), size=3) +
  geom_line(aes(y = BN.mu, colour="The proposed BN model based method"), lty=1, size=1) +
  px1 + py + ptheme+
  geom_text(aes(y = BN.mu, label = as.character(M.p)), vjust = -2)+
  labs(title = "(A)")+  
  scale_colour_manual(breaks = "The proposed BN model based method", values = "#e41a1c", guide = pguide)

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
  geom_text(aes(y = CN2.mu, label = as.character(M.p)), vjust = -2)+
  labs(title = "(B)")+  
  scale_colour_manual(breaks = c("The Copas-N method (only0)","The Copas-N method (all)"), 
                      values = c("#4daf4a","#b2df8a"), guide = pguide)

p3 = ggplot(tab1_all, aes(x = pnmin)) +
  geom_ribbon(aes(ymin = CH1.mu.lb, ymax = CH1.mu.ub), alpha = 0.1, fill = "#984ea3", na.rm = TRUE) + 
  geom_line(aes(y = CH1.mu.lb, colour="The Copas-Shi method (only0)"), lty=2, size=1) +
  geom_line(aes(y = CH1.mu.ub, colour="The Copas-Shi method (only0)"), lty=2, size=1) +
  geom_point(aes(y = CH1.mu, colour="The Copas-Shi method (only0)"), size=3) +
  geom_line(aes(y = CH1.mu, colour="The Copas-Shi method (only0)"), lty=1, size=1) +
  geom_text(aes(y = CH1.mu, label = as.character(M.c1)), vjust = 2)+
  geom_ribbon(aes(ymin = CH2.mu.lb, ymax = CH2.mu.ub), alpha = 0.1, fill = "#beaed4", na.rm = TRUE) + 
  geom_line(aes(y = CH2.mu.lb, colour="The Copas-Shi method (all)"), lty=2, size=1) +
  geom_line(aes(y = CH2.mu.ub, colour="The Copas-Shi method (all)"), lty=2, size=1) +
  geom_point(aes(y = CH2.mu, colour="The Copas-Shi method (all)"), size=3) +
  geom_line(aes(y = CH2.mu, colour="The Copas-Shi method (all)"), lty=1, size=1) +
  geom_text(aes(y = CH2.mu, label = as.character(M.c2)), vjust = -1.5)+
  px2 + py + ptheme+
  labs(title = "(C)")+  
  scale_colour_manual(breaks = c("The Copas-Shi method (only0)","The Copas-Shi method (all)"), 
                      values = c("#984ea3","#beaed4"), guide = pguide)

p4 = ggplot(res.t, aes(x = only0.p)) +
  geom_ribbon(aes(ymin = only0.bn.lower.1, ymax = only0.bn.upper.1), alpha = 0.1, fill = "#ff7f00", na.rm = TRUE) + 
  geom_line(aes(y = only0.bn.lower.1, colour="The t-statistic and BN model based method (only0)"), lty=2, size=1) +
  geom_line(aes(y = only0.bn.upper.1, colour="The t-statistic and BN model based method (only0)"), lty=2, size=1) +
  geom_point(aes(y = only0.bn.mean.1, colour="The t-statistic and BN model based method (only0)"), size=3) +
  geom_line(aes(y = only0.bn.mean.1, colour="The t-statistic and BN model based method (only0)"), lty=1, size=1) +
  geom_ribbon(aes(ymin = all.bn.lower.1, ymax = all.bn.upper.1), alpha = 0.1, fill = "#fdc086", na.rm = TRUE) + 
  geom_line(aes(y = all.bn.lower.1, colour="The t-statistic and BN model based method (all)"), lty=2, size=1) +
  geom_line(aes(y = all.bn.upper.1, colour="The t-statistic and BN model based method (all)"), lty=2, size=1) +
  geom_point(aes(y = all.bn.mean.1, colour="The t-statistic and BN model based method (all)"), size=3) +
  geom_line(aes(y = all.bn.mean.1, colour="The t-statistic and BN model based method (all)"), lty=1, size=1) +
  px3 + py2 + ptheme+
  geom_text(aes(y = all.bn.mean.1, label = as.character(only0.M.t)), vjust = -2)+
  labs(title = "(D)")+  
  scale_colour_manual(breaks = c("The t-statistic and BN model based method (only0)","The t-statistic and BN model based method (all)"), 
                      values = c("#ff7f00","#fdc086"), guide = pguide)

p=grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
ggsave(filename = "plot3.eps", plot = p, device = cairo_ps, width = 12, height = 12) 



