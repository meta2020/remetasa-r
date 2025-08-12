# Load required libraries
library(ggplot2)
library(grid)
library(gridBase)
library(kableExtra)
library(metafor)

# Create 4 ggplot2 plots
load("app4-prop.RData")
load("app4-t-prop.RData")

data = read.csv("../pritz1997.csv")
yvi1 = escalc(measure="PLO", xi=y1, ni=n1, data=data, to="only0")
yi1 = yvi1[,1]
vi1 = yvi1[,2]

yvi2 = escalc(measure="PLO", xi=y1, ni=n1, data=data, to="all")
yi2 = yvi2[,1]
vi2 = yvi2[,2]


## ggplot----

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
py= scale_y_continuous(limits = c(-3.5,0.5), name = "Log odds", n.breaks = 10)
py2= scale_y_continuous(limits = c(-3.5,0.5), name = "Log odds", n.breaks = 10)
px1=scale_x_reverse(n.breaks = 10, name="P(publishing studies with smallest sample size)")
px2=scale_x_reverse(n.breaks = 10, name="P(publishing studies with largest SE)")
px3=scale_x_reverse(n.breaks = 10, name="P(publishing studies from population)")
pguide=guide_legend(override.aes = list(lty = 1, size = 1))

p1 = ggplot(tab1_all, aes(x = pnmin)) +
  geom_ribbon(aes(ymin = BN.mu.lb, ymax = BN.mu.ub), alpha = 0, fill = "#d7191c", na.rm = TRUE) + 
  geom_line(aes(y = BN.mu.lb, colour="The proposed 1SBN model based method"), lty=2, size=1) +
  geom_line(aes(y = BN.mu.ub, colour="The proposed 1SBN model based method"), lty=2, size=1) +
  geom_point(aes(y = BN.mu, colour="The proposed 1SBN model based method"), size=3) +
  geom_line(aes(y = BN.mu, colour="The proposed 1SBN model based method"), lty=1, size=1) +
  px1 + py + ptheme+
  geom_text(aes(y = BN.mu, label = as.character(M.p)), vjust = -2)+
  labs(title = "C")+  
  scale_colour_manual(breaks = "The proposed 1SBN model based method", values = "#e41a1c", guide = pguide)

p2 = ggplot(tab1_all, aes(x = pnmin)) +
  geom_ribbon(aes(ymin = CN1.mu.lb, ymax = CN1.mu.ub), alpha = 0, fill = "#4daf4a", na.rm = TRUE) + 
  geom_line(aes(y = CN1.mu.lb, colour="The Copas-N method (only0)"), lty=2, size=1) +
  geom_line(aes(y = CN1.mu.ub, colour="The Copas-N method (only0)"), lty=2, size=1) +
  geom_point(aes(y = CN1.mu, colour="The Copas-N method (only0)"), size=3) +
  geom_line(aes(y = CN1.mu, colour="The Copas-N method (only0)"), lty=1, size=1) +
  geom_ribbon(aes(ymin = CN2.mu.lb, ymax = CN2.mu.ub), alpha = 0, fill = "#b2df8a", na.rm = TRUE) + 
  geom_line(aes(y = CN2.mu.lb, colour="The Copas-N method (all)"), lty=2, size=1) +
  geom_line(aes(y = CN2.mu.ub, colour="The Copas-N method (all)"), lty=2, size=1) +
  geom_point(aes(y = CN2.mu, colour="The Copas-N method (all)"), size=3) +
  geom_line(aes(y = CN2.mu, colour="The Copas-N method (all)"), lty=1, size=1) +
  px1 + py + ptheme+
  geom_text(aes(y = CN2.mu, label = as.character(M.p)), vjust = -2)+
  labs(title = "D")+  
  scale_colour_manual(breaks = c("The Copas-N method (only0)","The Copas-N method (all)"), 
                      values = c("#4daf4a","#b2df8a"), guide = pguide)

p3 = ggplot(tab1_all, aes(x = pnmin)) +
  geom_ribbon(aes(ymin = CH1.mu.lb, ymax = CH1.mu.ub), alpha = 0, fill = "#984ea3", na.rm = TRUE) + 
  geom_line(aes(y = CH1.mu.lb, colour="The Copas-Shi method (only0)"), lty=2, size=1) +
  geom_line(aes(y = CH1.mu.ub, colour="The Copas-Shi method (only0)"), lty=2, size=1) +
  geom_point(aes(y = CH1.mu, colour="The Copas-Shi method (only0)"), size=3) +
  geom_line(aes(y = CH1.mu, colour="The Copas-Shi method (only0)"), lty=1, size=1) +
  geom_ribbon(aes(ymin = CH2.mu.lb, ymax = CH2.mu.ub), alpha = 0, fill = "#beaed4", na.rm = TRUE) + 
  geom_line(aes(y = CH2.mu.lb, colour="The Copas-Shi method (all)"), lty=2, size=1) +
  geom_line(aes(y = CH2.mu.ub, colour="The Copas-Shi method (all)"), lty=2, size=1) +
  geom_point(aes(y = CH2.mu, colour="The Copas-Shi method (all)"), size=3) +
  geom_line(aes(y = CH2.mu, colour="The Copas-Shi method (all)"), lty=1, size=1) +
  px2 + py + ptheme+
  geom_text(aes(y = CH1.mu, label = as.character(M.c1)), vjust = 2)+
  geom_text(aes(y = CH2.mu, label = as.character(M.c2)), vjust = -1.5)+
  labs(title = "E")+  
  scale_colour_manual(breaks = c("The Copas-Shi method (only0)","The Copas-Shi method (all)"), 
                      values = c("#984ea3","#beaed4"), guide = pguide)

p4 = ggplot(res.t, aes(x = only0.p)) +
  geom_ribbon(aes(ymin = only0.bn.lower.1, ymax = only0.bn.upper.1), alpha = 0, fill = "#ff7f00", na.rm = TRUE) + 
  geom_line(aes(y = only0.bn.lower.1, colour="The t-statistic method (only0)"), lty=2, size=1) +
  geom_line(aes(y = only0.bn.upper.1, colour="The t-statistic method (only0)"), lty=2, size=1) +
  geom_point(aes(y = only0.bn.mean.1, colour="The t-statistic method (only0)"), size=3) +
  geom_line(aes(y = only0.bn.mean.1, colour="The t-statistic method (only0)"), lty=1, size=1) +
  geom_ribbon(aes(ymin = all.bn.lower.1, ymax = all.bn.upper.1), alpha = 0, fill = "#fdc086", na.rm = TRUE) + 
  geom_line(aes(y = all.bn.lower.1, colour="The t-statistic method (all)"), lty=2, size=1) +
  geom_line(aes(y = all.bn.upper.1, colour="The t-statistic method (all)"), lty=2, size=1) +
  geom_point(aes(y = all.bn.mean.1, colour="The t-statistic method (all)"), size=3) +
  geom_line(aes(y = all.bn.mean.1, colour="The t-statistic method (all)"), lty=1, size=1) +
  px3 + py2 + ptheme+
  geom_text(aes(y = all.bn.mean.1, label = as.character(only0.M.t)), vjust = -2)+
  labs(title = "F")+  
  scale_colour_manual(breaks = c("The t-statistic method (only0)","The t-statistic method (all)"), 
                      values = c("#ff7f00","#fdc086"), guide = pguide)




# Open EPS device----
setEPS()

postscript("plot4.eps", width = 12, height = 14)


# Start new page with 3x2 layout
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 2)))  # 3 rows, 2 columns

# Function to draw a base plot in a specified cell
draw_base_plot <- function(row, col, plot_fun) {
  pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
  par(fig=gridFIG())
  par(new=TRUE)
  plot_fun()
  popViewport()
}

# Function to draw ggplot2 plot in a specified cell
draw_ggplot <- function(row, col, ggplot_obj) {
  print(ggplot_obj, vp = viewport(layout.pos.row = row, layout.pos.col = col))
}

# --- First row: base R plots ---
draw_base_plot(1, 1, function() {
  res1 = rma(yi, vi, data=yvi1)
  funnel(trimfill(res1,estimator="L0"))
  abline(v=res1$beta)
  reg1 = regtest(res1, model="lm")
  rnk1 = ranktest(res1)
  mtext("A Continuity correction for studies with 0 cells", side = 3, line = 2, adj = 0)
  mtext(sprintf("Regression test: t = %.3f, p = %.3f", reg1$zval, reg1$pval),side = 3, line = 1, adj = 0)
  mtext(sprintf("Rank test: t = %.3f, p = %.3f", rnk1$tau, rnk1$pval),side = 3, adj = 0)
})
draw_base_plot(1, 2, function() {
  res2 = rma(yi, vi, data=yvi2)
  funnel(trimfill(res2,estimator="L0"))
  abline(v=res2$beta)
  reg2 = regtest(res2, model="lm")
  rnk2 = ranktest(res2)
  mtext("B Continuity correction for all studies", side = 3, line = 2, adj = 0)
  mtext(sprintf("Regression test: t = %.3f, p = %.3f", reg2$zval, reg2$pval),side = 3, line = 1, adj = 0)
  mtext(sprintf("Rank test: t = %.3f, p = %.3f", rnk2$tau, rnk2$pval),side = 3, adj = 0)
})

# --- Second and third rows: ggplot2 plots ---
draw_ggplot(2, 1, p1)
draw_ggplot(2, 2, p2)
draw_ggplot(3, 1, p3)
draw_ggplot(3, 2, p4)

# grid.newpage()
# Close EPS device
dev.off()

