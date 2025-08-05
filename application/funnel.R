#'
#' Create funnel plot
#'
#' Load R functions
rm(list=ls())


library(metafor)

# data = read.csv("niel-weise21.csv")
data = read.csv("thomas.csv")
# 
data=dat.pritz1997[,-2]
colnames(data)=c("study","y0","n1")
data$y1=data$n1-data$y0

#' Meta-analysis without PB ----------
#' Data

## meta-analysis of lnORs
#' Derive continuous outcomes (lnOR and se)
yvi1 = escalc(measure="OR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0, data=data, to="only0")
yi1 = yvi1[,1]
vi1 = yvi1[,2]

yvi2 = escalc(measure="OR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0, data=data, to="all")
yi2 = yvi2[,1]
vi2 = yvi2[,2]

# Funnel plot
postscript("funnel1.eps", width = 12, height = 12)
par(mfcol=c(2,2))

res1 = rma(yi, vi, data=yvi1, method="ML")
funnel(trimfill(res1,estimator="L0"), xlim = c(-6,4))
abline(v=res1$beta)
reg1 = regtest(res1, model="lm")
rnk1 = ranktest(res1)
mtext("A. L0 and continuity correction for only studies with 0 ", side = 3, adj = 0)
text(-1,0.05, sprintf("Regression test for funnel plot asymmetry: t = %.3f, p = %.3f", reg1$zval, reg1$pval))
text(-1,0.2, sprintf("Rank test for funnel plot asymmetry: t = %.3f, p = %.3f", rnk1$tau, rnk1$pval))

res2 = rma(yi, vi, data=yvi2, method="ML")
funnel(trimfill(res2,estimator="L0"), xlim = c(-6,4))
abline(v=res2$beta)
reg2 = regtest(res2, model="lm")
rnk2 = ranktest(res2)
mtext("B. L0 and continuity correction for all studies", side = 3, adj = 0)
text(-1,0.05, sprintf("Regression test for funnel plot asymmetry: t = %.3f, p = %.3f", reg2$zval, reg2$pval))
text(-1,0.2, sprintf("Rank test for funnel plot asymmetry: t = %.3f, p = %.3f", rnk2$tau, rnk2$pval))


res3 = rma(yi, vi, data=yvi1, method="ML")
funnel(trimfill(res3,estimator="R0"), xlim = c(-6,4))
abline(v=res3$beta)
reg3 = regtest(res3, model="lm")
rnk3 = ranktest(res3)
mtext("C. R0 and continuity correction for only studies with 0 ", side = 3, adj = 0)
text(-1,0.05, sprintf("Regression test for funnel plot asymmetry: t = %.3f, p = %.3f", reg3$zval, reg3$pval))
text(-1,0.2, sprintf("Rank test for funnel plot asymmetry: t = %.3f, p = %.3f", rnk3$tau, rnk3$pval))

res4 = rma(yi, vi, data=yvi2, method="ML")
funnel(trimfill(res4,estimator="R0"), xlim = c(-6,4))
abline(v=res4$beta)
reg4 = regtest(res4, model="lm")
rnk4 = ranktest(res4)
mtext("D. R0 and continuity correction for all studies", side = 3, adj = 0)
text(-1,0.05, sprintf("Regression test for funnel plot asymmetry: t = %.3f, p = %.3f", reg4$zval, reg4$pval))
text(-1,0.2, sprintf("Rank test for funnel plot asymmetry: t = %.3f, p = %.3f", rnk4$tau, rnk4$pval))

par(mfrow=c(1,1))
dev.off()


## meta-analysis of proportions
#' Derive continuous outcomes (yi and vi)
yvi1 = escalc(measure="PLO", xi=y1, ni=n1, data=data, to="only0")
yi1 = yvi1[,1]
vi1 = yvi1[,2]

yvi2 = escalc(measure="PLO", xi=y1, ni=n1, data=data, to="all")
yi2 = yvi2[,1]
vi2 = yvi2[,2]

# Funnel plot
postscript("funnel2.eps", width = 12, height = 12)
par(mfcol=c(2,2))

res1 = rma(yi1, vi1, data=yvi1)
funnel(trimfill(res1,estimator="L0"))
abline(v=res1$beta)
reg1 = regtest(res1, model="lm")
rnk1 = ranktest(res1)
mtext("A. L0 and continuity correction for only studies with 0 ", side = 3, adj = 0)
text(-3.5,0.05, sprintf("Regression test for funnel plot asymmetry: t = %.3f, p = %.3f", reg1$zval, reg1$pval))
text(-3.5,0.2, sprintf("Rank test for funnel plot asymmetry: t = %.3f, p = %.3f", rnk1$tau, rnk1$pval))

res2 = rma(yi2, vi2, data=yvi2)
funnel(trimfill(res2,estimator="L0"))
abline(v=res2$beta)
reg2 = regtest(res2, model="lm")
rnk2 = ranktest(res2)
mtext("B. L0 and continuity correction for all studies", side = 3, adj = 0)
text(-3.5,0.05, sprintf("Regression test for funnel plot asymmetry: t = %.3f, p = %.3f", reg2$zval, reg2$pval))
text(-3.5,0.2, sprintf("Rank test for funnel plot asymmetry: t = %.3f, p = %.3f", rnk2$tau, rnk2$pval))


res3 = rma(yi1, vi1, data=yvi1)
funnel(trimfill(res3,estimator="R0"))
abline(v=res3$beta)
reg3 = regtest(res3, model="lm")
rnk3 = ranktest(res3)
mtext("C. R0 and continuity correction for only studies with 0 ", side = 3, adj = 0)
text(-3.5,0.05, sprintf("Regression test for funnel plot asymmetry: t = %.3f, p = %.3f", reg3$zval, reg3$pval))
text(-3.5,0.2, sprintf("Rank test for funnel plot asymmetry: t = %.3f, p = %.3f", rnk3$tau, rnk3$pval))

res4 = rma(yi2, vi2, data=yvi2)
funnel(trimfill(res4,estimator="R0"))
abline(v=res4$beta)
reg4 = regtest(res4, model="lm")
rnk4 = ranktest(res4)
mtext("D. R0 and continuity correction for all studies", side = 3, adj = 0)
text(-3.5,0.05, sprintf("Regression test for funnel plot asymmetry: t = %.3f, p = %.3f", reg4$zval, reg4$pval))
text(-3.5,0.2, sprintf("Rank test for funnel plot asymmetry: t = %.3f, p = %.3f", rnk4$tau, rnk4$pval))

par(mfrow=c(1,1))
dev.off()