#'
#' Meta-analysis of OR without PB
#'

rm(list=ls())
file.sources = list.files("Rfn/")
sapply(paste0("Rfn/", file.sources), source)
library(metafor)
library(metadat)


## Data
data = read.csv("niel-weise21.csv")
# data = read.csv("thomas.csv")

## Meta-analysis of ORs ----
# Derive continuous outcomes (lnOR and se)
yvi1 = escalc(measure="OR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0, data=data, to="only0")
yi1 = yvi1$yi
vi1 = yvi1$vi

yvi2 = escalc(measure="OR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0, data=data, to="all")
yi2 = yvi2$yi
vi2 = yvi2$yi

#' NN model
lnOR_nn = NN_LMM(
  yi=yi1, vi=vi1, 
  parset=list(
    mu.bound = 10, 
    tau.bound = 5,
    eps = 1e-3,
    init.vals = NULL
  ))
lnOR_nn_lb = lnOR_nn$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*lnOR_nn$mu[2]
lnOR_nn_ub = lnOR_nn$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*lnOR_nn$mu[2]
sprintf("NN: theta (95CI, SE): %.3f (%.3f, %.3f; %.3f)", 
        lnOR_nn$mu[1], lnOR_nn_lb, lnOR_nn_ub, lnOR_nn$mu[2])
sprintf("NN: tau (SE): %.3f (%.3f)", 
        lnOR_nn$tau[1], lnOR_nn$tau[2])

# rma(yi1,vi1, method = "ML")

#' HN-GLMM (takes longer time)
lnOR_hn = HN_GLMM(
  y0 = data$y0, y1 = data$y1, n0 = data$n0, n1 = data$n1, 
  parset = list(
    mu.bound = 10,
    tau.bound = 5,
    eps = 1e-3,
    integ.limit = 10, 
    cub.tol = 1e-5,
    init.vals = NULL))
lnOR_hn_lb = lnOR_hn$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*lnOR_hn$mu[2]
lnOR_hn_ub = lnOR_hn$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*lnOR_hn$mu[2]
sprintf("HN: theta (95CI, SE): %.3f (%.3f, %.3f; %.3f)", 
        lnOR_hn$mu[1], lnOR_hn_lb, lnOR_hn_ub, lnOR_hn$mu[2])
sprintf("HN: tau (SE): %.3f (%.3f)", 
        lnOR_hn$tau[1], lnOR_hn$tau[2])

# rma.glmm(measure = "OR", ai=y1, bi=n1-y1,ci=y0,di=n0-y0, data = data, model = "CM.EL")

#' BN-GLMM
lnOR_bn = BN_GLMM(
  y0 = data$y0, y1 = data$y1, n0 = data$n0, n1 = data$n1, 
  parset = list(
    mu.bound = 10,
    tau.bound = 5,
    eps = 1e-3,
    integ.limit = 10, 
    cub.tol = 1e-5,
    init.vals = NULL))
lnOR_bn_lb = lnOR_bn$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*lnOR_bn$mu[2]
lnOR_bn_ub = lnOR_bn$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*lnOR_bn$mu[2]
sprintf("BN: theta (95CI, SE): %.3f (%.3f, %.3f; %.3f)", 
        lnOR_bn$mu[1], lnOR_bn_lb, lnOR_bn_ub, lnOR_bn$mu[2])
sprintf("BN: tau (SE): %.3f (%.3f)", 
        lnOR_bn$tau[1], lnOR_bn$tau[2])

# rma.glmm(measure = "OR", ai=y1, bi=n1-y1,ci=y0,di=n0-y0, data = data, model = "CM.AL")

# save(lnOR_nn,lnOR_hn,lnOR_bn, file = "res/app1-nopb.RData")
# save(lnOR_nn,lnOR_hn,lnOR_bn, file = "res/app2-nopb.RData")