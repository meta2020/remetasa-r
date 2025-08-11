#'
#' Meta-analysis of without PB
#'
#' Load R functions
rm(list=ls())
file.sources = list.files("Rfn/")
sapply(paste0("Rfn/", file.sources), source)
library(metafor)
library(metadat)


## Data
# data = read.csv("niel-weise21.csv")[,c(1,4,5)]
data = read.csv("pritz1997.csv")

## Meta-analysis of proportions ----
#' Derive continuous outcomes (logit-proportion and se)
yvi1 = metafor::escalc(measure="PLO", xi=y1, ni=n1, data=data, to = "only0")
yi1 = yvi1$yi
vi1 = yvi1$vi

yvi2 = metafor::escalc(measure="PLO", xi=y1, ni=n1, data=data, to = "all")
yi2 = yvi2$yi
vi2 = yvi2$vi


#' NN model
lgtP_nn = NN_LMM(
  yi=yi1, vi=vi1, 
  parset=list(
    mu.bound = 10, 
    tau.bound = 5,
    eps = 1e-3,
    init.vals = NULL
  ))
lgtP_nn_lb = lgtP_nn$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*lgtP_nn$mu[2]
lgtP_nn_ub = lgtP_nn$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*lgtP_nn$mu[2]
sprintf("NNprop: theta (95CI, SE): %.3f (%.3f, %.3f; %.3f)", 
        lgtP_nn$mu[1], lgtP_nn_lb, lgtP_nn_ub, lgtP_nn$mu[2])
sprintf("NNprop: tau (SE): %.3f (%.3f)", 
        lgtP_nn$tau[1], lgtP_nn$tau[2])

#' BN-GLMM
lgtP_bn = BN_GLMM_prop(
  y1 = data$y1, n1 = data$n1, 
  parset = list(
    mu.bound = 10,
    tau.bound = 5,
    eps = 1e-3,
    integ.limit = 10, 
    cub.tol = 1e-5,
    init.vals = NULL))
lgtP_bn_lb = lgtP_bn$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*lgtP_bn$mu[2]
lgtP_bn_ub = lgtP_bn$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*lgtP_bn$mu[2]
sprintf("1GBN: theta (95CI, SE): %.3f (%.3f, %.3f; %.3f)", 
        lgtP_bn$mu[1], lgtP_bn_lb, lgtP_bn_ub, lgtP_bn$mu[2])
sprintf("1GBN: tau (SE): %.3f (%.3f)", 
        lgtP_bn$tau[1], lgtP_bn$tau[2])


# save(lgtP_nn,lgtP_bn, file = "res/app3-nopb.RData")
