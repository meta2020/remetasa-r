#'
#' Meta-analysis of the logit-proportion (Example 2)
#'
#' Load R functions
rm(list=ls())
file.sources = list.files("Rfn/")
sapply(paste0("Rfn/", file.sources), source)


# data=dat.pritz1997[,-2]
# colnames(data)=c("study","y0","n1")
# data$y1=data$n1-data$y0

data = dat.axfors2021[,c(1,7,9,10)]
data = data[data$Published=="Published",c(1,4,3)]
colnames(data)=c("study","y1","n1")


yvi1 = metafor::escalc(measure="PLO", xi=y1, ni=n1, data=data, to = "only0")
yi1 = yvi1$yi
vi1 = yvi1$vi

yvi2 = metafor::escalc(measure="PLO", xi=y1, ni=n1, data=data, to = "all")
yi2 = yvi2$yi
vi2 = yvi2$vi

#' Proposed sensitivity analysis for PB in meta-analysis ----------
#' (Pnmax = 0.999, Pnmin = p_sa)
#' (Psemax = p_sa, Psemin = 0.999)
#' 
#' Proposal for BN-GLMM
p_sa = c(0.99, seq(0.9, 0.1, -0.1))
lgtP_COPAS_BNGLMM = vapply(
  p_sa, 
  function(p) {
    mod = COPAS_BNGLMM_prop(
      y1=data$y1, n1=data$n1, Pnmax = 0.999, Pnmin = p, 
      parset = list(
        mu.bound = 10,
        tau.bound = 5,
        estimate.rho = TRUE, 
        eps = 1e-3,
        integ.limit = 10, 
        cub.tol = 1e-5,
        init.vals = NULL)
    )
    mu_lb = mod$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*mod$mu[2]
    mu_ub = mod$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*mod$mu[2]
    c(mod$mu, mu_lb, mu_ub, mod$tau[1:2], mod$rho, mod$a, mod$opt$convergence)
  }, 
  c("mu"=0, "mu.se"=0, "mu.lb"=0, "mu.ub"=0, 
    "tau"=0, "tau.se"=0,"rho"=0, "rho.se"=0, 
    "a0"=0, "a1"=0, "converge"=0))
colnames(lgtP_COPAS_BNGLMM) = paste0("p = ", p_sa)
lgtP_COPAS_BNGLMM

#' Copas-Heckman-type selection model (2000)
lgtP_COPAS2000_1 = vapply(
  p_sa, 
  function(p) {
    mod = COPAS2000(
      yi = yi1, vi = vi1, Psemax = p, Psemin = 0.999, 
      parset = list(
        mu.bound = 10,
        tau.bound = 5,
        estimate.rho = TRUE, 
        eps = 1e-3,
        init.vals = NULL))
    mu_lb = mod$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*mod$mu[2]
    mu_ub = mod$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*mod$mu[2]
    c(mod$mu, mu_lb, mu_ub, mod$tau[1:2], mod$rho, mod$gamma, mod$opt$convergence)
  }, 
  c("mu"=0, "mu.se"=0, "mu.lb"=0, "mu.ub"=0, 
    "tau"=0, "tau.se"=0,"rho"=0, "rho.se"=0, 
    "gamma0"=0, "gamma1"=0, "converge"=0))
colnames(lgtP_COPAS2000_1) = paste0("p = ", p_sa)
lgtP_COPAS2000_1

lgtP_COPAS2000_2 = vapply(
  p_sa, 
  function(p) {
    mod = COPAS2000(
      yi = yi2, vi = vi2, Psemax = p, Psemin = 0.999, 
      parset = list(
        mu.bound = 10,
        tau.bound = 5,
        estimate.rho = TRUE, 
        eps = 1e-3,
        init.vals = NULL))
    mu_lb = mod$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*mod$mu[2]
    mu_ub = mod$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*mod$mu[2]
    c(mod$mu, mu_lb, mu_ub, mod$tau[1:2], mod$rho, mod$gamma, mod$opt$convergence)
  }, 
  c("mu"=0, "mu.se"=0, "mu.lb"=0, "mu.ub"=0, 
    "tau"=0, "tau.se"=0,"rho"=0, "rho.se"=0, 
    "gamma0"=0, "gamma1"=0, "converge"=0))
colnames(lgtP_COPAS2000_2) = paste0("p = ", p_sa)
lgtP_COPAS2000_2

#' Copas-N-type selection model (2000)
lgtP_COPAS1999_1 = vapply(
  p_sa, 
  function(p) {
    mod = COPAS1999(
      yi = yi1, ni = data$n1, Pnmax = 0.999, Pnmin = p, 
      parset = list(
          mu.bound = 10,
          tau.bound = 5,
          estimate.rho = TRUE, 
          eps = 1e-3,
          init.vals = NULL))
    mu_lb = mod$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*mod$mu[2]
    mu_ub = mod$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*mod$mu[2]
    c(mod$mu, mu_lb, mu_ub, mod$tau[1:2], mod$rho, mod$gamma, mod$opt$convergence)
  }, 
  c("mu"=0, "mu.se"=0, "mu.lb"=0, "mu.ub"=0, 
    "tau"=0, "tau.se"=0,"rho"=0, "rho.se"=0, 
    "gamma0"=0, "gamma1"=0, "converge"=0))
colnames(lgtP_COPAS1999_1) = paste0("p = ", p_sa)
lgtP_COPAS1999_1

lgtP_COPAS1999_2 = vapply(
  p_sa, 
  function(p) {
    mod = COPAS1999(
      yi = yi2, ni = data$n1, Pnmax = 0.999, Pnmin = p, 
      parset = list(
        mu.bound = 10,
        tau.bound = 5,
        estimate.rho = TRUE, 
        eps = 1e-3,
        init.vals = NULL))
    mu_lb = mod$mu[1] + qnorm((1-0.95)/2, lower.tail = TRUE)*mod$mu[2]
    mu_ub = mod$mu[1] + qnorm((1-0.95)/2, lower.tail = FALSE)*mod$mu[2]
    c(mod$mu, mu_lb, mu_ub, mod$tau[1:2], mod$rho, mod$gamma, mod$opt$convergence)
  }, 
  c("mu"=0, "mu.se"=0, "mu.lb"=0, "mu.ub"=0, 
    "tau"=0, "tau.se"=0,"rho"=0, "rho.se"=0, 
    "gamma0"=0, "gamma1"=0, "converge"=0))
colnames(lgtP_COPAS1999_2) = paste0("p = ", p_sa)
lgtP_COPAS1999_2

#' NUMBER OF THE UNPUBLISHED
M_n = sapply(p_sa, function(p) {
  
  P_max = 0.999
  P_min = p

  n_min = min(data$n1) 
  n_max = max(data$n1)
  
  a1 = (qnorm(P_max)-qnorm(P_min))/(sqrt(n_max)-sqrt(n_min))
  a0 = qnorm(P_max)-a1*sqrt(n_max)
  sum((1 - pnorm(a0+a1*sqrt(data$n1)))/pnorm(a0+a1*sqrt(data$n1))) 
  
})


M_copas1 = sapply(p_sa, function(p) {
  
  P_max = p
  P_min = 0.999

  se_min_inv = 1/sqrt(min(vi1)) 
  se_max_inv = 1/sqrt(max(vi1))
  
  gamma1 = (qnorm(P_max)-qnorm(P_min))/(se_max_inv-se_min_inv)
  gamma0 = qnorm(P_max)-gamma1*se_max_inv
  sum((1 - pnorm(gamma0+gamma1/sqrt(vi1)))/pnorm(gamma0+gamma1/sqrt(vi1)))
  
})

M_copas2 = sapply(p_sa, function(p) {
  
  P_max = p
  P_min = 0.999
  
  se_min_inv = 1/sqrt(min(vi2)) 
  se_max_inv = 1/sqrt(max(vi2))
  
  gamma1 = (qnorm(P_max)-qnorm(P_min))/(se_max_inv-se_min_inv)
  gamma0 = qnorm(P_max)-gamma1*se_max_inv
  sum((1 - pnorm(gamma0+gamma1/sqrt(vi2)))/pnorm(gamma0+gamma1/sqrt(vi2)))
  
})

#' Table 1
tab1 = data.frame(
  BN = t(lgtP_COPAS_BNGLMM[c(1,3,4),]),
  CH1 = t(lgtP_COPAS2000_1[c(1,3,4),]), # copas-heckman
  CH2 = t(lgtP_COPAS2000_2[c(1,3,4),]), # copas-heckman
  CN1 = t(lgtP_COPAS1999_1[c(1,3,4),]), # copas-n
  CN2 = t(lgtP_COPAS1999_2[c(1,3,4),]), 
  M.c1 = round(M_copas1), 
  M.c2 = round(M_copas2), 
  M.p = round(M_n)
)

# tab1_all = rbind("p = 1" = tab1_p1, tab1)
tab1_all = tab1
tab1_all$pnmin = p_sa
tab1_all$pnmax = rep(0.999, 10)

#' SAVE RESULTS1
#' 
save(lgtP_COPAS_BNGLMM, 
     lgtP_COPAS2000_1,lgtP_COPAS2000_2,
     lgtP_COPAS1999_1,lgtP_COPAS1999_2,
     tab1_all,
     file = "res/app4-prop.RData")


