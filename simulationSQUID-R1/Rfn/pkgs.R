library(dplyr)
library(MASS)
library(mnormt)
library(foreach)
library(doRNG)
library(metafor)

s = c(15,50)

set.eps = 1e-4
set.int.lmt = 10
set.cub.tol = 1e-4
set.tau.bound = 2
set.mu.bound = abs(-2)*3

rtimes=1000


