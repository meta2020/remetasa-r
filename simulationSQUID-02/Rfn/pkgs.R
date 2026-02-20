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

set.cutoff2 = c(0.01,0.03,0.05,0.1)
set.wi2 = c(0.9, 0.8, 0.7, 0.5, 0.5)