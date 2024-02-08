#!/usr/bin/Rscript
#  tri_vs_vor_timing.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 02.01.2024

source("R/vornoi_cands.R")
source("tricands/R/tricands.R")

set.seed(123)

P <- 10
N <- 100
X <- matrix(runif(N*P, min = 0.1, max = 0.9), ncol=P)

norm = 'l2'

Nc <- 2000
tt <- Sys.time()
Xcand <- vorwalkcands(X, ncand = Nc, st = 'unif', norm = norm, half2bound = F, in_iters = 20)$Xs
print(Sys.time() - tt)

tt <- Sys.time()
Xcand <- tricands(X, max = Nc)
print(Sys.time() - tt)
