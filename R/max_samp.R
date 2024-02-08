## This makes the figure with various candidate schemes.

library(colorRamps)
library(SobolSequence)
library(lhs)
library(laGP)

source("tricands/R/tricands.R")
source("R/vornoi_cands.R")

set.seed(123)

func <- 'ackley2'
source("R/sim_settings.R")

f <- function(x) sum(sin(2*2*pi*(x-0.5))^2)

source("R/optim.R")

f2 <- 10
N <- 32 # good but lengthscales from GP are wrong.

P <- 2

X <- randomLHS(N, P)

Nx <- 10000

norm <- 'linf'
Xcand_std <- vorwalkcands(X, ncand = Nx, st = 'unif', norm =norm, half2bound = FALSE)$Xs

pdf("std.pdf")
plot(NA,NA,xlim=c(0,1),ylim=c(0,1))
points(X[,1],X[,2])
points(Xcand_std[,1],Xcand_std[,2])
dev.off()

n <- 1
ns <- rep(n,4)

U <- matrix(0,2*1*P,P)
for (p in 1:P) {
    up <- (1:nrow(U)-1) %% (2*P) == 2*(p-1)
    down <- (1:nrow(U)-1) %% (2*P) == 2*(p-1)+1
    U[up,p] <- 1
    U[down,p] <- -1
}

Xcand_max <- vorwalk(X, ns, U, norm = norm, in_iters = 10, half2bound = FALSE)$Xs
#duplicated(Xcand_max)


N_rem <- Nx - ncand_init
nblocks <- ceiling(N_rem / N)
#for (i in 1:nblocks) {
for (i in 1:N_rem) {
    dists <- RANN.Linf::nn2(Xcand_max,Xcand_max,2)
    maxb <- which.max(dists$nn.dists[,2])
    j1 <- maxb
    j2 <- dists$nn.idx[maxb,2]
    mid <- (Xcand_max[j1,,drop=F] + Xcand_max[j2,drop=F])/2

    u <- mid - X[n,,drop=F] 
    u <- sqrt(P) * u / sqrt(sum(u^2))

    xnew <- vorwalk(X, n, u, norm = norm, in_iters = 10, half2bound = FALSE)$Xs
    Xcand_max <- rbind(Xcand_max, xnew)
}

pdf("max.pdf")
plot(NA,NA,xlim=c(0,1),ylim=c(0,1))
points(X[,1],X[,2], pch = 4, col = 'red')
points(Xcand_max[,1],Xcand_max[,2])
#points(Xcand_max[j1,1],Xcand_max[j1,2], col = 'cyan')
#points(Xcand_max[j2,1],Xcand_max[j2,2], col = 'cyan')
#points(mid[1], mid[2], col = 'orange')
#points(xnew[1], xnew[2], col = 'blue')
dev.off()

