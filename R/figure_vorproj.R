#!/usr/bin/Rscript
#  figure_vorproj.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 02.01.2024

source("R/vornoi_cands.R")

set.seed(123)

P <- 2
N <- 15
X <- matrix(runif(N*P, min = 0.1, max = 0.9), ncol=P)

norm = 'l2'
#Nc <- 500
Nc <- 50000
Xcand <- vorwalkcands(X, ncand = Nc, st = 'unif', norm = norm, half2bound = F, in_iters = 20)$Xs

# Single proj
cind <- 1
rad <- 0.15
X[cind,]
ns <- cind

u <- c(-0.04, -0.02)
lhs_point <- X[cind,] + u
U <- matrix(u, nrow = 1)
U <- U / sqrt(sum(U^2))

lhs_proj <- vorwalk(X, cind, U, norm = 'l2', in_iters = 20)$Xs

arlen <- 0.05

png("figures/single_proj.png", width = 4, height = 4, unit = 'in', res = 1200)
par(mar=c(1,1,1,0)+0.1)
par(mgp=c(0,0,0)+0.1)
radm <- 0.15
radp <- 0.1
plot(NA,NA,xlim=c(X[cind,1]-radm,X[cind,1]+radp),ylim=c(X[cind,2]-radm,X[cind,2]+radp), main = 'Point Projection', xlab = '', ylab = '', tck = -0.005)
#plot(NA,NA,xlim=c(X[cind,1]-rad,X[cind,1]+rad),ylim=c(X[cind,2]-rad,X[cind,2]+rad), main = 'Bisection Search', xlab = '', ylab = '', tck = -0.005)
arrows(X[cind,1], X[cind,2], lhs_proj[1], lhs_proj[2], col = 'blue', lwd = 2)
points(Xcand[,1], Xcand[,2], col = 'slategray')
points(X[,1], X[,2], col = 'red', pch = 4, lwd = 6)
points(lhs_point[1], lhs_point[2], col = 'orange', pch = 2, lwd = 7.5)
points(lhs_proj[1], lhs_proj[2], col = 'black', pch = 5, lwd = 6+6)
points(lhs_proj[1], lhs_proj[2], col = 'white', pch = 5, lwd = 6+6/2)
points(lhs_proj[1], lhs_proj[2], col = 'cyan', pch = 5, lwd = 6)
dev.off()


# Multi proj
cind <- 10
rad <- 0.15
X[cind,]
ns <- cind

Nc <- 100

Xlhs <- lhs::randomLHS(Nc,P)
knn <- RANN::nn2(X, Xlhs, k = 1) # NOTE: Assumes l2 norm!
ns <- c(knn$nn.idx)
U <- matrix(NA, nrow = Nc, ncol = P)
for (i in 1:Nc) {
    U[i,] <- Xlhs[i,] - X[ns[i],]
}
U <- U / sqrt(rowSums(U^2))

lhs_proj <- vorwalk(X, ns, U, norm = 'l2', in_iters = 20, half2bound = F)$Xs

arlen <- 0.05

png("figures/multi_proj.png", width = 4, height = 4, unit = 'in', res = 1200)
par(mar=c(1,1,1,0)+0.1)
par(mgp=c(0,0,0)+0.1)
radm <- 0.3
radp <- radm
plot(NA,NA,xlim=c(X[cind,1]-radm,X[cind,1]+radp),ylim=c(X[cind,2]-radm,X[cind,2]+radp), main = 'Proj of LHS', xlab = '', ylab = '', tck = -0.005)
#plot(NA,NA,xlim=c(X[cind,1]-rad,X[cind,1]+rad),ylim=c(X[cind,2]-rad,X[cind,2]+rad), main = 'Bisection Search', xlab = '', ylab = '', tck = -0.005)
#arrows(X[cind,1], X[cind,2], lhs_proj[1], lhs_proj[2], col = 'blue', lwd = 2)
for (i in 1:Nc) {
    points(c(X[ns[i],1],lhs_proj[i,1]), c(X[ns[i],2],lhs_proj[i,2]), type = 'l', lty = 'dotted')
}
points(Xcand[,1], Xcand[,2], col = 'slategray')
points(X[,1], X[,2], col = 'red', pch = 4, lwd = 3)
points(Xlhs[,1], Xlhs[,2], col = 'orange', pch = 2, lwd = 2)
points(lhs_proj[,1], lhs_proj[,2], col = 'black', pch = 5, lwd = 4)
points(lhs_proj[,1], lhs_proj[,2], col = 'white', pch = 5, lwd = 3)
points(lhs_proj[,1], lhs_proj[,2], col = 'cyan', pch = 5, lwd = 2)
dev.off()

