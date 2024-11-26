#!/usr/bin/Rscript
#  change_coords.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 11.21.2024
source("R/vornoi_cands.R")

set.seed(123)

P <- 2
N <- 15
X <- matrix(runif(N*P, min = 0.1, max = 0.9), ncol=P)

cind <- 1
ns <- cind
Xb <- X[ns,]

norm = 'l2'
Nc <- 50000
#Nc <- 500
Xcand <- vorwalkcands(X, ncand = Nc, st = 'unif', norm = norm, half2bound = F, in_iters = 20)$Xs

im <- which.min(rowSums(t(t(Xcand)-c(0.4,0.9))^2))
xc <- Xcand[im,]

#im <- which.min(rowSums(t(t(Xcand)-c(Xb[1],))^2))
cost <- function(x) {
    1000*abs(x[2]-Xb[2]) + abs(x[1]-0.45)
}
#im2 <- which.min(rowSums(t(t(Xcand)-c(Xb[1],))^2))
im2 <- which.min(c(apply(Xcand, 1, cost)))
xc2 <- Xcand[im2,]

#uscale <- 2.5
uscale <- sqrt(sum((xc-Xb)^2))
U <- matrix(c(1,1)/sqrt(2),nrow=1)*uscale
arlen <- 1.

#uscale <- 2.5
uscale2 <- sqrt(sum((xc2-Xb)^2))
U2 <- matrix(c(1,0),nrow=1)*uscale
arlen <- 1.

## Change of coords figure.

pdf("figures/coords.pdf", width = 4, height = 4)
par(mar=c(1,1,1,0)+0.1)
par(mgp=c(0,0,0)+0.1)
radm <- 0.1
radp <- 0.15
plot(NA,NA,xlim=c(X[cind,1]-radm,X[cind,1]+radp),ylim=c(X[cind,2]-radm,X[cind,2]+radp), main = 'Change of Coordinates', xlab = '', ylab = '', tck = -0.005)
points(Xcand[,1], Xcand[,2], col = 'slategray')
points(Xcand[im,1], Xcand[im,2], col = 'black', lwd = 8, pch = 4)
text(Xcand[im,1]-0.06, Xcand[im,2], 'x=(0.4,0.9)', col = 'black', cex = 1.5)

arrow_end <- Xb+U
arrows(Xb[1], Xb[2], arrow_end[1], arrow_end[2], col = 'blue', lwd = 2)

arrow_end2 <- Xb+U2
#arrows(Xb[1], Xb[2], arrow_end2[1], arrow_end2[2], col = 'gray', lwd = 1, length=0)
lines(c(Xb[1], arrow_end2[1]), c(Xb[2], arrow_end2[2]), col = 'black', lwd = 2, lty='dashed')

ns <- 10000
circ <- matrix(rnorm(ns*2),ncol=2)
circ <- circ / sqrt(rowSums(circ^2))
isok <- circ[,1] > cos(pi/4) & circ[,2] < sin(pi/4) & circ[,2] > 0
circ <- t(t(circ*0.05) + Xb)
points(circ[isok,1],circ[isok,2], lwd = 0.3, cex = 0.5)

text(Xb[1]+0.03, Xb[2]+0.015, expression(theta), col = 'blue', cex = 2)

points(X[,1], X[,2], col = 'red', pch = 4, lwd = 4)
text(arrow_end[1], arrow_end[2]-0.02, 'u', col = 'blue')
dev.off()

