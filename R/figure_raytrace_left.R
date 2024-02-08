source("R/vornoi_cands.R")

set.seed(123)

P <- 2
N <- 15
X <- matrix(runif(N*P, min = 0.1, max = 0.9), ncol=P)

norm = 'l2'
Nc <- 50000
Xcand <- vorwalkcands(X, ncand = Nc, st = 'unif', norm = norm, half2bound = F, in_iters = 20)$Xs

### Left Figure
cind <- 1
rad <- 0.15
X[cind,]
ns <- cind

#uscale <- 2.5
uscale <- 3.2
U <- matrix(c(1,1)/sqrt(2),nrow=1)/uscale

## fragment from vornoi source.
ncand <- 1
in_iters <- 10
Xb <- X[ns,]
ts <- rep(1,ncand)

lb <- rep(0,ncand)
ub <- rep(1,ncand)

T_trace <- matrix(NA, in_iters, ncand)
for (it in 1:in_iters) {
    mid <- (lb + ub) / 2
    Xs <- getXs(mid, Xb, U)

    # Condition 1 for positive sign (true): must be closer to anchor point than any other point in design.
    fmid <- getsign(X, Xs, ns, norm = norm)
    dp <- mean(fmid$sign)
    in_upper <- fmid$sign
    nd1 <- sum(!in_upper)

    # Condition 2 must be closer to anchor point than to boundary/exterior of the hybercube.
    outofbounds <- apply(abs(Xs-0.5)>0.5, 1, any)
    in_upper <- in_upper & (!outofbounds)

    lb[in_upper] <- mid[in_upper]
    ub[!in_upper] <- mid[!in_upper]
    T_trace[it,] <- mid
}
## end fragment

arlen <- 0.05

png("figures/bisect.png", width = 4, height = 4, unit = 'in', res = 1200)
par(mar=c(1,1,1,0)+0.1)
par(mgp=c(0,0,0)+0.1)
radm <- 0.1
radp <- 0.15
plot(NA,NA,xlim=c(X[cind,1]-radm,X[cind,1]+radp),ylim=c(X[cind,2]-radm,X[cind,2]+radp), main = 'Bisection Search', xlab = '', ylab = '', tck = -0.005)
#plot(NA,NA,xlim=c(X[cind,1]-rad,X[cind,1]+rad),ylim=c(X[cind,2]-rad,X[cind,2]+rad), main = 'Bisection Search', xlab = '', ylab = '', tck = -0.005)
points(Xcand[,1], Xcand[,2], col = 'slategray')
arrow_end <- Xb+c(arlen*U*uscale)
#arrows(Xb[1],Xb[2], Xb[1]+arlen*U[,1],Xb[2]+arlen*U[,2])
arrows(Xb[1], Xb[2], arrow_end[1], arrow_end[2], col = 'blue', lwd = 2)
points(X[,1], X[,2], col = 'red', pch = 4, lwd = 4)
for (it in 1:in_iters) {
    xit <- getXs(T_trace[it,],Xb,U)
    points(xit[1], xit[2], pch = 3)
    if (it <= 4) {
        text(xit[1]+0.006, xit[2]-0.006, it)
    }
}
text(arrow_end[1], arrow_end[2]-0.02, 'u', col = 'blue')
dev.off()

### Mid and right
for (st in c('unif','rect')) {
    nvizcand <- 100
    ret <- vorwalkcands(X, ncand = nvizcand, st = st, norm = norm, half2bound = F, in_iters = 20)
    Xcand2 <- ret$Xs
    ns <- ret$ns
    nvizcand <- length(ns)

    cind <- 10

    if (st=='unif') {
        title <- "Uniform Angles"
    } else if (st=='rect') {
        title <- "Axis-Aligned Angles"
    } else stop("Unknown st")

    png(paste("figures/vor_",st,".png", sep=''), width = 4, height = 4, unit = 'in', res = 1200)
    par(mar=c(1,1,1,0)+0.1)
    par(mgp=c(0,0,0)+0.1)
    radp <- 0.29
    radm <- radp
    plot(NA,NA,xlim=c(X[cind,1]-radm,X[cind,1]+radp),ylim=c(X[cind,2]-radm,X[cind,2]+radp), main = title, xlab = '', ylab = '', tck = -0.005)
    for (i in 1:nvizcand) {
        #points(X[ns[i],1],X[ns[i],2],Xcand2[i,1],Xcand2[i,2], type = 'l')
        #abline(a = c(X[ns[i],1],X[ns[i],2]), b = c(Xcand2[i,1],Xcand2[i,2]))
        points(c(X[ns[i],1],Xcand2[i,1]), c(X[ns[i],2],Xcand2[i,2]), type = 'l', lty = 'dotted')
    }
    points(Xcand[,1], Xcand[,2], col = 'slategray')
    points(Xcand2[,1], Xcand2[,2], col = 'black', pch = 5, lwd = 4)
    points(Xcand2[,1], Xcand2[,2], col = 'white', pch = 5, lwd = 3)
    points(Xcand2[,1], Xcand2[,2], col = 'cyan', pch = 5, lwd = 2)
    points(X[,1], X[,2], col = 'red', pch = 4, lwd = 4)
    dev.off()
}
