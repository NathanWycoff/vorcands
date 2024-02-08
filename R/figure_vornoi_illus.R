library(latex2exp)

source("R/vornoi_cands.R")

set.seed(123)

P <- 2
N <- 15
X <- matrix(runif(N*P, min = 0.1, max = 0.9), ncol=P)

for (norm in c("l1","l2","linf")) {
    Nc <- 50000
    Xcand <- vorwalkcands(X, ncand = Nc, st = 'unif', norm = norm, half2bound = F)$Xs

    if (norm=='l1') {
        tit <- TeX(r'(1-norm)')
    } else if (norm=='l2') {
        tit <- TeX(r'(2-norm)')
    } else if (norm=='linf') {
        tit <- TeX(r'($\infty$-norm)')
    } else {
        stop("Unknown norm.")
    }

    #fn <- paste('figures/vorill_',norm,'.pdf',sep='')
    fn <- paste('figures/vorill_',norm,'.png',sep='')

    #pdf(fn, width=4, height=4)
    png(fn, width=4, height=4, units = 'in', res = 1200)
    par(mar=c(1,1,1,0)+0.1)
    par(mgp=c(0,0,0)+0.1)
    plot(NA,NA,xlim=c(0,1),ylim=c(0,1), main = tit, xlab = '', ylab = '', tck = -0.005)
    points(Xcand[,1],Xcand[,2], col = 'slategray')
    points(X[,1],X[,2], col = 'red', pch = 4, lwd = 4)
    dev.off()
}
