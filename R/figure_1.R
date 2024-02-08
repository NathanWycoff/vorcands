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
y_raw <- apply(X, 1, f)
mu_y <- mean(y_raw)
sig_y <- sd(y_raw)
y <- (y_raw-mu_y) / sig_y

# Fit GP
dmin=0.01; dmax=0.1
da <- darg(list(mle=TRUE, start = (dmin+dmax)/2, min=dmin, max=dmax), randomLHS(1000, ncol(X)))
gpi <- newGPsep(X, y, d=da$start, g=1e-6, dK=TRUE)
mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)

cnt <<- 0
obj <- function(x, fmin, gpi, noise) - EI(gpi, x, fmin, predGPsep, noise)
check <- function(gpi, x, fmin, noise, tol) { EI(gpi, x, fmin, predGPsep, noise) <= tol }

obj(runif(P), min(y), gpi, FALSE)

noise <- FALSE
multi.start <- 1000
fmin <- min(y)
eps <- sqrt(.Machine$double.eps)
tol <- eps

## Find local optima.
library(lhs)
start <- randomLHS(multi.start, ncol(X))
xnew <- matrix(NA, nrow=nrow(start), ncol=ncol(X)+1)
for(i in 1:nrow(start)) {
    if (check(gpi, start[i,], fmin, noise, tol)) {
        out <- list(value=Inf) 
    }
    else {
        out <- optim(start[i,], obj, method="L-BFGS-B", lower=0, upper=1, gpi=gpi, fmin=fmin, noise=noise) 
        #print(out$message)
        if (substr(out$message,1,11)!='CONVERGENCE') {
            out$value <- Inf
        }
    }
    xnew[i,] <- c(out$par, -out$value)
}

## Create grid for EI surface viz.
ng <- 100
N <- ng*ng
P <- 2
Xg <- matrix(NA, nrow=N, ncol=P)
xgrid <- seq(0,1,length.out=ng)
Yg <- matrix(NA,nrow=ng,ncol=ng)
for (n1 in 1:ng) {
    for (n2 in 1:ng) {
        ind <- (n1-1)*ng+n2
        Xg[ind,1] <- xgrid[n1]
        Xg[ind,2] <- xgrid[n2]
        Yg[n1,n2] <- obj(Xg[ind,], fmin, gpi, noise)
    }
}
Yg <- -log10(-Yg+1e-4)
vals <- xnew[,P+1]

# Drop bad optima.
newtol <- 5e-8

lm <- xnew[is.finite(xnew[,P+1]) & vals>newtol,1:P,drop=F]

Nx <- 100000
norm <- 'linf'
Xcand <- vorwalkcands(X, ncand = Nx, st = 'unif', norm =norm, half2bound = FALSE)$Xs

#### 
## Make toy acquisition near other point. 
ll_quad <- Xcand[,1]<0.4 & Xcand[,1]>0.3 & Xcand[,2] < 0.3 & Xcand[,2]>0.1
CE <- Xcand[ll_quad,]
ei_ce <- rep(NA,nrow(CE))
for (i in 1:nrow(CE)) {
    ei_ce[i] <- obj(CE[i,], fmin, gpi, noise)
}
it_min <- which.min(ei_ce)

CE[it_min,]

Xnew <- rbind(X, CE[it_min,])

## Get some candidates on the boundary, not evenly distributed.
#Nx2 <- 100000
Nx2 <- 100000
ret <- vorwalkcands(Xnew, ncand = Nx2, st = 'unif', norm =norm, half2bound = FALSE, lb = 0, ub = 0.5, in_iters = 20)
dtn <- apply(ret$Xs, 1, function(x) max(abs(x-Xnew[nrow(Xnew),])))
dt2 <- apply(ret$Xs, 1, function(x) max(abs(x-Xnew[2,])))
ours <- (ret$ns == nrow(Xnew)) | (ret$ns == 2 & (abs(dtn-dt2)<1e-4))
Xcand2 <- ret$Xs[ours,]

## Coarsen to the largest distance to evenly distribute.
nn <- RANN.Linf::nn2(Xcand2, Xcand2, 2)
mnd <- max(nn$nn.dists[,2])
XC <- Xcand2
ndropped <- 0
for (i in 1:nrow(Xcand2)) {
    ind <- i-ndropped
    nni <- RANN.Linf::nn2(XC,XC[ind,,drop=F],2)
    ndi <- nni$nn.dists[,2]
    if (ndi<mnd) {
        XC <- XC[-ind,]
        ndropped <- ndropped+1
    }
}
Xcand2 <- XC

## Find corners 
nn <- RANN.Linf::nn2(Xcand2, Xcand2, 3)$nn.idx[,2:3]
sv2s <- rep(NA, nrow(Xcand2))
for (i in 1:nrow(Xcand2)) {
    seg1 <- Xcand2[i,] - Xcand2[nn[i,1],]
    seg2 <- Xcand2[i,] - Xcand2[nn[i,2],]
    seg1 <- seg1 / sqrt(sum(seg1^2))
    seg2 <- seg2 / sqrt(sum(seg2^2))
    sv2s[i] <- svd(rbind(seg1,seg2))$d[2]
}
#corners <- which(rank(-sv2s)<=11)
corners <- which(rank(-sv2s)<=20)

## Now coarsen the overall candidates just a little bit to help with the pdf size.
nn <- RANN.Linf::nn2(Xcand, Xcand, 2)
#mnd <- max(nn$nn.dists[,2])
todrop = 0.8*nrow(Xcand)
XC <- Xcand
#at_a_time <- 1000
at_a_time <- 100
iters <- ceiling(todrop/at_a_time)
for (i in 1:iters) {
    nni <- RANN.Linf::nn2(XC,XC,2)
    ndi <- nni$nn.dists[,2]
    #XC <- XC[-which.min(ndi),]
    XC <- XC[-which(rank(ndi)<=at_a_time),]
}
Xcand_cut <- XC

## Order corners by angle made with central point.
Xcorn <- Xcand2[corners,]
xo <- CE[it_min,]
XD <- t(t(Xcorn) - xo)
angle <- apply(XD, 1, function(x) atan2(x[1],x[2]))
ord <- order(angle)
Xcorn <- Xcorn[ord,]

lm_u <- lm
D <- as.matrix(dist(lm_u))

todrop <- c()
lm_thresh <- 1e-4
for (i in 1:(nrow(lm_u)-1)) {
    if (any(D[i,1:(i-1)] < lm_thresh)) {
        todrop <- c(todrop, i)
    }
}
lm_u <- lm_u[-todrop,]

## Plot everything, finally.
pdf('expo_it3.pdf', width=5, height = 5)
par(mar=c(1.5,1.5,0.4,0.4)+0.1)
par(mgp=c(1.5,0.5,0))

image(Yg, xlab = '', ylab = '')
points(X[,1],X[,2], pch = 4, lwd = 8, col = 'white')
points(X[,1],X[,2], pch = 4, lwd = 6)

for (i in 1:length(corners)) {
    if (i<length(corners)) {
        ind <- i:(i+1)
    } else {
        ind <- c(length(corners),1)
    }
    lines(Xcorn[ind,1], Xcorn[ind,2], lty = 2, col = 'slategray', lwd = 4)
}
points(Xcand_cut[,1], Xcand_cut[,2], col = 'slategray', cex=0.4)
lines(c(0,0),c(0,1), col = 'slategray', lwd = 4)
lines(c(0,1),c(1,1), col = 'slategray', lwd = 4)
lines(c(1,1),c(1,0), col = 'slategray', lwd = 4)
lines(c(1,0),c(0,0), col = 'slategray', lwd = 4)

#points(Xcand2[,1], Xcand2[,2], col = 'brown', cex=0.4)

points(CE[it_min,1], CE[it_min,2], pch = 3, lwd = 10, col = 'white')
points(CE[it_min,1], CE[it_min,2], pch = 3, lwd = 8, col = 'blue')
points(lm_u[,1], lm_u[,2], col = 'white', pch = 5, lwd = 7)
points(lm_u[,1], lm_u[,2], col = 'deeppink', pch = 5, lwd = 5)

legend('topright', legend = c("Design Points","EI Local Optima"), col = c("Black","deeppink"), 
    pch = c(4,5), pt.lwd = c(6,5), bg = 'white')

dev.off()

## Obtained from https://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)

    ticks_lab <- round(ticks, 1)

    ylab = substitute(paste(bold("log Expected Improvement")))
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab=ylab, main=title)
    axis(2, ticks, ticks_lab, las=1, cex.axis = 0.75)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
}

#pdf("cb.pdf",width=1.35, height=5)
#pdf("cb.pdf",width=1.35*0.6, height=5)
pdf("cb.pdf",width=0.875, height=5)
par(mar=c(0,2.5,0,0)+0.1, mgp = c(1.6,0.6,0))
color.bar(hcl.colors(100, "YlOrRd"), min(-Yg), max(-Yg), nticks = 7)
dev.off()

5 / (0.37/0.1)
5 / (0.40/0.07)

5 / (x/y) = 0.81
1 / (x/y) = 0.81 / 5
(x/y) = 1/(0.81 / 5)
(x/y) = 6.17284
x+y = 0.47
#x = 0.4
