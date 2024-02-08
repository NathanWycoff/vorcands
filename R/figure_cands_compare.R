## This makes the figure with various candidate schemes.

library(colorRamps)
source("tricands/R/tricands.R")
source("R/vornoi_cands.R")

library(SobolSequence)
library(lhs)

set.seed(1234)

f2 <- 10
Nx <- 2^f2
N <- 15

#cols <- primary.colors(N)

for (P in c(2,10)) {
    X <- matrix(runif(N*P), ncol=P)

    if (P==10) {
        cols <- rgb.tables(N)
    } else {
        #cols <- rgb.tables(4)
        cols <- primary.colors(6)[2:5]
    }

    for (target in c("Sobol","Trust_Region","Tricands","Vorcands")) {
        if (target=='Sobol') {
            Xcand <- sobolSequence.points(dimR=P, dimF2=f2, count=Nx)
        } else if (target=='Trust_Region') {
            #Xcand <- randomLHS(Nx,P)
            Xcand <- sobolSequence.points(dimR=P, dimF2=f2, count=Nx)
            #optind <- which.min(apply(X,1,function(x) sum((x-0.5)^2)))
            #optind <- which.min(apply(X,1,function(x) max(abs(x-0.5))))
            #optind <- which.min(apply(X,1,function(x) sum((x[1:2]-c(0.75,0.5))^2)))
            optind <- which.min(apply(X,1,function(x) sum((x[1:2]-c(0.70,0.70))^2)))
            #optind <- which.min(apply(X,1,function(x) sum((x[1:2]-c(0.75,0.6))^2)))
            x <- X[optind,]
            rad <- c(0.1, 0.2) # NOTE: is recycled 5 times for P=10 case.
            ub <- pmin(x+rad,1)
            lb <- pmax(x-rad,0)
            Xcand <- t(apply(Xcand, 1, function(x) x*(ub-lb)+lb))
        } else if (target=='Tricands') {
            Xcand <- tricands(X, max = Nx)
        } else if (target=='Vorcands') {
            Xcand <- vorwalkcands(X, ncand = Nx, st = 'unif', norm = 'linf')$Xs
        } else {
            stop("Unknown Candidates.")
        }
        tit <- paste(gsub('_',' ', target), " P=",P,sep='')

        #if (P==10) {
        if (TRUE) {
            same_orthant_as <- rep(NA, Nx)
            for (nc in 1:nrow(Xcand)) {
                orth_c <- as.numeric(Xcand[nc,]>0.5)
                for (n in 1:N) {
                    orth_x <- as.numeric(X[n,]>0.5)
                    if (all(orth_x==orth_c)) {
                        if (is.na(same_orthant_as[nc])) {
                            same_orthant_as[nc] <- n
                        } else {
                            if (P==10)
                                stop("Double orthant occupation.")
                        }
                    }
                } 
            }

            if (P==2) {
                pointcols <- rep(NA, length(N))
                pointcols[1] <- cols[1]
                usedcols <- 1
                for (n in 2:N) {
                    xn <- X[n,]
                    poa <- rep(NA, n-1)
                    for (ni in 1:(n-1)) {
                        poa[ni] <- all((xn>0.5)==(X[ni,]>0.5))
                    }
                    if (any(poa)) {
                        pointcols[n] <- pointcols[which(poa)[1]]
                    } else {
                        usedcols <- usedcols + 1
                        pointcols[n] <- cols[usedcols]
                    }
                }
            } else {
                pointcols <- cols
            }

            vals <- unique(same_orthant_as)
            neighs <- vals[!is.na(vals)]
            cand_cols <- rep(NA, nrow(Xcand))
            cand_lwd <- rep(NA, nrow(Xcand))
            for (n in 1:nrow(Xcand)) {
                if (!is.na(same_orthant_as[n])) {
                    cand_cols[n] = pointcols[same_orthant_as[n]]
                    cand_lwd[n] = 2
                } else {
                    cand_cols[n] = 'gray'
                    cand_lwd[n] = 1
                }
            }
        } else {
            cand_cols <- 'gray'
            cand_lwd <- 1
        }

        pdf(paste("figures/",target,"_example_",P,".pdf",sep=''), width = 4, height = 4)
        par(mar=c(1,1,1,0)+0.1)
        par(mgp=c(0,0,0)+0.1)
        plot(NA,NA,xlim=c(0,1),ylim=c(0,1), main = tit, xlab = '', ylab = '', tck = -0.005)
        points(Xcand[,1],Xcand[,2], col = cand_cols, lwd = cand_lwd)
        #
        points(X[,1],X[,2], col = 'white', pch = 4, lwd = 12)
        points(X[,1],X[,2], col = 'black', pch = 4, lwd = 8)
        points(X[,1],X[,2], col = pointcols, pch = 4, lwd = 4)
        #
        dev.off()

    }
}
