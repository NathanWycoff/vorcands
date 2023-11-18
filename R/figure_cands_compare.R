## This makes the figure with various candidate schemes.

library(colorRamps)
source("tricands/R/tricands.R")
source("R/vornoi_cands.R")

library(SobolSequence)
library(lhs)

set.seed(123)

f2 <- 10
Nx <- 2^f2
N <- 15

#cols <- primary.colors(N)
cols <- rgb.tables(N)

for (P in c(2,10)) {
    X <- matrix(runif(N*P), ncol=P)

    for (target in c("Sobol","LHS","Tricands","Vorcands")) {
        if (target=='Sobol') {
            Xcand <- sobolSequence.points(dimR=P, dimF2=f2, count=Nx)
        } else if (target=='LHS') {
            Xcand <- randomLHS(Nx,P)
        } else if (target=='Tricands') {
            Xcand <- tricands(X, max = Nx)
        } else if (target=='Vorcands') {
            Xcand <- vorwalkcands(X, ncand = Nx, st = 'unif', norm = 'linf')$Xs
        } else {
            stop("Unknown Candidates.")
        }
        tit <- paste(target, " P=",P,sep='')

        if (P==10) {
            same_orthant_as <- rep(NA, Nx)
            for (nc in 1:nrow(Xcand)) {
                orth_c <- as.numeric(Xcand[nc,]>0.5)
                for (n in 1:N) {
                    orth_x <- as.numeric(X[n,]>0.5)
                    if (all(orth_x==orth_c)) {
                        if (is.na(same_orthant_as[nc])) {
                            same_orthant_as[nc] <- n
                        } else {
                            stop("Double orthant occupation.")
                        }
                    }
                } 
            }

            vals <- unique(same_orthant_as)
            neighs <- vals[!is.na(vals)]
            cand_cols <- rep(NA, nrow(Xcand))
            cand_lwd <- rep(NA, nrow(Xcand))
            for (n in 1:nrow(Xcand)) {
                if (!is.na(same_orthant_as[n])) {
                    cand_cols[n] = cols[same_orthant_as[n]]
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
        points(X[,1],X[,2], col = 'black', pch = 4, lwd = 6)
        points(X[,1],X[,2], col = cols, pch = 4, lwd = 4)
        #
        dev.off()

    }
}