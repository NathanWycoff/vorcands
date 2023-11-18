source("R/vornoi_cands.R")

set.seed(123)

#lowres <- TRUE
lowres <- FALSE

if (lowres) {
    Nc <- 50
} else {
    Nc <- 500
}

sts <- c('rect','unif','lhs')
Ps <- c(2,10,100)
Ns <- c(10, 100, 500, 1000)
norms <- c("l1","l2","linf")
if (lowres) {
    reps <- 2
} else {
    reps <- 30
}

resdf <- as.data.frame(matrix(0, nrow = length(sts)*length(norms)*length(Ps)*length(Ns)*reps, ncol = 5))
colnames(resdf) <- c("POB","P","N","norm",'st')

#st <- 'rect' 
ind <- 1
for (st in sts) {
    eps <- 1e-8

    tt <- Sys.time()
    for (pi in 1:length(Ps)) {
        P <- Ps[pi]
        for (ni in 1:length(Ns)) {
            N <- Ns[ni]
            for (rep in 1:reps) {
                X <- matrix(runif(N*P, min = 0.1, max = 0.9), ncol=P)

                for (norm in norms) {
                    Xcand <- vorwalkcands(X, ncand = Nc, st = st, norm = norm, half2bound = F, in_iters = 30)$Xs
                    res <- mean(dist2bound(Xcand) <= eps)
                    resdf[ind,1] <- res
                    resdf[ind,2] <- P
                    resdf[ind,3] <- N
                    resdf[ind,4] <- norm
                    resdf[ind,5] <- st
                    ind <- ind + 1
                }
            }
        }
    }
    tt1 <- Sys.time()
    print(tt1-tt)

}
save(resdf, file="Rdata/bound.RData")