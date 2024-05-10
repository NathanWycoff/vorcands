library(Rcpp)

getXs <- function(ts, Xb, U) {
    Xb + ts * U
}

#N1 <- 8
#N2 <- 3
#P <- 3
#X1 <- matrix(runif(N1*P), ncol = P)
#X2 <- matrix(runif(N2*P), ncol = P)
#
#dist_unnorm(X1, X2, norm = 1)
#as.matrix(dist(rbind(X1,X2),method='manhattan'))[1:nrow(X1),(nrow(X1)+1):(nrow(X1)+nrow(X2))]
#
#dist_unnorm(X1, X2, norm = 2)
#as.matrix(dist(rbind(X1,X2),method='euclidean'))[1:nrow(X1),(nrow(X1)+1):(nrow(X1)+nrow(X2))]^2
#
#dist_unnorm(X1, X2, norm = 3)
#as.matrix(dist(rbind(X1,X2),method='maximum'))[1:nrow(X1),(nrow(X1)+1):(nrow(X1)+nrow(X2))]

#cppFunction( "
#NumericMatrix dist_unnorm(NumericMatrix X1, NumericMatrix X2, int norm) {
#    //std::cout << \"Ye, boi!\" << std::endl;
#    int N1 = X1.nrow();
#    int N2 = X2.nrow();
#    int P = X1.ncol();
#    if (X2.ncol() != P) {
#        Rcpp::stop(\"Column sizes in X1,X2 must match.\");
#    }
#    NumericMatrix D(N1,N2);
#    double d;
#
#    if ((norm!=1) && (norm!=2) && (norm!=3)) {
#        Rcpp::stop(\"norm must be 1, 2 or 3.\");
#    }
#
#    for (int n1=0; n1<N1; n1++) {
#        for (int n2=0; n2<N2; n2++) {
#            for (int p=0; p<P; p++) {
#                d = fabs(X1(n1,p) - X2(n2,p));
#                if (norm==1) {
#                    D(n1,n2) += d;
#                } else if (norm==2) {
#                    D(n1,n2) += d*d;
#                } else if (norm==3) {
#                    if (d > D(n1,n2)) {
#                        D(n1,n2) = d;
#                    }
#                }
#            }
#        }
#    }
#    return(D);
#} ")
#
#get_nn <- function(X1, X2, norm = 'l2') {
#    if (norm=='l1') {
#        cnorm <- 1
#    } else if (norm=='l2') {
#        cnorm <- 2
#    } else if (norm=='linf') {
#        cnorm <- 3
#    } else {
#        stop("Bad norm.")
#    }
#    D <- dist_unnorm(X1,X2,cnorm)
#    neigh <- apply(D, 2, which.min)
#
#    ret <- list()
#    ret$nn.idx <- neigh
#    ret$dist <- sapply(1:nrow(X2), function(i) D[neigh[i],i])
#    return(ret)
#}
#
#getsign <- function(X, Xs, ns, norm = 'l2') {
##    if (norm=='l1') {
##        cnorm <- 1
##    } else if (norm=='l2') {
##        cnorm <- 2
##    } else if (norm=='linf') {
##        cnorm <- 3
##    } else {
##        stop("Bad norm.")
##    }
##    D <- dist_unnorm(X,Xs,cnorm)
##    neigh <- apply(D, 2, which.min)
##
##    ret <- list()
#    nn <- get_nn(X, Xs, norm = norm)
#    ret <- list()
#    ret$sign <- nn$nn.idx==ns
#    if (norm=='l2') {
#        nn$dist <- sqrt(nn$dist)
#    }
#    ret$dist <- nn$dist
#
#    return(ret)
#}

get_nn <- function(X1, X2, K, norm = 'l2') {
    if (norm=='l2') {
        neigh <- RANN::nn2(X1, X2, K)
        #stop("Not implemented for portability reasons.")
    } else if (norm=='l1') {
        neigh <- RANN.L1::nn2(X1, X2, K)
        #stop("Not implemented for portability reasons.")
    } else if (norm=='linf') {
        neigh <- RANN.Linf::nn2(X1, X2, K)
        ##D <- matrix(NA, nrow=nrow(X1),ncol=nrow(X2))
        ##for (n in 1:nrow(X1)) {
        ##    for (nc in 1:nrow(X2)) {
        ##        D[n,nc] <- max(abs(X1[n,]-X2[nc,]))
        ##    }
        ##}
        #D <- as.matrix(dist(rbind(X1,X2),method='maximum'))[1:nrow(X1),(nrow(X1)+1):(nrow(X1)+nrow(X2))]
        #if (K>1) stop("K>1 not implemented for linf norm.")
        ##apply(D, 2, rank)
        #neigh <- list()
        #neigh$nn.idx <- matrix(apply(D, 2, which.min), ncol = 1)
        #neigh$nn.dists <- matrix(apply(D, 2, min), ncol = 1)
    } else {
        stop(paste("Unknown norm:",norm))
    }
    return(neigh)
}

getsign <- function(X, Xs, ns, norm = 'l2') {
    neigh <- get_nn(X, Xs, K=1, norm=norm)
    ret <- list()
    ret$sign <- neigh$nn.idx[,1]==ns
    ret$dist <- neigh$nn.dists[,1]

    # Account for the case with many ties
    initdists <- apply(X[ns,]-Xs, 1, function(x) lpnorm(x, norm = norm))
    ret$sign <- ret$sign | (ret$dist==initdists)
    return(ret)
}

lpnorm <- function(x, norm = 'l2') {
    if (norm=='l2') {
        nf <- sqrt(sum(x^2))
    } else if (norm=='l1') {
        nf <- sum(abs(x))
    } else if (norm=='linf') {
        nf <- max(abs(x))
    } else {
        stop("Bad norm.")
    }
    return(nf)
}

dist2bound <- function(X) {
    return(apply(X, 1, function(x) min(pmin(x,1-x))))
}

#N <- 4
#P <- 3
#X <- matrix(runif(N*P), ncol = P)
#U <- matrix(rnorm(N*P), ncol = P)

cppFunction("
NumericVector dist2boundalong_C(NumericMatrix X, NumericMatrix U, int norm) {
    int N = X.nrow();
    int P = X.ncol();

    double eps = 1e-10;
    double d,d0,d1,ma;

    static double BIG_DOUBLE = 1000000.;
    NumericVector travel(N, BIG_DOUBLE);
    NumericVector dist(N);

    for (int n=0; n<N; n++) {
        // Compute distance along line segment to travel
        for (int p=0; p<P; p++) {
            if (fabs(U(n,p)) > eps) {
                d0 = (1-X(n,p)) / U(n,p);
                d1 = -X(n,p) / U(n,p);
                ma = std::max(d0,d1);
                if (ma < travel(n)) {
                    travel(n) = ma;
                }
            }
        }
        // Compute lp norm of that vector.
        for (int p=0; p<P; p++) {
            d = fabs(travel(n) * U(n,p));
            if (norm==1) {
                dist(n) += d;
            } else if (norm==2) {
                dist(n) += d*d;
            } else if (norm==3) {
                if (d > dist(n)) {
                    dist(n) = d;
                }
            }
        }
        if (norm==2) {
            dist(n) = sqrt(dist(n));
        }
    }

    // If anybody is outside of the unit cube, we define them to have zero distance to boundary.
    NumericVector is_oob(N);
    for (int n=0; n<N; n++){
        for (int p=0; p<P; p++) {
            if (X(n,p) < 0. || X(n,p) > 1.) {
                is_oob(n) = 1.;
            }
        }
        if (is_oob(n) > 0) {
            dist(n) = 0;
        }
    }

    return(dist);
}")

dist2boundalong <- function(X, U, norm = 'l2') {
    if (norm=='l1') {
        cnorm <- 1
    } else if (norm=='l2') {
        cnorm <- 2
    } else if (norm=='linf') {
        cnorm <- 3
    } else {
        stop("Bad norm.")
    }

    return(dist2boundalong_C(X, U, cnorm))
}

## How far along U[i,] must X[i,] travel to reach the boundary?
##X <- Xs
#dist2boundalong <- function(X, U, norm = 'l2') {
#
#    # TODO: double check this na.rm
#    travel <- sapply(1:nrow(X), function(n) min(pmax((1-X[n,])/U[n,], -X[n,]/U[n,]), na.rm = T))
#    dists <- sapply(1:nrow(X), function(n) lpnorm(travel[n]*U[n,], norm=norm))
#
#    oob <- (apply(X, 1, min) < 0) | (apply(X, 1, max) > 1)
#    dists[oob] <- 0
#
#    return(dists)
#}

#vorwalk <- function(X, ns, U, norm = 'l2', in_iters = 20, half2bound = TRUE) {
vorwalk <- function(X, ns, U, norm = 'l2', in_iters = 10, half2bound = TRUE, plot = FALSE) {
    ncand <- length(ns)

    if (nrow(U)!=ncand) stop("Nonconforming ns and U arguments.")

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
        if (half2bound) {
            bound_dists <- dist2boundalong(Xs, U, norm = norm)
            in_upper <- in_upper & (bound_dists >= fmid$dist)
            nd2 <- sum(!in_upper) - nd1
            if (it == in_iters) {
                bdp <- 100*nd2/(nd1+nd2)
                if (nd1+nd2==0) {
                    print("Warn: no decreases! Should be about 50-50 at convergence.")
                }
            }
        } else {
            outofbounds <- apply(abs(Xs-0.5)>0.5, 1, any)
            in_upper <- in_upper & (!outofbounds)
        }

        lb[in_upper] <- mid[in_upper]
        ub[!in_upper] <- mid[!in_upper]
        T_trace[it,] <- mid
    }

    rownames(Xs) <- ns

    if (plot) {
        pdf("t_trace.pdf")
        plot(NA,NA,xlim=c(0,nrow(T_trace)+1), ylim = c(0,1))
        for (nc in 1:ncand) {
            points(T_trace[,nc], type = 'l')
        }
        dev.off()
    }

    #changed_mind <- expand_point & !expand_bound
    Xs <- pmax(pmin(Xs,1),0)
    ret <- list(Xs=Xs,ns=ns)
    return(ret)
}

#vorwalkcands <- function(X, ncand = 1000, y = NULL, st = 'rect', in_iters = 20, norm = 'l2', half2bound = TRUE) {
vorwalkcands <- function(X, ncand = 1000, y = NULL, st = 'rect', norm = 'linf', lb = rep(0, ncol(X)), ub = rep(1, ncol(X)), ...) {
    X <- t((t(X) - lb) / (ub-lb))

    if ('matrix' %in% class(st)) {
        matin <- TRUE
        ncand <- nrow(st)
        if (!ncol(X)==ncol(st)) stop("Points to project must be in same space.")
    } else {
        matin <- FALSE
    }

    N <- nrow(X)
    P <- ncol(X)
    
    # Let LHS determine our starting points and moving directions.
    if (matin || st=='lhs') {
        if (matin) {
            Xin <- st
        } else {
            Xin <- lhs::randomLHS(ncand, k = P)
        }
        ns <- get_nn(X, Xin, 1, norm = norm)$nn.idx

        U <- matrix(NA, nrow = ncand, ncol = P)
        for (nc in 1:ncand) {
            U[nc,] <- Xin[nc,] - X[ns[nc],]
        }
        U <- U / sqrt(rowSums(U^2))
    } else {
        ## Select points to start from.
        if (st=='rect') {
            ns <- rep(1:N, rep(2*P,N))
            if (ncand >= 2*N*P) {
                ncand <- 2*N*P
            } else {
                nopt <- 2*P
                if (!is.null(y)) {
                    best <- which.min(y)
                    optinds <- which(ns==best)
                    ni <- sample(setdiff(1:(2*P*N), optinds),ncand-length(optinds),replace=FALSE)
                    ni <- c(optinds, ni)
                } else {
                    ni <- sample(1:(2*P*N),ncand,replace=FALSE)
                }
                ns <- ns[ni]
            }
        } else {
            if (!is.null(y)) {
                best <- which.min(y)
                nopt <- floor(ncand*0.25)
                ns <- rep(best,nopt)
                if (N > 1) {
                    ns <- c(ns, sample(setdiff(1:N,best),ncand-nopt,replace=T))
                } else {
                    ns <- c(ns, rep(best, ncand-nopt))
                }
            } else {
                ns <- sample(N,ncand,replace=T)
            }
        }

        ## Select direction to move in.
        if (st=='rect') {
            U <- matrix(0,2*N*P,P)
            for (p in 1:P) {
                up <- (1:nrow(U)-1) %% (2*P) == 2*(p-1)
                down <- (1:nrow(U)-1) %% (2*P) == 2*(p-1)+1
                U[up,p] <- 1
                U[down,p] <- -1
            }
            if (ncand < 2*N*P) {
                U <- U[ni,]
            }
        } else if (st == 'unif') {
            U <- matrix(rnorm(ncand*P), nrow = ncand)
            U <- U / sqrt(rowSums(U^2))
            if (st=='active') {
                browser()
                Rdimtools::do.sir
            }
        } else {
            stop("Bad st!")
        }
    }

    ret <- vorwalk(X, ns, U, norm, ...)
    ret$Xs <- t(t(ret$Xs) * (ub-lb) + lb)
    return(ret)
}
