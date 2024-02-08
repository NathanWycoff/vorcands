#source("R/vorlhs.R")
#source("R/vornoi_fun.R")
source("R/vornoi_cands.R")
library(stringr)
library(hetGP)
library(akima)

##
## begin laGP helper functions
##

library(laGP)

## EI:
## 
## return Expected improvement on candidate(s) x

EI <- function(gpi, x, fmin, pred=predGPsep, noise=FALSE)
 {
  if(is.null(nrow(x))) x <- matrix(x, nrow=1)
  cnt <<- cnt + nrow(x)
  p <- pred(gpi, x, lite=TRUE, nonug=noise)
  d <- fmin - p$mean
  sigma <- sqrt(p$s2)
  dn <- d/sigma
  ei <- d*pnorm(dn) + sigma*dnorm(dn)
  return(ei)
}


get_cands <- function(X, m, ncands, cands = 'lhs', cand_params = NULL, y = NULL) {
  if(substr(cands,1,3)=='tri') {
    Xcand <- tricands(X, max=ncands, best=m, cand_params=cand_params)
  } else if (substr(cands,1,3)=='vor') {

    if (cand_params$style=='rect' || (cand_params$style=='alt' && nrow(X) %% 2 == 0)) {
      Xcand <- vorwalkcands(X, ncands, y=y, st = cand_params$vor_st, norm = cand_params$vor_norm, half2bound = T)$Xs
      rownames(Xcand) <- paste('rect.',rownames(Xcand),sep='')
    } else if (cand_params$style=='lhs' || (cand_params$style=='alt' && nrow(X) %% 2 == 1)) {
      #Xcand <- vorwalkcands(X, ncands, y=y, st = 'lhs', norm = 'l2', half2bound = T)$Xs
      if (nrow(X)==100) {
        print("Using norm for both LHS and Rect.")
        print(cand_params$vor_st)
      }
      Xcand <- vorwalkcands(X, ncands, y=y, st = 'lhs', norm = cand_params$vor_norm, half2bound = T)$Xs
      rownames(Xcand) <- paste('lhs.',rownames(Xcand),sep='')
    } else {
      stop("Unrecognized candidate style.")
    }

  } else if (cands=='hitri') {

    Xc1 <- hitricands(X, ncands/2, best = m, astrat = 'unif', half2bound = T)
    Xc2 <- hitricands(X, ncands/2, best = m, astrat = 'nn', half2bound = T)
    Xcand <- rbind(Xc1, Xc2)

  } else if (cands=='corner') {
    Xcand <- matrix(sample(c(0,1),ncands*ncol(X), replace=T), ncol = ncol(X))
    Xcand <- Xcand[!duplicated(Xcand),]
  } else if (cands=='lhs') {
    Xcand <- randomLHS(ncands, ncol(X))
  } else {
    stop("Unrecognized candidate name.")
  }
  return(Xcand)
}

## EI.cands:
## 
## new EI fucntion that evaluates EI on triangulated 
## gap-filling candidates

EI.cands <- function(X, y, gpi, ncands=100*ncol(X),
  cands='tri', noise=FALSE, pred=predGPsep, cand_params = cand_params )
 {
  if(noise) ym <- pred(gpi, X, lite=TRUE)$mean
  else ym <- y
  m <- which.min(ym)
  fmin <- ym[m]
  Xcand <- get_cands(X=X, m=m, ncands=ncands, cands=cands, cand_params = cand_params, y = y)
  solns <- data.frame(Xcand, Xcand, EI(gpi, Xcand, fmin, pred=pred, noise=noise))
  names(solns) <- c(paste0("s", 1:ncol(X)), paste0("x", 1:ncol(X)), "val") 
  return(solns)
}

## EY.cands:
## 
## new EY fucntion that evaluates EY on triangulated 
## gap-filling candidates

EY.cands <- function(X, y, gpi, ncands=ncol(X)*100, 
  cands='tri', noise=FALSE, pred=predGPsep, cand_params = cand_params)
 {
  if(noise) ym <- pred(gpi, X, lite=TRUE)$mean
  else ym <- y
  m <- which.min(ym)
  ## fmin <- ym[m]
  Xcand <- get_cands(X=X, m=m, ncands=ncands, cands=cands, cand_params = cand_params, y = y)
  cnt <<- cnt + nrow(Xcand)
  solns <- data.frame(Xcand, Xcand, - pred(gpi, Xcand, lite=TRUE, nonug=noise)$mean)
  names(solns) <- c(paste0("s", 1:ncol(X)), paste0("x", 1:ncol(X)), "val") 
  return(solns)
}

## PI:
##
## return probability of improvement on candidate(s) x

PI <- function(gpi, x, fmin, pred=predGPsep, noise=FALSE)
 {
  if(is.null(nrow(x))) x <- matrix(x, nrow=1)
  cnt <<- cnt + nrow(x)
  p <- pred(gpi, x, lite=TRUE, nonug=noise)
  d <- fmin - p$mean
  sigma <- sqrt(p$s2) 
  dn <- d/sigma
  pi <- pnorm(dn) 
  return(pi)
 }

## PI.cands:
## 
## new PI fucntion that evaluates PI on triangulated 
## gap-filling candidates

PI.cands <- function(X, y, gpi, ncands=ncol(X)*100, 
  cands='tri', noise=FALSE, pred=predGPsep, cand_params = cand_params)
 {
  if(noise) ym <- pred(gpi, X, lite=TRUE)$mean
  else ym <- y
  m <- which.min(ym)
  fmin <- ym[m]
  Xcand <- get_cands(X=X, m=m, ncands=ncands, cands=cands, cand_params = cand_params, y = y)
  solns <- data.frame(Xcand, Xcand, PI(gpi, Xcand, fmin, pred=pred, noise=noise))
  names(solns) <- c(paste0("s", 1:ncol(X)), paste0("x", 1:ncol(X)), "val") 
  return(solns)
}

## optim.crit
##
## generic BO criteria optimizer, coppied with minor modification
## from Chapter 7 of Surrogates homework solutions

eps <- sqrt(.Machine$double.eps)
library(lhs)
optim.crit <- function(X, y, gpi, obj, check, noise=FALSE, multi.start=2*ncol(X)+1, pred=predGPsep, tol=eps, ...)
 {
  if(noise) ym <- pred(gpi, X, lite=TRUE)$mean
  else ym <- y
  m <- which.min(ym)
  fmin <- ym[m]
  ## might not want to guarantee that we search from the best value
  ## For optim start
  #start <- matrix(X[m,], nrow=1)
  #if(multi.start > 1) 
  #  start <- rbind(start, randomLHS(multi.start-1, ncol(X)))
  ## For optim start
  ## For dumb start
  start <- randomLHS(multi.start, ncol(X))
  ## For dumb start
  xnew <- matrix(NA, nrow=nrow(start), ncol=ncol(X)+1)
  for(i in 1:nrow(start)) {
    if(check(gpi, start[i,], fmin, noise, tol)) { out <- list(value=Inf) }
    else { out <- optim(start[i,], obj, method="L-BFGS-B", 
      lower=0, upper=1, gpi=gpi, fmin=fmin, noise=noise, ...) }
    xnew[i,] <- c(out$par, -out$value)
  }
  solns <- data.frame(cbind(start, xnew))
  names(solns) <- c(paste0("s", 1:ncol(X)), paste0("x", 1:ncol(X)), "val") 
  solns <- solns[solns$val > tol,]
  
  ## random search if empty  
  if(nrow(solns) == 0) {
    start <- matrix(runif(ncol(X)), nrow=1)
    solns[1,] <- c(start, start, -obj(start, gpi=gpi, fmin=fmin, noise=noise, ...))
  }
  return(solns)
}


## thompson.samp
##
## Thomposon Sapling for BO, coppied with minor modification
## from Chapter 7 of Surrogates homework solutions

#library(Rfast) ## Rfast has rmvnorm too
thompson.samp <- function(gpi, X, y, noise=FALSE, ncands=100*ncol(X), close=0, pred=predGPsep)
 {
  ## find fmin
  if(noise) ym <- pred(gpi, X, lite=TRUE)$mean
  else ym <- y
  m <- which.min(ym)

  ## set candidates for random evaluation
  if(close >= 0) Xcand <- randomLHS(ncands-close, ncol(X))
  else Xcand <- tricands(X, max=ncands, best=m)

  ## augment with nearby candidates in close-best box
  if(close > 0) {

    Xc <- X[ym <= nth(ym, 5),]
    r <- apply(Xc, 2, range)
    Xcc <- matrix(runif(close*ncol(X)), ncol=ncol(X))
    for(j in 1:ncol(Xcc)) 
      Xcc[,j] <- Xcc[,j]*(r[2,j] - r[1,j]) + r[1,j]
    Xcand <- rbind(Xcand, Xcc)
  } 

  ## keep track of number of candidates
  cnt <<- cnt + nrow(Xcand)

  ## predict on the candidate set
  p <- pred(gpi, Xcand, nonug=noise)
  Ycand <- rmvnorm(1, p$mean, p$Sigma)

  ## report the best as a dummy data.frame
  m <- which.min(Ycand)
  solns <- as.data.frame(matrix(c(rep(0, ncol(X)), Xcand[m,], Ycand[m]), nrow=1))
  names(solns) <- c(paste0("s", 1:ncol(X)), paste0("x", 1:ncol(X)), "val") 
  return(solns)
}


## optim.surr:
##
## borrowed from Chapter 7 of Surrogates homework solution(s), modified to 
## use custom initialization, and to allow tricands candidates

library(laGP)
# ssm = "Start Skipping MLE": sample size at which we don't update like every iter.
# ume = "Update MLE Every": Once we do start skipping, how often do we update?
optim.surr <- function(f, ninit, m, end, X=NULL, sur = 'gp',
    criteria=c("EI", "PI", "EY", "TS", "IECI"), cands="opt",
    ncands=100*m, close=0, noise=FALSE, dmin=0.001, dmax=0.5, cand_params=NULL, viz_design = FALSE,
    ssm = 250, ume = 25)
 {
  ## run on initial design
  if(is.null(X)) X <- randomLHS(ninit, m)
  if(nrow(X) != ninit || ncol(X) != m) stop("nrow(X) must match ninit and m")
  y_raw <- f(X)
  mu_y <- mean(y_raw)
  sig_y <- sd(y_raw)
  y <- (y_raw-mu_y) / sig_y

  if (sur=='gp') {
    ## select an optimization criterion
    #cands <- match.arg(cands)
    criteria <- match.arg(criteria)
    ## this next bit is irrelevant if cands != "opt"
    if(criteria == "EI") {
      obj <- function(x, fmin, gpi, noise) - EI(gpi, x, fmin, predGPsep, noise)
      check <- function(gpi, x, fmin, noise, tol) { EI(gpi, x, fmin, predGPsep, noise) <= tol }
      tol <- eps
    } else if(criteria == "PI") {
      obj <- function(x, fmin, gpi, noise) - PI(gpi, x, fmin, predGPsep, noise)
      check <- function(gpi, x, fmin, noise, tol) { PI(gpi, x, fmin, predGPsep, noise) <= tol }
      tol <- eps
    } else if(criteria == "EY") {
      obj <- function(x, fmin, gpi, noise) { cnt <<- cnt + 1; predGPsep(gpi, matrix(x, nrow=1), lite=TRUE, nonug=noise)$mean }
      check <- function(gpi, x, fmin, noise, tol) FALSE
      tol <- -Inf
    } else if(criteria == "IECI") {
      obj <- function(x, fmin, gpi, noise, Xref) ieciGPsep(gpi, matrix(x, nrow=1), fmin, nonug=noise, Xref=Xref)
      check <- function(gpi, x, fmin, noise, tol) FALSE
      tol <- -Inf
    }

    ## fit initial GP model, differentially depending on whether or not noise is modeled
    ### DUP BLOCK A
    da <- darg(list(mle=TRUE, start = (dmin+dmax)/2, min=dmin, max=dmax), randomLHS(1000, ncol(X)))
    if(noise) {
      ga <- garg(list(mle=TRUE), y)
      #gpi <- newGPsep(X, y, d=0.01, g=ga$start, dK=TRUE)
      gpi <- newGPsep(X, y, d=da$start, g=ga$start, dK=TRUE)
      mleGPsep(gpi, param="both", tmin=c(da$min, ga$min), tmax=c(da$max, ga$max), ab=c(da$ab, ga$ab))
    } else {
      gpi <- newGPsep(X, y, d=da$start, g=1e-6, dK=TRUE)
      mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)
    }
    ### DUP BLOCK A

    ## optimization loop
    cnt <<- 0
    #best <- c()
    crits <- rep(0, end-ninit)
    for(i in (ninit+1):end) {

      ## slight variations on calls for the different criteria
      if(criteria == "EI" && cands != "opt") solns <- EI.cands(X, y, gpi, noise, ncands=ncands, cands=cands, cand_params=cand_params) 
      else if(criteria == "EY" && cands != "opt") solns <- EY.cands(X, y, gpi, noise, ncands=ncands, cands=cands, cand_params=cand_params) 
      else if(criteria == "PI" && cands != "opt") solns <- PI.cands(X, y, gpi, noise, ncands=ncands, cands=cands, cand_params=cand_params) 
      else if(criteria == "IECI") solns <- optim.crit(X, y, gpi, obj, check, noise, Xref=randomLHS(100, m), tol=tol)
      else if(criteria != "TS") solns <- optim.crit(X, y, gpi, obj, check, noise, tol=tol)
      else solns <- thompson.samp(gpi, X, y, noise, ncands=ncands, close=close)

      Xc <- solns[,1:ncol(X)]
      #print(dim(Xc))
      wm <- which.max(solns$val)
      crits[i-ninit] <- max(solns$val)
      #print(paste("Max ei is:", max(solns$val)))
      camefrom <- as.numeric(sapply(strsplit(substr(str_extract(rownames(Xc), "X[0-9]+\\.*"), 2, 10000L),'\\.'), function(x) x[[1]]))
      cf <- camefrom[wm]

      xnew <- as.matrix(solns[wm,(ncol(X)+1):(2*ncol(X))])
      ynew_raw <- f(xnew)
      ynew <- (ynew_raw - mu_y) / sig_y
      updateGPsep(gpi, xnew, ynew)

      ## keep track of grand design and runs    
      X <- rbind(X, xnew)
      y_raw <- c(y_raw, ynew_raw)
      y <- c(y, ynew)

      ## Renormalize y and refit GP every now and then.
      if (nrow(X) <= ssm || nrow(X) %% ume == 0) {
        #if(noise) mleGPsep(gpi, param="both", tmin=c(da$min, ga$min), tmax=c(da$max, ga$max), ab=c(da$ab, ga$ab)) 
        #else mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)

        mu_y <- mean(y_raw)
        sig_y <- sd(y_raw)
        y <- (y_raw-mu_y) / sig_y

        ### DUP BLOCK A
        da <- darg(list(mle=TRUE, start = (dmin+dmax)/2, min=dmin, max=dmax), randomLHS(1000, ncol(X)))
        if(noise) {
          ga <- garg(list(mle=TRUE), y)
          #gpi <- newGPsep(X, y, d=0.01, g=ga$start, dK=TRUE)
          gpi <- newGPsep(X, y, d=da$start, g=ga$start, dK=TRUE)
          mleGPsep(gpi, param="both", tmin=c(da$min, ga$min), tmax=c(da$max, ga$max), ab=c(da$ab, ga$ab))
        } else {
          gpi <- newGPsep(X, y, d=da$start, g=1e-6, dK=TRUE)
          mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)
        }
        ### DUP BLOCK A
      }

    }

    ## clean up and return
    deleteGPsep(gpi)

    #return(list(X=X, y=y, best=best, cnt=cnt))
    return(list(X=X, y=y_raw, cnt = cnt, crits=crits))
  } else if (sur=='tgp') {
   for(i in (ninit+1):end) {
     #out <- optim.step.tgp.simple(f, X = X, Z = y, cands=cands, verb = 0, NN=ncands, nug.p=0, gd=c(1e-5,2), cand_params=cand_params)
     out <- optim.step.tgp.simple(f, X = X, Z = y, cands=cands, verb = 0, NN=ncands, nug.p=0, cand_params=cand_params)
     X <- rbind(X, out)
     y <- c(y, f(out))
   }

   return(list(X=X, y=y, crits=crits))
  } else if (sur=='deepgp') {
   for(i in (ninit+1):end) {
     out <- optim.step.deepgp.simple(f, X = X, Z = y, cands=cands, verb=FALSE, NN=ncands, true_g=1e-7, cand_params=cand_params)
     X <- rbind(X, out)
     y <- c(y, f(out))
   }

   return(list(X=X, y=y, crits=crits))
  } else {
    stop("Unknown Surrogate!")
  }
 }


## bov:
##
## Calculatates best observed value for BO progress
## copied from Surrogates

bov <- function(y, end=length(y))
 {
  prog <- rep(min(y), end)
  prog[1:min(end, length(y))] <- y[1:min(end, length(y))]
  for(i in 2:end) { 
    if(is.na(prog[i]) || prog[i] > prog[i-1]) 
      prog[i] <- prog[i-1]
  }
  return(prog)
 }
