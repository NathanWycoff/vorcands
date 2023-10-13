## main tricands (tri) v optim toggle on line 54

library(R.matlab)
library(hetGP)
source("tricands.R")

## ATO:
##
## Stand-alone ATO R interface function
## uses matlab global variable set via code above;
## don't forget to close matlab when you're done

ATO <- function(Xnew, time=20, seed) 
{
  ## for normalizing output z-values
  zm <- 43
  zs <- 35
  
  if(is.null(nrow(Xnew))) Xnew <- matrix(Xnew, nrow=1)
  if(max(Xnew) <= 1) Xnew <- Xnew*19 + 1
  out <- rep(NA, nrow(Xnew))
  for(i in 1:nrow(Xnew)){
    input <- as.numeric(Xnew[i,])
    setVariable(matlab, input=input)
    setVariable(matlab, seed=seed[i])
    setVariable(matlab, time=time)
    evaluate(matlab, "output=ATO(input', time, seed);")
    output <- getVariable(matlab, "output")
    out[i] <- output
  }
  out <- (as.numeric(out)-zm)/zs
  return(out)
}


## bov: best observed value for keeping track of BO
## progress

bov <- function(y, end=length(y))
 {
  prog <- rep(min(y), end)
  prog[1:min(end, length(y))] <- y[1:min(end, length(y))]
  for(i in 2:end) 
    if(is.na(prog[i]) || prog[i] > prog[i-1]) prog[i] <- prog[i-1]
  return(prog)
 }


## problem specification
ninit <- 80
ntot <- 300
dim <- 8
tri <- 2 ## TRUE, FALSE or 2
h <- -1  ## h>=0 leads to too many replicates

## open matlab server
Matlab$startServer()
matlab <- Matlab()
isOpen <- open(matlab)
if (!isOpen) throw("MATLAB server is not running: waited 30 seconds.")
evaluate(matlab, "rng('default');")
evaluate(matlab, "rng(1);")

## run multiple re-starts (j) with identical initial designs in
## one go

for(j in 76:100){

  ## calculate initial design
  set.seed(j)
  Xinit <- matrix(runif(ninit*dim), ncol=dim)
  Xinit <- (round(Xinit*19+1)-1)/19
  
  ## evaluate initial design via matlab server
  randMax <- 10000
  Y <- - ATO(Xinit, time=20, seed=1 + round(runif(nrow(Xinit))*randMax))
  
  ## fit initial hetGP model
  noiseControl <- list(g_min=1e-6, g_bounds=c(1e-6, 1), lowerDelta=log(1e-6))
  settings <- list(linkThetas="none", initStrategy="smoothed")
  mod <- mleHetGP(X=Xinit, Z=Y, lower=rep(0.01, dim), upper=rep(10,dim), 
                  maxit=10^5, noiseControl=noiseControl, settings=settings)
  X <- Xinit
  
  ## iterations of acquisition via BO
  for(i in 1:(ntot - nrow(X))) {
    
    ## choose between tricands and ordinary inner-optimization EI
    if(tri) {
      Xcand <- tricands(mod$X0) ##, best=which.min(predict(mod, mod$X0)$mean))
      Xcand <- round(19*Xcand + 1)
      Xcand <- Xcand[!duplicated(Xcand),] 
      Xcand <- (Xcand - 1)/19     
      opt <- crit_optim(mod, crit="crit_EI", h=h, Xcand=Xcand)
      if(tri >= 2) {
        opt2 <- crit_optim(mod, crit="crit_EI", h=h, control=list(Xstart=matrix(opt$par, nrow=1))) 
        opt2$par <- (round(opt2$par*19+1)-1)/19
        opt2$val <- crit_EI(opt2$par, mod) 
        if(opt2$val > opt$val) opt <- opt2
      }
    } else {
      opt <- crit_optim(mod, crit="crit_EI", h=h)
      opt$par <- (round(opt$par*19+1)-1)/19
    }

    ## ATO----
    randMax <- 10000
    Ynew <- - ATO(opt$par, time=20, seed=1 + round(runif(1)*randMax))
    ## ATO----
    
    ## augment data and updated hetGP fit
    X <- rbind(X, opt$par)
    Y <- c(Y, Ynew)
    mod <- mleHetGP(X=X, Z=Y, lower=rep(0.01,dim), upper=rep(10,dim), 
                  maxit=10^5, noiseControl=noiseControl, settings=settings)

    ## keep track of progress
    p <- predict(mod, X)
    ybest <- bov(p$mean)

    ## save intermediate results
    ## see ato_collect.R for visualization and comparison
    if(tri) {
      if(tri >= 2) save.image(paste0("ato_rdata/ato_seed_", j, "both.RData"))
      else save.image(paste0("ato_rdata/ato_seed_", j, "tricands.RData"))
    } else {
      save.image(paste0("ato_rdata/ato_seed_", j, "optim.RData"))
    }
    cat("rep:", j, " acq:", ninit+i, " ybest:", round(ybest[length(ybest)], 5), "\n")
  }
  
  Sys.time()                    
}

## close matlab server
close(matlab)
print(matlab)
