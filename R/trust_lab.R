tt <- Sys.time()

args <- commandArgs(trailingOnly=TRUE)
func <- 'ackley2'
seed <- 1
set.seed(seed)

source("R/sim_settings.R")
source("R/optim.R")

#end <- 400
end <- 40

print("Seed:")
print(seed)

## Stuff to look closer at:
#ncands <- max(50,100*m)

## Indicator formulation of box-constrained optimization.
f.prime <- function(x)
 { 
  if(any(x < 0) || any(x > 1)) ynew <- Inf
  else ynew <- f(x)
  y <<- c(y, ynew)
  return(ynew)
 }

#X <- matrix(runif(ninit*m), ncol=m) ## 
X <- randomLHS(ninit, m)

comp <- 'gp.ei.tr'

tt <- Sys.time()
csplit <- strsplit(comp,"\\.")[[1]]
sur <- csplit[1]

## Surrogate based methods
acq <- csplit[2]
os <- csplit[3]

options <- list(sur=sur, f=f, ninit=ninit, m=m, end=end, X=X, criteria=toupper(acq), ncands=ncands, cands=os, ssm = 25)
source("R/optim.R")
ret <- do.call(optim.surr, options)
