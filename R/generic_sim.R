tt <- Sys.time()
.libPaths( '/home/nwycoff_umass_edu/R/x86_64-pc-linux-gnu-library/4.4')
args <- commandArgs(trailingOnly=TRUE)
#func <- 'pomp10'
#func <- 'ackley10'
#seed <- 2
#func <- args[1]
#seed <- args[2]
ind <- args[1]

source("R/parset.R")

print("Func:")
print(func)
print("Seed:")
print(seed)
set.seed(seed)

source("R/sim_settings.R")
source("R/optim.R")

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

cntl <- lapply(competitors, function(x) 0)
names(cntl) <- competitors
prog <- lapply(competitors, function(x) rep(NA, end))
names(prog) <- competitors
crits <- lapply(competitors, function(x) rep(NA, end-ninit))
names(crits) <- competitors
opt <- lapply(competitors, function(x) rep(NA, m))
names(opt) <- competitors

#X <- matrix(runif(ninit*m), ncol=m) ## 
X <- randomLHS(ninit, m)

dir.create('./sim_inits/', showWarnings = FALSE)
write.csv(X, file = paste('./sim_inits/',func,'_',seed,'.csv',sep=''))

#apply(X,1,f)
#f((param.tab[,'mle1']-param.tab[,'box_min'])/(param.tab[,'box_max']-param.tab[,'box_min']))

# Dacca current opt: 8.229755

options(error=recover)

#competitors <- 'gp.ei.voralti'
#options(error=recover)
#for (i in 1:100) {
#    print("Oh no")
#}

times <- c()
for (comp in competitors) {
  print(comp)
  tt <- Sys.time()
  csplit <- strsplit(comp,"\\.")[[1]]
  sur <- csplit[1]

  ## Surrogate based methods
  if (sur%in% c('gp','hgp')) {
    acq <- csplit[2]
    ot <- csplit[3]

    options <- list(sur=sur, f=f, ninit=ninit, m=m, end=end, X=X, criteria=toupper(acq), ncands=ncands, cands=ot, ssm = 25)
    #print(sur)
    if (sur=='hgp') {
        options$pack <- 'hetGP'
    } else if (sur=='gp') {
        options$pack <- 'lagp'
    } else {print(sur);stop("Unknown sur")}

    options$cand_params <- list()
    if (substr(ot,1,3)=='vor') {
      if (grepl('1',ot,fixed=TRUE)) {
        options$cand_params$vor_norm = 'l1'
      } else if (grepl('2',ot,fixed=TRUE)) {
        options$cand_params$vor_norm = 'l2'
      } else if (grepl('i',ot,fixed=TRUE)) {
        options$cand_params$vor_norm = 'linf'
      } else {
        stop(paste("Unknown norm."))
      }
      #if (grepl('R',ot,fixed=TRUE)) {
      #  options$cand_params$vor_st = 'rect'
      #} else if (grepl('U',ot,fixed=TRUE)) {
      #  options$cand_params$vor_st = 'unif'
      #} else {
      #  stop("Unknown st.")
      #}
      #if (grepl("IS",ot,fixed=TRUE)) {
      #  options$cand_params$style <- 'rect'
      #} else if (grepl("LS",ot,fixed=TRUE)) {
      #  options$cand_params$style <- 'lhs'
      #} else if (grepl("AS",ot,fixed=TRUE)) {
      #  options$cand_params$style <- 'alt'
      #} else {
      #  stop("Unknown style.")
      #}
    }

    # Method specific settings
    if (acq=='ts') {
      istri <- ot=='tri'
      if (istri) options$close <- -1
      isvor <- ot=='vor'
      stopifnot(!isvor)
    }

    if (ot=='tri') {
        if (m >= 8) {
            options$end <- pmin(100, options$end)
        }
    } 

    os <- do.call(optim.surr, options)
    
    prog[[comp]] <- bov(os$y)
    crits[[comp]] <- os$crits
    opt[[comp]] <- os$X[which.min(os$y),]
    cntl[[comp]] <- cntl[[comp]] + os$cnt

  } else { ## Baseline methods
    if (comp=='bfgs') {
      ## try is necessary because optim with L-BFGS-B occasionally fails to start
      for(i in 1:ninit) {
        y <- c()
        to <- try(os <- optim(X[i,], f.prime, lower=0, upper=1, method="L-BFGS-B"), silent=TRUE)
        if(class(to) == "try-error") next;
        prog[[comp]] <- bov(y, end)
        break;
      }
    } else if (comp=='nm') {
      y <- c()
      os <- optim(X[1,], f.prime)
      prog[[comp]] <- bov(y, end)
    } else {
      stop(paste("Unknown competitor",comp))
    }
  }
  #times[comp] <- Sys.time() - tt
  times[comp] <- difftime(Sys.time(), tt, units='secs')
}

prog[['gp.ei.tri']] <- c(prog[['gp.ei.tri']], rep(NA, length(prog[[1]])-length(prog[['gp.ei.tri']])))
crits[['gp.ei.tri']] <- c(crits[['gp.ei.tri']], rep(NA, length(crits[[1]])-length(crits[['gp.ei.tri']])))

pdf <- data.frame(prog)
dir.create(sim_path, showWarnings = FALSE)
write.csv(pdf, paste(sim_path,func,'_',seed,'.csv',sep=''))

tdf <- data.frame(times)
dir.create(time_path, showWarnings = FALSE)
write.csv(tdf, paste(time_path,func,'_',seed,'.csv',sep=''))

cdf <- data.frame(crits)
dir.create(crits_path, showWarnings = FALSE)
write.csv(cdf, paste(crits_path,func,'_',seed,'.csv',sep=''))

#warnings()
