tt <- Sys.time()

args <- commandArgs(trailingOnly=TRUE)
#func <- 'pomp10'
#func <- 'ackley10'
func <- 'ackley2'
seed <- 1
set.seed(seed)

source("R/sim_settings.R")
source("R/optim.R")

#end <- 400
end <- 100

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

times <- c()
ret <- list()
for (comp in c("hgp.ei.opt", "gp.ei.opt")) {
    tt <- Sys.time()
    csplit <- strsplit(comp,"\\.")[[1]]
    sur <- csplit[1]

    ## Surrogate based methods
    acq <- csplit[2]
    os <- csplit[3]

    options <- list(sur=sur, f=f, ninit=ninit, m=m, end=end, X=X, criteria=toupper(acq), ncands=ncands, cands=os, ssm = 25)
    print(sur)
    if (sur=='hgp') {
        options$pack <- 'hetGP'
    } else if (sur=='gp') {
        options$pack <- 'lagp'
    } else {print(sur);stop("Unknown sur")}

    options$cand_params <- list()
    if (substr(os,1,3)=='vor') {
      if (grepl('1',os,fixed=TRUE)) {
        options$cand_params$vor_norm = 'l1'
      } else if (grepl('2',os,fixed=TRUE)) {
        options$cand_params$vor_norm = 'l2'
      } else if (grepl('i',os,fixed=TRUE)) {
        options$cand_params$vor_norm = 'linf'
      } else {
        stop("Unknown norm.")
      }
      if (grepl('R',os,fixed=TRUE)) {
        options$cand_params$vor_st = 'rect'
      } else if (grepl('U',os,fixed=TRUE)) {
        options$cand_params$vor_st = 'unif'
      } else {
        stop("Unknown st.")
      }
      if (grepl("IS",os,fixed=TRUE)) {
        options$cand_params$style <- 'rect'
      } else if (grepl("LS",os,fixed=TRUE)) {
        options$cand_params$style <- 'lhs'
      } else if (grepl("AS",os,fixed=TRUE)) {
        options$cand_params$style <- 'alt'
      } else {
        stop("Unknown style.")
      }
    }

    # Method specific settings
    if (acq=='ts') {
      istri <- os=='tri'
      if (istri) options$close <- -1
      isvor <- os=='vor'
      stopifnot(!isvor)
    }

    #options$pack <- 'lagp'
    #options$pack <- 'hetGP'

    ret[[comp]] <- do.call(optim.surr, options)
    times[comp] <- Sys.time() - tt
}

pdf("temp.pdf")
plot(NA,NA,xlim=c(0,end), ylim = c(0,22))
for (comp in competitors) {
    points(ret[[comp]]$y, col = cols[[comp]])
}
dev.off()

