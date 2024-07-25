source("R/functions/function_lib.R")

test <- TRUE
#test <- FALSE

if (test) {
  print("TESTING ONLY!")
}
if (!exists('func')) stop("Need to define func before calling sim_settings.R")


problems <- c("ackley10","levy10","rosen10","lunar","push","rover","pomp10","dacca")

if (func=='goldprice') {
  f <- goldprice
  m <- 2
  end <- 500
} else if (substr(func,1,nchar('ackley')) =='ackley') {
  f <- ackley
  m <- as.numeric(gsub('ackley','',func))
  ackley_u <- runif(m,-0.5,0.5)
  end <- 800
} else if (substr(func,1,nchar('edgy_ackley')) =='edgy_ackley') {
  f <- ackley
  m <- as.numeric(gsub('edgy_ackley','',func))
  ackley_u <- runif(m,-0.5,0.5)
  i <- sample(m,1)
  ackley_u[i] <- sample(c(0,1),1)
  end <- 800
} else if (substr(func,1,nchar('pompy_ackley')) =='pompy_ackley') {
  f <- ackley
  m <- as.numeric(gsub('pompy_ackley','',func))
  if (m!=10) stop("pompy ackley only works in 10D")
  ackley_u <- unname(c(read.csv('pomp10_opt.csv')))[[1]]
  end <- 800
} else if (substr(func,1,nchar('lunar_ackley')) =='lunar_ackley') {
  f <- ackley
  m <- 12
  ackley_u <- c(read.csv('tri_lunar.csv')[,1])
  end <- 800
} else if (substr(func,1,nchar('rosen')) =='rosen') {
  f <- rosen
  m <- as.numeric(gsub('rosen','',func))
  end <- 500
} else if (substr(func,1,nchar('levy')) =='levy') {
  f <- levy
  m <- as.numeric(gsub('levy','',func))
  end <- 500
} else if (func=='pomp10') {
    f <- pomp10_R
    m <- 10
    end <- 800
} else if (func=='pomp10log') {
    f <- pomp10_R_log
    m <- 10
    end <- 800
} else if (func=='dacca') {
    f <- dacca_R
    m <- 23
    #end <- 400
    end <- 800
} else if (func=='pde') {
    f <- pde_R
    m <- 4
    end <- 100
} else if (func=='lunar') {
    f <- lunar_R
    m <- 12
    end <- 400
} else if (func=='borehole') {
    f <- borehole
    m <- 8
    end <- 800
} else if (func=='push') {
    f <- push_R
    m <- 14
    end <- 1000
} else if (func=='rover') {
    f <- rover_R
    m <- 60
    #end <- 1000
    #end <- 600
    end <- 800
} else {
  stop("Unknown target function.")
}

ninit <- max(3*m,12)

if (test) {
  end <- 2*ninit
}

##ninit <- 100
#print("Small sim!")
#end <- 100
##reps <- 4

ncands <- min(5000,100*m)

root <- './sim_out/'
dir.create(root, showWarnings = FALSE)
sim_path <- paste(root,func,'/',sep='')
time_path <- paste(substr(sim_path,1,nchar(sim_path)-1),'_time/',sep='')
crits_path <- paste(substr(sim_path,1,nchar(sim_path)-1),'_crits/',sep='')
#opt_path <- paste(substr(sim_path,1,nchar(sim_path)-1),'_opt/',sep='')

#competitors <- c('gp.ei.vorsmRi',
#               'gp.ei.vorsmR1',
#               'gp.ei.vorsmR2',
#               'gp.ei.vorsmUi',
#               'gp.ei.vorsmU1',
#               'gp.ei.vorsmU2',
#               'gp.ei.vorsmRUi',
#               'gp.ei.vorsmRU1',
#               'gp.ei.vorsmRU2',
#               'gp.ei.vorUi',
#               'gp.ei.vorU1',
#               'gp.ei.vorU2',
#               'gp.ei.voralti',
#               'gp.ei.voralt1',
#               'gp.ei.voralt2')

#print("Only sobol!!!")
#competitors <- c('gp.ei.sobol')
# NOTE: The first method should be somthing like nm, gp.ei.lhs or anything except for gp.ei.tri because I rely on this having the normal length later.
competitors <- c('nm','bfgs','gp.ei.opt','gp.ei.lhs','gp.ei.sobol','gp.ei.tri', 'gp.ei.tr','gp.ei.voralti')
#competitors <- c('gp.ei.voralti')
#competitors <- c('gp.ei.tr')

#competitors <- c('gp.ei.vorsmRi',
#               'gp.ei.vorsmUi',
#               'gp.ei.vorsmRUi',
#               'gp.ei.vorUi',
#               'gp.ei.voralti')

#competitors <- c('gp.ei.vorsmR1',
#               'gp.ei.vorsmU1',
#               'gp.ei.vorsmRU1',
#               'gp.ei.vorU1',
#               'gp.ei.voralt1')

#competitors <- c('gp.ei.vorsmR2',
#               'gp.ei.vorsmU2',
#               'gp.ei.vorsmRU2',
#               'gp.ei.vorU2',
#               'gp.ei.voralt2')


if (func=='rover') {
    competitors <- competitors[competitors!='bfgs']
    competitors <- competitors[competitors!='gp.ei.opt']
} 

#pretty_method_names <- list(
#  nm = 'Nelder-Meade',
#  bfgs = 'BFGS',
#  gp.ei.opt = 'Optimized GP',
#  gp.ei.lhs = 'LHS Cands GP',
#  gp.ei.corner = 'Corner GP',
#  gp.ei.voriRIS = 'Vor Cands Walk GP',
#  gp.ei.voriRLS = 'Vor Cands Proj GP',
#  gp.ei.voriRAS = 'Vor Cands Alt GP'
#)

short_method_names <- list(
  nm = 'NM',
  bfgs = 'BFGS',
  gp.ei.opt = 'GP-Opt',
  gp.ei.lhs = 'GP-LHS',
  gp.ei.voralti = 'GP-Vor',
  gp.ei.tri = 'GP-Tri',
  gp.ei.tr = 'GP-TR'
)

pretty_sim_names <- list(
  ackley10 = 'Random Ackley (P=10)',
  levy10 = 'Levy (P=10)',
  rosen10 = 'Rosen (P=10)',
  lunar = 'Pygame Lunar (P=12)',
  push = 'Pygame Push (P=14)',
  rover = 'Pygame Rover (P=60)',
  pomp10log = 'Pomp10 Simulator (P=10)',
  dacca = 'Dacca Cholera Simulator (P=24)'
)

ltys <- sapply(competitors, function(comp) 1)
cols <- sapply(1:length(competitors), function(i) rainbow(length(competitors))[i])

names(cols) <- competitors
names(ltys) <- competitors
short_comp <- sapply(competitors, function(comp) {
    if (comp=='nm') return('nm')
    if (comp=='bfgs') return('bfgs')
    csplit <- strsplit(comp,"\\.")[[1]]
    os <- csplit[3]
    return(os)
})
print(short_comp)
