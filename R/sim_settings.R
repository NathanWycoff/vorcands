source("R/functions/function_lib.R")

#test <- TRUE
test <- FALSE

if (test) {
  print("TESTING ONLY!")
}
if (!exists('func')) stop("Need to define func before calling sim_settings.R")

## Sim settings.
if (test) {
  reps <- 2
} else {
  #reps <- 90
  print("Only 30 reps!")
  reps <- 30
}
#end <- 250

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

#ninit <- 100
#end <- 101
#reps <- 4

ncands <- min(5000,100*m)

root <- './sim_out/'
dir.create(root, showWarnings = FALSE)
sim_path <- paste(root,func,'/',sep='')
time_path <- paste(substr(sim_path,1,nchar(sim_path)-1),'_time/',sep='')
crits_path <- paste(substr(sim_path,1,nchar(sim_path)-1),'_crits/',sep='')
#opt_path <- paste(substr(sim_path,1,nchar(sim_path)-1),'_opt/',sep='')

# Prefix:
# 'pi' - Probability of Improvement.
# 'ei' - Expected Improvement.
# 'ts' - Thompson Sampling.

# Suffix:
# 's' - BFGS optim
# 't' - Tricands
# 'l' - LocoCands
# 'r' - LHS cands
#
# "Traditional" methods:
# 'bfgs' - BFGS + finite differencing
# 'nm' - Nelder-Meade Simplex
if (test) {
  competitors <- c(
    "gp.ei.opt",
    "gp.ei.voriRAS"
  )
  #competitors <- c('gp.ei.opt','hgp.ei.opt')
} else {
  competitors <- c(
    "nm",
    "bfgs",
    "gp.ei.opt",
    "gp.ei.lhs",
    #"gp.ei.corner",
    "gp.ei.voriRIS",
    "gp.ei.voriRLS",
    "gp.ei.voriRAS"
  )
  #competitors <- c('gp.ei.opt','hgp.ei.opt')
}

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
  hgp.ei.opt = 'hetGP-Opt',
  gp.ei.lhs = 'GP-LHS',
  gp.ei.corner = 'GP-Corner',
  gp.ei.voriRIS = 'GP-VWalk',
  gp.ei.voriRLS = 'GP-VProj',
  gp.ei.voriRAS = 'GP-VAlt'
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

#competitors <- c(
#                'gp.ei.voriRAS',
#                'gp.ei.vor1RAS',
#                'gp.ei.vor2RAS'
#                )

ltys <- sapply(competitors, function(comp) {
    csplit <- strsplit(comp,"\\.")[[1]]
    acq <- csplit[2]
    os <- csplit[3]

    if (comp %in% c("bfgs","nm")) {
        return(1)
    } else if (acq=='ei') {
      if (substr(os,1,3)=='vor') {
        return(2)
      } else {
        return(1)
      }
    } else {
      return(1)
        }
    #} else if (acq=='ts') {
    #    return(3)
    #} else if (acq=='pi') {
    #    return(4)
    #}
 })
 cols <- sapply(competitors, function(comp) {
    csplit <- strsplit(comp,"\\.")[[1]]
    acq <- csplit[2]
    os <- csplit[3]

    if (comp == "bfgs") {
        return('black')
    } else if (comp=="nm") {
        return("gray")
    } else if (os=='opt') {
        if (csplit[1]=='gp') {
            return('blue')
        } else {
            return('purple')
        }
    #} else if (os=='tri') {
    #    return('green')
    #} else if (substr(os,1,4)=='loco') {
    } else if (substr(os,1,3)=='tri') {
        args <- substr(os,4,nchar(os))
        if (args=='BS') return('red')
        if (args=='BL') return('purple')
    } else if (os=='lhs') {
        return('cyan')
    } else if (substr(os,1,3)=='vor') {
      args <- substr(os,4,nchar(os))
      #if (args=='1U') return('yellow')
      #if (args=='1R') return('orange')
      #if (args=='2U') return('gray')
      #if (args=='2R') return('blue')
      #if (args=='iU') return('green')
      if (args=='iRAS') return('red')
      if (args=='iRLS') return('pink')
      if (args=='iRIS') return('green')
      if (args=='1RAS') return('orange')
      if (args=='2RAS') return('green')
      if (grepl("fast",os,fixed=TRUE)) return('orange')
    } else if (os=='corner') {
      return('orange')
    }
 })
