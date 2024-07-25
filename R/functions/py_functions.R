library(reticulate)
use_python("/home/r41x461/.conda/envs/nebs05/bin/python", required = T)
#use_virtualenv("./.venv/", required = T)
use_condaenv("nebs05", required = T)
source_python("./python/bnn_fun.py")
#source_python("./R/functions/push.py")
source_python("./python/push.py")
source_python("./python/rover.py")

## If these are max problems, negation should be on the python side.

lunar_R <- function(X) {
  if (length(dim(X)) > 1) {
    return(apply(X, 1, lunar))#DRY1
  } else {
    return(lunar(X))#DRY1
  }
}

pde_R <- function(X) {
  if (length(dim(X)) > 1) {
    return(apply(X, 1, Ppde))
  } else {
    return(Ppde(X))
  }
}

push_R <- function(X) {
  if (length(dim(X)) > 1) {
    return(apply(X, 1, push_py2R))#DRY1
  } else {
    return(push_py2R(X))#DRY1
  }
}

rover_R <- function(X) {
  if (length(dim(X)) > 1) {
    return(apply(X, 1, rover_py2R))#DRY1
  } else {
    return(rover_py2R(X))#DRY1
  }
}
