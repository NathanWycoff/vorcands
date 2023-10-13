library(pomp)

## load dacca cholera model from pomp package
po <- dacca()

## incase you wanna see what we working with data-wise
# plot(po)

## roughly the best parameter set they've found for dacca
baseline_params <- coef(po)

## wrapper to return marginal log-likelihood values from dacca model
dacca_bo <- function(filter, params, Np = 2000){
  ## set fixed parameter values (??)
  params['rho'] = 0
  params['Y_0'] = 0
  params['delta'] = .02
  params['clin'] = 1
  params['alpha'] = 1
  return(logLik(pfilter(filter, params = params, Np = Np)))
}

## upper and lower bound infos for the 28 dacca parameters

param.tab <- as.data.frame(read.table(row.names=1,header=TRUE,text="
                  mle1 box_min box_max
gamma      20.800000000   10.00   40.00
eps        19.100000000    0.20   30.00
deltaI      0.060000000    0.03    0.60
beta_trend -0.004980000   -0.01    0.00
logbeta1   0.747000000   -4.00    4.00
logbeta2   6.380000000    0.00    8.00
logbeta3  -3.440000000   -4.00    4.00
logbeta4   4.230000000    0.00    8.00
logbeta5   3.330000000    0.00    8.00
logbeta6   4.550000000    0.00    8.00
logomega1 -1.692819521  -10.00    0.00
logomega2 -2.543383579  -10.00    0.00
logomega3 -2.840439389  -10.00    0.00
logomega4 -4.691817993  -10.00    0.00
logomega5 -8.477972478  -10.00    0.00
logomega6 -4.390058806  -10.00    0.00
sd_beta     3.130000000    1.00    5.00
tau         0.230000000    0.10    0.50
S_0         0.621000000    0.00    1.00
I_0         0.378000000    0.00    1.00
R1_0        0.000843000    0.00    1.00
R2_0        0.000972000    0.00    1.00
R3_0        0.000000116    0.00    1.00
"))

## function to generate random parameter set, with All drawn from uniform
dacca.rparam.unif <- function(hyperparams, ...)
{
  r <- runif(
    n=length(hyperparams$box_min),
    min=hyperparams$box_min,
    max=hyperparams$box_max
  )
  names(r) <- row.names(hyperparams)
  return(r)
}

## test on best value
#dacca_bo(filter = po, params = baseline_params)

## test on randomly generated boys
#dacca_bo(filter = po, params = dacca.rparam.unif(param.tab))

# Good value is about -3750

# [0,1] valued
dacca_R_scalar <- function(xx) {
    reseed <- sample(1e5,1)
    set.seed(123)
    xxab <- (param.tab$box_max-param.tab$box_min) * xx + param.tab$box_min
    names(xxab) <- rownames(param.tab)
    yr <- log(-dacca_bo(filter = po, params = xxab))
    set.seed(reseed)
    return(yr)
}


dacca_R <- function(X) {
  if (length(dim(X)) > 1) {
    return(apply(X, 1, dacca_R_scalar))
  } else {
    return(dacca_R_scalar(X))
  }
}