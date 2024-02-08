# Voronoi Candidates for Bayesian Optimization

This R repository contains code accompanying the paper Voronoi Candidates for Bayesian Optimization, available at https://arxiv.org/abs/2402.04922

## License Notice
Please note that this directory contains a modified project not authored by the repository owners:
1) RANN.Linf

RANN.Linf is version of the RANN library modified to use the Linf norm.

Please see the original versions of those libraries for licensing information. 
The license information for Vorcands does not cover those libraries.

## To get started using Voronoi Candidates in your projects:
1) Clone this repo
2) Install the R package RANN.Linf, provided here, maybe by setting your working directory to be the root of this repo, and then running (in your command line) R CMD INSTALL RANN.Linf
3) Install CRAN R packages via running (in R) install.packages(c("Rcpp","lhs"))

Then, you can get started using candidates by running something like the below code:

<R>
    
    source("R/vornoi_cands.R")
    N <- 10
    P <- 2
    
    X <- matrix(runif(N*P), ncol = P)
    
    Xcands <- vorwalkcands(X, st = 'unif')$Xs # For random vorwalk
    #Xcands <- vorwalkcands(X, st = 'rect')$Xs # For axis-aligned vorwalk
    #Xcands <- vorwalkcands(X, st = 'lhs')$Xs # For vorproj
    
    plot(NA,NA,xlim=c(0,1),ylim=c(0,1))
    points(X[,1],X[,2])
    points(Xcands[,1],Xcands[,2], col = 'red')

In practice, we recommend alternating between "rect" and "lhs" when dealing with high dimensional problems (but "unif" looks much cooler for a 2D plot).

## To reproduce the results in the paper:

### Installation Instructions:
Since the numerical experiments use functions written in Python, there are many more installation requirements to reproduce the results than simply to use the methodology.

- Install RANN.Linf: 
    R CMD INSTALL ./RANN.Linf
- Make a python virtual environment named ".venv" (R is expecting this name to be used, and for it to exist in the root of this directory.)
    maybe using "python3 -m venv .venv"
- After loading this venv, install python packages:
    pip install pygame gym[box2d] botorch rpy2 matplotlib pandas
- Install CRAN R packages 
    install.packages(c("stringr","hetGP","akima","laGP","Rcpp","pomp","reticulate","lhs"))

### Reproducing Illustrative Figures
To reproduce the illustrative figures in the paper, simply run R files which begin with the string "figure", e.g. Rscript R/figure_1.R

### Reproducing Numerical Experiments
To reproduce the numerical experiments, first  install gnu parallel on your system. On Ubuntu, this can be done via apt install parallel.

Then, give permission to the script run_experiments.sh to run and execute it.
