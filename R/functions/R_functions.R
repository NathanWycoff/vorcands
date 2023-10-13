## f:
##
## Goldstein-Price test function
goldprice <- function(X)
 {
  if(is.null(nrow(X))) X <- matrix(X, nrow=1)
  m <- 8.6928
  s <- 2.4269
  x1 <- 4*X[,1] - 2
  x2 <- 4*X[,2] - 2
  a <- 1 + (x1 + x2 + 1)^2 * 
    (19 - 14*x1 + 3*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2)
  b <- 30 + (2*x1 - 3*x2)^2 * 
    (18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2)
  f <- log(a*b)
  f <- (f - m)/s
  return(f) 
 }

 

ackley <- function(xx, a=20, b=0.2, c=2*pi)
{
  ##########################################################################
  #
  # ACKLEY FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., xd)
  # a = constant (optional), with default value 20
  # b = constant (optional), with default value 0.2
  # c = constant (optional), with default value 2*pi
  #
  ##########################################################################
	
  if ('matrix' %in% class(xx)) {
    d <- ncol(xx)
  } else {
    d <- length(xx)
  }

  # Optimum in the middle of the space is cheating.
  xx <- xx-ackley_u
  
  # Scale to approximately [0,1].
  xx <- 2*32.768*(xx-0.5)

  if ('matrix' %in% class(xx)) {
    sum1 <- rowSums(xx^2)
    sum2 <- rowSums(cos(c*xx))
  } else {
    sum1 <- sum(xx^2)
    sum2 <- sum(cos(c*xx))
  }

  term1 <- -a * exp(-b*sqrt(sum1/d))
  term2 <- -exp(sum2/d)

  y <- term1 + term2 + a + exp(1)
  return(y)
}

borehole <- function(xx)
{
  ##########################################################################
  #
  # BOREHOLE FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # OUTPUT AND INPUT:
  #
  # y  = water flow rate
  # xx = c(rw, r, Tu, Hu, Tl, Hl, L, Kw)
  #
  ##########################################################################

  lb <- c(0.05, 100, 63070, 990, 63.1, 700, 1120, 9855)
  ub <- c(0.15, 50000, 115600, 1110, 116, 820, 1680, 12045)

  if ('matrix' %in% class(xx)) {
    xx <- t(t(xx) * (ub-lb) + lb)

    rw <- xx[,1]
    r  <- xx[,2]
    Tu <- xx[,3]
    Hu <- xx[,4]
    Tl <- xx[,5]
    Hl <- xx[,6]
    L  <- xx[,7]
    Kw <- xx[,8]
  } else {
    xx <- xx * (ub-lb) + lb

    rw <- xx[1]
    r  <- xx[2]
    Tu <- xx[3]
    Hu <- xx[4]
    Tl <- xx[5]
    Hl <- xx[6]
    L  <- xx[7]
    Kw <- xx[8]
  }

  frac1 <- 2 * pi * Tu * (Hu-Hl)

  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)

  y <- frac1 / frac2
  return(y)
}

rosen <- function(xx, small = FALSE)
{
  ##########################################################################
  #
  # ROSENBROCK FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUT:
  #
  # xx = c(x1, x2, ..., xd)
  #
  ##########################################################################
  xismat <- 'matrix' %in% class(xx)

  d <- ifelse(xismat, ncol(xx), length(xx))
  if (small) {
    lb <- -2.048
    ub <- 2.048
  } else {
    lb <- -5
    ub <- 10
  }

  xx <- xx * (ub-lb) + lb
  if (xismat) {

    xi <- xx[,1:(d-1),drop=F]
    xnext <- xx[,2:d,drop=F]
    y <- rowSums(100*(xnext-xi^2)^2 + (xi-1)^2)
  } else {
    xx <- xx * (ub-lb) + lb

    xi <- xx[1:(d-1)]
    xnext <- xx[2:d]
    y <- sum(100*(xnext-xi^2)^2 + (xi-1)^2)
  }

  return(y)
}


#N <- 10
#P <- 10
#xx <- matrix(runif(N*P),ncol=P)

levy <- function(xx)
{
  ##########################################################################
  #
  # LEVY FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUT:
  #
  # xx = c(x1, x2, ..., xd)
  #
  ##########################################################################

  xismat <- 'matrix' %in% class(xx)

  lb <- -10
  ub <- 10
  xx <- xx * (ub-lb) + lb

  d <- ifelse(xismat, ncol(xx), length(xx))
  
  if (xismat) { 
    w <- 1 + (xx - 1)/4

    term1 <- (sin(pi*w[,1,drop=F]))^2 
    term3 <- (w[,d,drop=F]-1)^2 * (1+1*(sin(2*pi*w[,d,drop=F]))^2)

    wi <- w[,1:(d-1),drop=F]
    ss <- rowSums((wi-1)^2 * (1+10*(sin(pi*wi+1))^2))
  } else {
    w <- 1 + (xx - 1)/4

    term1 <- (sin(pi*w[1]))^2 
    term3 <- (w[d]-1)^2 * (1+1*(sin(2*pi*w[d]))^2)

    wi <- w[1:(d-1)]
    ss <- sum((wi-1)^2 * (1+10*(sin(pi*wi+1))^2))
  }

  y <- term1 + ss + term3
  return(y)
}
