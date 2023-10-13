source("R/optim.R")

## adapted from VLSE with vectorized and coded inputs, and normalized outputs

f <- function(X)
{

  if(is.null(nrow(X))) X <- matrix(X, nrow=1)
  X <- X*20 - 10 ## assumes coded inputs
  
  d <- ncol(X)
  w <- 1 + (X - 1)/4
  
  term1 <- (sin(pi*w[,1]))^2 
  term3 <- (w[,d]-1)^2 * (1+1*(sin(2*pi*w[,d]))^2)
  
  wi <- w[,1:(d-1),drop=FALSE]
  sum <- rowSums((wi-1)^2 * (1+10*(sin(pi*wi+1))^2))
  
  y <- term1 + sum + term3
  return(y/80)
  
}

rect <- matrix(c(0, 1), nrow=5, ncol=2, byrow=TRUE) 

start <- 12
end <- 75
reps <- 100

prog.tri <- prog <- matrix(NA, nrow=reps, ncol=end)
prog.dgp <- prog.dgptri <- prog.s <- prog.stri <- prog

for(r in 1:reps) {

  X <- matrix(runif(start*nrow(rect)), ncol=nrow(rect)) 
  Z <- f(X)
  
  out <- NULL

  for(i in 1:(end-start)) {
     out <- optim.step.tgp.new(f, X = X, Z = Z, rect = rect, prev = out, 
      improv=c(1,1), cands="lhs", verb = 0, NN=200, nug.p=0, gd=c(1e-5,2))
     X <- rbind(X, out$X)
     Z <- c(Z, f(out$X))
   }
   prog[r,] <- bov(Z)

   X <- X[1:start,]
   Z <- Z[1:start]
   out <- NULL

   for(i in 1:(end-start)) {
      out <- optim.step.tgp.new(f, X = X, Z = Z, rect = rect, prev = out, 
        improv=c(1,1), cands="tri", verb = 0, NN=200, nug.p=0, gd=c(1e-5,2))
      X <- rbind(X, out$X)
      Z <- c(Z, f(out$X))
    }
   prog.tri[r,] <- bov(Z)

   X <- X[1:start,]
   Z <- Z[1:start]

   for(i in 1:(end-start)) {
     out <- optim.step.tgp.simple(f, X = X, Z = Z, cands="lhs", verb=0, 
      NN=200, nug.p=0, gd=c(1e-5,2))
     X <- rbind(X, out)
     Z <- c(Z, f(out))
   }
   prog.s[r,] <- bov(Z)

   X <- X[1:start,]
   Z <- Z[1:start]

   for(i in 1:(end-start)) {
      out <- optim.step.tgp.simple(f, X = X, Z = Z, cands="tri", verb = 0, 
        NN=200, nug.p=0, gd=c(1e-5,2))
      X <- rbind(X, out)
      Z <- c(Z, f(out))
    }
   prog.stri[r,] <- bov(Z)

   X <- X[1:start,]
   Z <- Z[1:start]

   for(i in 1:(end-start)) {
     out <- optim.step.deepgp.simple(f, X = X, Z = Z, cands="lhs", verb=FALSE, 
      NN=200, true_g=1e-7)
     X <- rbind(X, out)
     Z <- c(Z, f(out))
   }
   prog.dgp[r,] <- bov(Z)

   X <- X[1:start,]
   Z <- Z[1:start]

   for(i in 1:(end-start)) {
      out <- optim.step.deepgp.simple(f, X = X, Z = Z, cands="tri", verb=FALSE, 
        NN=200, true_g=1e-7)
      X <- rbind(X, out)
      Z <- c(Z, f(out))
    }
   prog.dgptri[r,] <- bov(Z)

   save.image(file="levy.RData")
}

colMedians <- function(x, ...) { apply(x, 2, median, ...) }

## par(mfrow=c(1,3))

pdf("levy_prog.pdf", width=6, height=6)
plot(colMedians(prog, na.rm=TRUE), type="l", col=2, main="Levy 5d, no noise",
  ylab="median best observed value", xlab="n: blackbox evaluations", ylim=c(0.01, 0.2))
lines(colMedians(prog.tri, na.rm=TRUE), lty=2, lwd=2, col=2)
lines(colMedians(prog.s, na.rm=TRUE))
lines(colMedians(prog.stri, na.rm=TRUE), lty=2, lwd=2)
lines(colMedians(prog.dgp, na.rm=TRUE), col=3)
lines(colMedians(prog.dgptri, na.rm=TRUE), lty=2, lwd=2, col=3)
abline(v=start, lty=2)
legend("topright", c("TGP-lhs", "TGP-tri", "hyb-lhs", "tgp2-tri", "DGP-lhs", "DGP-tri"), 
  lty=c(1,2,1,2,1,2), col=c(1,1,2,2,3,3), lwd=c(1,2,1,2,1,2))
dev.off()

pdf("levy_50.pdf", width=4, height=6)  ## NOTE: ylim cuts off several zero-values
w <- 50
boxplot(prog.s[,w], prog.stri[,w], prog[,w], prog.tri[,w], prog.dgp[,w], prog.dgptri[,w], 
  xaxt="n",
  border=c(1,1,2,2,3,3), lwd=c(1,2,1,2,1,2), ylim=c(0, 0.15), 
  cex.axis=0.75, ylab="best observed value", xlab="", main=paste0(w, "th evaluation"))
axis(1, 1:6, c("TGP-lhs", "TGP-tri", "hyb-lhs", "hyb-tri", "DGP-lhs", "DGP-tri"), las=2)
dev.off()

pdf("levy_75.pdf", width=4, height=6)
w <- end
boxplot(prog.s[,w], prog.stri[,w], prog[,w], prog.tri[,w], prog.dgp[,w], prog.dgptri[,w], 
  xaxt="n",
  border=c(1,1,2,2,3,3), lwd=c(1,2,1,2,1,2), ylim=c(0, 0.15),
  cex.axis=0.75, ylab="best observed value", xlab="", main="final evaluation")
axis(1, 1:6, c("TGP-lhs", "TGP-tri", "hyb-lhs", "hyb-tri", "DGP-lhs", "DGP-tri"), las=2)
dev.off()