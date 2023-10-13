source("optim.R")

library(DiceOptim)
hartman6vec <- function(X)
 {
  if(is.matrix(X)) {
    y <- apply(X, 1, hartman6)
  } else y <- hartman6(X)
  return(y)
 }

hartman6.prime <- function(x)
 { 
  if(any(x < 0) || any(x > 1)) ynew <- Inf
  else ynew <- hartman6(x)
  y <<- c(y, ynew)
  return(ynew)
 }

reps <- 100
ninit <- 12
end <- 50
cnt.pi <- cnt.pit <- cnt.pir <- cnt.ei <- cnt.eit <- cnt.eir <- cnt.ts <- cnt.tst <- 0

prog.ei <- prog.eit <- prog.eir <- prog.pi <- prog.pit <- prog.pir <- prog.bfgs <- prog.nm <- 
  prog.ts <- prog.tst <- matrix(NA, nrow=reps, ncol=end)
for(r in 1:reps) {

  ## surrogate-based comparators
  X <- matrix(runif(ninit*6), ncol=6) ## randomLHS(ninit, 6)
  os <- optim.surr(hartman6vec, ninit, 6, end, X, "EI")
  prog.ei[r,] <- bov(os$y) 
  cnt.ei <- cnt.ei + os$cnt
  os <- optim.surr(hartman6vec, ninit, 6, end, X, "EI", cands="tri")
  prog.eit[r,] <- bov(os$y) 
  cnt.eit <- cnt.eit + os$cnt
  os <- optim.surr(hartman6vec, ninit, 6, end, X, "EI", cands="lhs")
  prog.eir[r,] <- bov(os$y) 
  cnt.eir <- cnt.eir + os$cnt
  os <- optim.surr(hartman6vec, ninit, 6, end, X, "PI")
  prog.pi[r,] <- bov(os$y) 
  cnt.pi <- cnt.pi + os$cnt
  os <- optim.surr(hartman6vec, ninit, 6, end, X, "PI", cands="tri")
  prog.pit[r,] <- bov(os$y)
  cnt.pit <- cnt.pit + os$cnt
  os <- optim.surr(hartman6vec, ninit, 6, end, X, "PI", cands="lhs")
  prog.pir[r,] <- bov(os$y)
  cnt.pir <- cnt.pir + os$cnt
  os <- optim.surr(hartman6vec, ninit, 6, end, X, "TS")
  prog.ts[r,] <- bov(os$y) 
  cnt.ts <- cnt.ts + os$cnt
  os <- optim.surr(hartman6vec, ninit, 6, end, X, "TS", close=-1)
  prog.tst[r,] <- bov(os$y)
  cnt.tst <- cnt.tst + os$cnt

  ## try is necessary because optim with L-BFGS-B occasionally fails to start
  for(i in 1:ninit) {
    y <- c()
    to <- try(os <- optim(X[i,], hartman6.prime, lower=0, upper=1, method="L-BFGS-B"), silent=TRUE)
    if(class(to) == "try-error") next;
    prog.bfgs[r,] <- bov(y, end)
    break;
  }
  
  y <- c()
  os <- optim(X[1,], hartman6.prime)
  prog.nm[r,] <- bov(y, end)
}

colMedians <- function(x) { apply(x, 2, median) }

pdf("hartmann6_prog.pdf", width=6, height=6)
plot(colMedians(prog.ei), col=1, type="l", main="Hartmann 6, no noise", ylim=c(-3.25, 0),
  xlab="n: blackbox evaluations", ylab="median best observed value")
lines(colMedians(prog.eit), col=1, lty=2, lwd=2)
lines(colMedians(prog.eir), col=1, lty=3)
lines(colMedians(prog.bfgs), col=2)
lines(colMedians(prog.nm), col=2, lty=2)
lines(colMedians(prog.ts), col=4 ,lty=3)
lines(colMedians(prog.tst), col=4, lty=2, lwd=2)
abline(v=ninit, lty=2)
legend("bottomleft", lty=2, col=1, legend="initial")
legend("top", col=2, lty=c(2,1), legend=c("Nelder-Mead", "BFGS"), bty="n")
legend("topright", c("EI", "EI-tri", "EI-lhs", "TS-tri", "TS-lhs"), 
  col=c(1, 1, 1, 4, 4), lty=c(1,2,3,2,3), lwd=c(1,2,1,2,1), bty="n")
dev.off()

pdf("hartmann6_30.pdf", width=4, height=6)
w <- 30
boxplot(prog.ei[,w], prog.eit[,w], prog.eir[,w],
  prog.tst[,w], prog.ts[,w], prog.nm[,w], prog.bfgs[,w], 
  ylim=c(-3.25, -0.75), border=c(1, 1, 1, 4, 4, 2, 2), lwd=c(1,2,1,2,1,1,1), xaxt="n",
  ylab="best observed value", xlab="", main="30th evaluation", cex.lab=0.75)
axis(1, 1:7, c("EI", "EI-tri", "EI-lhs", "TS-tri", "TS-lhs", "NM", "BFGS"), las=2)
dev.off()

pdf("hartmann6_50.pdf", width=4, height=6)
w <- end
boxplot(prog.ei[,end], prog.eit[,end], prog.eir[,end],
  prog.tst[,end], prog.ts[,end], prog.nm[,end], prog.bfgs[,end], 
  ylim=c(-3.25, 0), border=c(1, 1, 1, 4, 4, 2, 2), lwd=c(1,2,1,2,1,1,1), xaxt="n", 
  ylab="best observed value", xlab="", main="final evaluation")
text(x=3.1, y=-0.1, labels="# of candidates/criteria evals", cex=0.75)
text(x=c(1:5), y=rep(-0.25,5), labels=round(c(cnt.ei, cnt.eit, cnt.eir,
  cnt.tst, cnt.ts)/reps), cex=0.65)
axis(1, 1:7, c("EI", "EI-tri", "EI-lhs", "TS-tri", "TS-lhs", "NM", "BFGS"), las=2, cex=0.75)
dev.off()

save.image("hartmann6.RData")