## First run ato.R with tri alternately 
## set to TRUE and FALSE

## read in tricands runs (tri <- TRUE)
tcf <- list.files("ato_rdata", pattern="tricands")
tprog <- matrix(NA, nrow=length(tcf), ncol=300)
for(fi in 1:length(tcf)) {
  load(paste0("ato_rdata/", tcf[fi]))
  tprog[fi,1:length(ybest)] <- ybest
}

## w <- c(1:8,10:82, 84:99)
w <- 1:100
tprog <- tprog[w,]

## plot idividual tricands runs
matplot(t(tprog), col=1, type="l", lty=3, xlim=c(0,300), 
  ylim=c(-2.75, -0.5), xlab="n: blackbox evaluations",
  ylab="model predicted best observed value",
  main="Assemble to Order, input-dependent noise")

## summarize middle
## tavg <- colMeans(tprog, na.rm=TRUE)
tavg <- apply(tprog, 2, median, na.rm=TRUE)
lines(tavg, lwd=3)

## then read in inner-opt runs (tri <- FALSE)
ocf <- list.files("ato_rdata", pattern="optim")
oprog <- matrix(NA, nrow=length(ocf), ncol=300)
for(fi in 1:length(ocf)) {
  load(paste0("ato_rdata/", ocf[fi]))
  oprog[fi,1:length(ybest)] <- ybest
}
oprog <- oprog[w,]

## overlay individual optim runs
matlines(t(oprog), col=2, lty=3)

## summarize middle
## oavg <- colMeans(oprog, na.rm=TRUE)
oavg <- apply(oprog, 2, median, na.rm=TRUE)
lines(oavg, lwd=3, col=2)

## then read in inner-opt runs (tri <- FALSE)
bcf <- list.files("ato_rdata", pattern="both")
bprog <- matrix(NA, nrow=length(bcf), ncol=300)
for(fi in 1:length(bcf)) {
  load(paste0("ato_rdata/", bcf[fi]))
  bprog[fi,1:length(ybest)] <- ybest
}
bprog <- bprog[w,]

## overlay individual optim runs
matlines(t(bprog), col=3, lty=3)

## summarize middle
## bavg <- colMeans(bprog, na.rm=TRUE)
bavg <- apply(bprog, 2, median, na.rm=TRUE)
lines(bavg, lwd=3, col=3)

## initialization
abline(v=ninit, col="gray", lty=2)

## add legend
legend("topright", c("tricands", "bfgs", "hybrid", "each", "median"),
   lty=c(1,1,1,2,1), lwd=c(1,1,1,1,2), col=c(1,2,3,1,1), bty="n")

## figures for the paper

pdf("ato_prog.pdf", width=6, height=6)
plot(oavg, type="l", col=1, lwd=1, main="Assemble to Order, input-dependent noise",
  ylab="median (estimated) best observed value", xlab="n: blackbox evaluations", ylim=c(-2.75, -0.5))
lines(tavg, lty=2, col=1, lwd=2)
lines(bavg, lty=3, col=2, lwd=2)
abline(v=80, lty=2)
legend("topright", c("EI", "EI-tri", "EI-hyb"), lty=1:3, col=c(1,1,2), lwd=c(1,2,2), bty="n")
dev.off()

pdf("ato_150.pdf", width=4, height=6) 
w <- 150
boxplot(oprog[,w], tprog[,w], bprog[,w],
  xaxt="n",
  border=c(1,1,2), ylim=c(-2.75, -2), lwd=c(1,2,2),
  cex.axis=0.75, ylab="(estimated) best observed value", xlab="", main=paste0(w, "th evaluation"))
axis(1, 1:3, c("EI", "EI-tri", "EI-hyb"), las=2)
dev.off()

pdf("ato_300.pdf", width=4, height=6)
w <- 300
boxplot(oprog[,w], tprog[,w], bprog[,w],
  xaxt="n",
  border=c(1,1,2), ylim=c(-2.75, -2), lwd=c(1,2,2), 
  cex.axis=0.75, ylab="(estimated) best observed value", xlab="", main="final evaluation")
axis(1, 1:3, c("EI", "EI-tri", "EI-hyb"), las=2)
dev.off()

## testing

wilcox.test(tprog[,w], oprog[,w], paired=TRUE, alternative="less")
