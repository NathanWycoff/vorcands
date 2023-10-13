## Opt-fringecands?

#library(Rfast)

#addturbo <- FALSE
addturbo <- TRUE

func <- commandArgs(trailingOnly=TRUE)[1]
source("R/sim_settings.R")

ltys[['TuRBO']] <- 1
ltys[['vor.TuRBO']] <- 1
cols[['TuRBO']] <- 'orange'
cols[['vor.TuRBO']] <- 'blue'

if (addturbo) {
  competitors <- c(competitors, 'TuRBO','vor.TuRBO')
}

### Progress
files <- list.files(sim_path, full.names = TRUE)

dfs <- list()
for (f in files) {
  dfs[[length(dfs)+1]] <-read.csv(f)
}

reps <- length(dfs)
progs <- lapply(competitors, function(i) matrix(0,end,reps))
names(progs) <- competitors
for (comp in competitors) {
  for (i in 1:reps) {
    progs[[comp]][,i] <- dfs[[i]][,comp]
  }
}

#dim(progs[['sur.ei.vor']])
#progs[['sur.ei.vor']][800,]
#min(progs[['sur.ei.triBS']][800,])

pdf(paste(func,'_norm.pdf',sep=''))
#ul <- max(sapply(progs, max))
ul <- max(sapply(progs, function(x) quantile(x, 0.8)))
#ul <- quantile(progs[['sur.ei.opt']],0.8)
ll <- min(sapply(progs, function(a) min(a, na.rm=T)))
lwd <- 5
plot(NA,NA,xlim=c(1,end), ylim=c(ll,ul), main = func, xlab='n: blackbox evaluations', ylab='median best observed value')

qu <- pmax(0,pmin(1,0.5+1.96*sqrt(0.25/reps)))
ql <- pmax(0,pmin(1,0.5-1.96*sqrt(0.25/reps)))

max_intervals <- Inf

for (comp in competitors) {
  print(comp)
  print(cols[[comp]])
  if (any(is.na(progs[[comp]]))) {
    print("Skipping nan competitor:")
    print(comp)
  } else {
    #points(1:end, rowMedians(progs[[comp]]), type = 'l', lty = ltys[comp], col = cols[comp], lwd = lwd)
    points(1:end, apply(progs[[comp]], 1, function(x) median(x)), type = 'l', lty = ltys[comp], col = cols[comp], lwd = 0.5*lwd)
    if (length(competitors) <= max_intervals) {
      points(1:end, apply(progs[[comp]], 1, function(x) quantile(x, ql)), type = 'l', lty = ltys[comp], col = cols[comp], lwd = 0.5*lwd)
      points(1:end, apply(progs[[comp]], 1, function(x) quantile(x, qu)), type = 'l', lty = ltys[comp], col = cols[comp], lwd = 0.5*lwd)
    }
  }
}
legend('topright', legend = competitors, lty=ltys,col=cols, lwd = lwd)
legend('bottomleft', legend = c(ql,qu), lty=ltys[competitors[[1]]],col=cols[competitors[[1]]], lwd = 0.5*lwd)
dev.off()

pdf(paste(func,'_box_norm.pdf',sep=''), width = 9, height = 5)
par(mfrow=c(1,2))
end <- dim(progs[[1]])[1]
half <- round(end/2)
labs <- rep(names(progs), rep(ncol(progs[[1]]),length(names(progs))))

mids <- do.call(c, lapply(progs, function(p) p[half,]))
boxplot(mids~labs, main = paste("Iteration",half), xlab = '', las = 2)

ends <- do.call(c, lapply(progs, function(p) p[end,]))
boxplot(ends~labs, main = paste("Iteration",end), xlab = '', las = 2)
dev.off()

### Time
#if (!addturbo) {
files <- list.files(time_path, full.names = TRUE)

dfs <- list()
for (f in files) {
  dfs[[length(dfs)+1]] <-read.csv(f)
}

names <- dfs[[1]][,1]
tdf <- sapply(dfs,function(x) x[,'times'])
rownames(tdf) <- names

pdf(paste(func,'_time_norm.pdf',sep=''))
boxplot(t(tdf))
dev.off()
#}
