## Opt-fringecands?
#library(Rfast)

.libPaths( '/home/nwycoff_umass_edu/R/x86_64-pc-linux-gnu-library/4.4')

func <- commandArgs(trailingOnly=TRUE)[1]
source("R/sim_settings.R")

if ('gp.ei.corner' %in% competitors) {
    competitors <- competitors[-which(competitors=='gp.ei.corner')]
}

#pretty_comp <- sapply(competitors, function(x) pretty_method_names[[x]])
#short_comp <- sapply(competitors, function(x) short_method_names[[x]])

### Progress
files <- list.files(sim_path, full.names = TRUE)

dfs <- list()
for (f in files) {
  dfs[[length(dfs)+1]] <-read.csv(f)
}

reps <- length(dfs)
progs <- lapply(competitors, function(i) matrix(NA,end,reps))
names(progs) <- competitors
#torem <- c('bfgs')
torem <- c()
for (comp in competitors) {
  if ((comp %in% colnames(dfs[[1]]))) {
    for (i in 1:reps) {
      progs[[comp]][,i] <- dfs[[i]][,comp]
    }
  } else {
    torem <- c(torem, comp)
  }
}

competitors <- competitors[!(competitors %in% torem)]
progs <- progs[!(names(progs) %in% torem)]
#pretty_comp <- pretty_comp[!(names(pretty_comp) %in% torem)]
short_comp <- short_comp[!(names(short_comp) %in% torem)]
cols <- cols[!(names(cols) %in% torem)]
ltys <- ltys[!(names(ltys) %in% torem)]

if (func=='dacca') {
  for (comp in competitors) {
    progs[[comp]] <- log10(exp(progs[[comp]]) - 3750)
  }
} else if (func=='pomp10log') {
  for (comp in competitors) {
    progs[[comp]] <- log10(10^(progs[[comp]]) - 630)
  }
}

#dim(progs[['sur.ei.vor']])
#progs[['sur.ei.vor']][800,]
#min(progs[['sur.ei.triBS']][800,])

#### Main progress pdf.
pdf(paste(func,'.pdf',sep=''), width = 4, height = 4)
par(mar=c(3.0,3,1,0.1)+0.1)
par(mgp=c(1.5,0.6,0)+0.1)
##ul <- max(sapply(progs, max))
#ul <- max(sapply(progs, function(x) quantile(x, 0.8)))
##ul <- quantile(progs[['sur.ei.opt']],0.8)
#ll <- min(sapply(progs, function(a) min(a, na.rm=T)))

qu <- pmax(0,pmin(1,1.1*(0.5+1.96*sqrt(0.25/reps))))
ql <- pmax(0,pmin(1,0.9*(0.5-1.96*sqrt(0.25/reps))))

budget <- nrow(progs[[1]])
#start_at = budget / 2
pcomp <- grep("^gp.*", competitors, value = TRUE)
ul_end <- max(sapply(pcomp, function(x) quantile(progs[[x]][budget,], qu, na.rm = T)), na.rm = T)
ll_end <- min(sapply(pcomp, function(x) quantile(progs[[x]][budget,], ql, na.rm = T)), na.rm = T)

if (substr(func,1,nchar('rosen'))=='rosen') {
    sa <- 100
} else {
    sa <- ninit
}

ul_start <- max(sapply(pcomp, function(x) quantile(progs[[x]][sa,], qu, na.rm = T)), na.rm = T)
ll_start <- min(sapply(pcomp, function(x) quantile(progs[[x]][sa,], ql, na.rm = T)), na.rm = T)
ul <- max(c(ul_end,ul_start))
ll <- min(c(ll_end,ll_start))

lwd <- 5
if (func %in% c("pomp10log","dacca")) {
  ylab <- 'Log10 Optimality Gap'
} else {
  ylab <- 'Best Observed Value'
}
plot(NA,NA,xlim=c(1,end), ylim=c(ll,ul), main = pretty_sim_names[[func]], xlab='n: blackbox evaluations', ylab=ylab, font.lab = 2)


max_intervals <- Inf

for (comp in competitors) {
  if (all(is.na(progs[[comp]]))) {
    print("Skipping nan competitor:")
    print(comp)
  } else {
    if (length(competitors) <= max_intervals) {
      prog <- progs[[comp]]
      prog <- prog[apply(prog, 1, function(x) !all(is.na(x))),]
      il = apply(prog, 1, function(x) quantile(x, ql, na.rm = T))
      iu = apply(prog, 1, function(x) quantile(x, qu, na.rm = T))
      polygon(c(1:length(il),length(il):1), c(il,rev(iu)), col = adjustcolor(cols[comp], alpha.f = 0.2), border = NA)
    }
  }
}

for (comp in competitors) {
  print(comp)
  print(cols[[comp]])
  if (all(is.na(progs[[comp]]))) {
    print("Skipping nan competitor:")
    print(comp)
  } else {
    #points(1:end, rowMedians(progs[[comp]]), type = 'l', lty = ltys[comp], col = cols[comp], lwd = lwd)
    points(1:end, apply(progs[[comp]], 1, function(x) median(x, na.rm = T)), type = 'l', col = cols[comp], lwd = 0.5*lwd)
    if (length(competitors) <= max_intervals) {
      il = apply(progs[[comp]], 1, function(x) quantile(x, ql, na.rm = T))
      iu = apply(progs[[comp]], 1, function(x) quantile(x, qu, na.rm = T))
      points(1:end, il, type = 'l', lty = 'dotted', col = cols[comp], lwd = 0.5*lwd)
      points(1:end, iu, type = 'l', lty = 'dotted', col = cols[comp], lwd = 0.5*lwd)
    }
  }
}
dev.off()

#pdf("legend.pdf", width = 6, height = 2)
if (func=='ackley10') {
    pdf("legend.pdf", width = 6, height = 4)

    plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    #legend('bottomleft', legend = pretty_comp, lty=ltys,col=cols, lwd = lwd, bg = 'white', cex = 0.5, ncol = 4)
    legend('bottomleft', legend = short_comp, lty=ltys,col=cols, lwd = lwd, bg = 'white', cex = 0.4, ncol = 7)
    dev.off()
}

pdf(paste(func,'_box.pdf',sep=''), width = 7, height = 4)
par(mar=c(4.7,2,1,0)+0.1)
par(mgp=c(1,0.6,0)+0.1)
par(mfrow=c(1,2))
end <- dim(progs[[1]])[1]
half <- round(end/2)
labs <- unname(rep(sapply(names(progs), function(x) short_method_names[[x]]), rep(ncol(progs[[1]]),length(names(progs)))))
#factor(labs, levels = unique(labs))
#labs <- factor(labs, levels = c("nm","bfgs","gp.ei.opt","gp.ei.lhs","gp.ei.voriRIS","gp.ei.voriRLS","gp.ei.voriRAS"))
labs <- factor(labs, levels = short_method_names)

mids <- do.call(c, lapply(progs, function(p) p[half,]))
if (func=='rosen10') {
  mids <- log10(mids)
}
yl <- range(mids[!(labs %in% c('NM','BFGS'))], na.rm = T)
boxplot(mids~labs, main = paste("Iteration",half), xlab = '', las = 2, ylab='', ylim = yl)
#bp <- boxplot(mids~labs, main = paste("Iteration",half), xlab = '', las = 2, xaxt = "n")
#tick <- seq_along(bp$names)
#axis(1, at = tick, labels = FALSE)
#text(tick-0.3, par("usr")[3] - 0.5, bp$names, srt = 45, xpd = TRUE)

ends <- do.call(c, lapply(progs, function(p) p[end,]))
if (func=='rosen10') {
  ends <- log10(ends)
}
yl <- range(ends[!(labs %in% c('NM','BFGS'))], na.rm = T)
boxplot(ends~labs, main = paste("Iteration",end), xlab = '', las = 2, angle = 45, ylab='', ylim = yl)
dev.off()

### Time
#if (!addturbo) {
timecomps <- competitors[!(competitors%in% c("nm","bfgs"))]
files <- list.files(time_path, full.names = TRUE)

dfs <- list()
for (f in files) {
  dfs[[length(dfs)+1]] <-read.csv(f)
}

names <- dfs[[1]][,1]
tdf <- data.frame(sapply(dfs,function(x) x[,'times']))
rownames(tdf) <- names
tdf <- log10(tdf[timecomps,])

rownames(tdf) <- sapply(rownames(tdf), function(n) short_method_names[[n]])

pdf(paste(func,'_time.pdf',sep=''), width = 4, height = 4)
par(mar=c(4.7,3,1,0)+0.1)
par(mgp=c(2,0.6,0)+0.1)
#boxplot(t(tdf))
boxplot(t(tdf), main = 'Elapsed Real Time', xlab = '', las = 2, ylab='log10 Seconds')
dev.off()
#}


### Acquisition
#if (!addturbo) {
files <- list.files(crits_path, full.names = TRUE)

dfs <- list()
for (f in files) {
  dfs[[length(dfs)+1]] <-read.csv(f)
}

reps <- length(dfs)
progs <- lapply(competitors, function(i) matrix(NA,end-ninit,reps))
names(progs) <- competitors
torem <- c()
for (comp in competitors) {
  if ((comp %in% colnames(dfs[[1]]))) {
    for (i in 1:reps) {
      progs[[comp]][,i] <- dfs[[i]][,comp]
    }
  } else {
    torem <- c(torem, comp)
  }
}
pdf(paste(func,'_crits.pdf',sep=''), width = 4, height = 4)
par(mar=c(3.0,3,1,0.1)+0.1)
par(mgp=c(1.5,0.6,0)+0.1)

qu <- pmax(0,pmin(1,0.5+1.96*sqrt(0.25/reps)))
ql <- pmax(0,pmin(1,0.5-1.96*sqrt(0.25/reps)))

budget <- nrow(progs[[1]])
pcomp <- grep("^gp.*", competitors, value = TRUE)
ul <- 1.1*max(sapply(pcomp, function(x) quantile(progs[[x]][1,], qu, na.rm = T)), na.rm = T)
ll <- 0.9*min(sapply(pcomp, function(x) quantile(progs[[x]][budget,], ql, na.rm = T)), na.rm = T)

lwd <- 5
ylab <- 'Acquired Expected Improvement'
plot(NA,NA,xlim=c(1,end), ylim=c(ll,ul), main = pretty_sim_names[[func]], xlab='n: blackbox evaluations', ylab=ylab, font.lab = 2)


max_intervals <- Inf

for (comp in competitors) {
  print(comp)
  print(cols[[comp]])
  if (any(is.na(progs[[comp]]))) {
    print("Skipping nan competitor:")
    print(comp)
  } else {
    points((ninit+1):end, apply(progs[[comp]], 1, function(x) median(x)), type = 'l', lty = ltys[comp], col = cols[comp], lwd = 0.5*lwd)
    if (length(competitors) <= max_intervals) {
      points((ninit+1):end, apply(progs[[comp]], 1, function(x) quantile(x, ql)), type = 'l', lty = ltys[comp], col = cols[comp], lwd = 0.5*lwd)
      points((ninit+1):end, apply(progs[[comp]], 1, function(x) quantile(x, qu)), type = 'l', lty = ltys[comp], col = cols[comp], lwd = 0.5*lwd)
    }
  }
}
legend('topleft', legend = short_comp, lty=ltys,col=cols, lwd = lwd, bg = 'white', cex = 0.5)
dev.off()
