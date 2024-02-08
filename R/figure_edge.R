load("Rdata/bound.RData")

adf <- aggregate(resdf, by = POB ~ P+N+norm+st, mean)

#sts <- c("unif","rect","lhs")


#styles <- 1:length(Ps)
cols <- list(l1='red',l2='blue',linf='green')

for (st in c('rect','unif','lhs')) {
    sdf <- adf[adf$st==st,]

    Ps <- unique(sdf$P)
    Ns <- unique(sdf$N)
    norms <- unique(sdf$norm)

    pdf(paste('figures/borders_',st,".pdf",sep=''), height = 2, width = 6)
    par(mfrow=c(1,3))
    par(mar=c(3.0,3,1,0.1)+0.1)
    par(mgp=c(1.5,0.6,0)+0.1)
    for (P in Ps) {
        tit <- paste(st, '; dim = ', P, sep = '')
        plot(NA,NA,xlim=range(log10(Ns)),ylim=c(0,1), main = tit, ylab = 'Proportion on Boundary', xlab = 'Log N')
        for (norm in norms) {
            col <- cols[[norm]]
            ssdf <- sdf[sdf$P==P & sdf$norm==norm,]
            lines(log10(ssdf$N), ssdf$POB, col = col)
        }
        if (P==Ps[1]) {
            legend('topleft', legend = norms, fill = sapply(cols, function(x) x))
        }
    }
    dev.off()
}