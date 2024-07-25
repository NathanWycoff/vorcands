ind <- as.numeric(ind)

## Sim settings.
reps <- 2
#reps <- 100
problems <- c("ackley10","levy10","rosen10","lunar","push","rover","pomp10","dacca")

seed <- ind %% reps
pind <- ceiling(ind/reps)
func <- problems[pind]

