ind <- as.numeric(ind)

## Sim settings.
#reps <- 2
reps <- 100
problems <- c("ackley10","levy10","rosen10","lunar","push","rover","pomp10log","dacca")

seed <- (ind-1) %% reps +1
pind <- floor((ind-1)/reps)+1
func <- problems[pind]

