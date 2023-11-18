func <- commandArgs(trailingOnly=TRUE)[1]

source("R/sim_settings.R")

## Clear dir
#dir.create(file.path(root, func), showWarnings = FALSE)
dir.create(sim_path, showWarnings = FALSE)
do.call(file.remove, list(list.files(sim_path, full.names = TRUE)))

dir.create(time_path, showWarnings = FALSE)
do.call(file.remove, list(list.files(time_path, full.names = TRUE)))

dir.create(crits_path, showWarnings = FALSE)
do.call(file.remove, list(list.files(crits_path, full.names = TRUE)))

## Make Args.
aa <- list()
for (r in 1:reps) {
    aa[r] <- paste(func,r)
}

write(paste(aa,collapse='\n'), 'sim_args.txt')