start=`date +%s`
source .venv/bin/activate
func=$1
Rscript R/clean_dir_make_args.R $func
nice -n 10 parallel --jobs 15 --colsep ' ' --will-cite -a sim_args.txt Rscript R/generic_sim.R
nice -n 10 parallel --jobs 15 --colsep ' ' --will-cite -a sim_args.txt python python/turbo_run.py
Rscript R/plot_generic.R $func
end=`date +%s`
runtime=$((end-start))
echo $runtime
