#start=`date +%s`
#func=$1
#echo $func
Rscript R/clean_dir_make_args.R $func
#if [ "$func" == "rover" ]; then
#    echo "using fewer jobs for rover"
#    #nice -n 10 parallel --jobs 15 --colsep ' ' --will-cite -a sim_args.txt Rscript R/generic_sim.R
#    nice -n 10 parallel --jobs 10 --colsep ' ' --will-cite -a sim_args.txt Rscript R/generic_sim.R
#else 
#    nice -n 10 parallel --jobs 20 --colsep ' ' --will-cite -a sim_args.txt Rscript R/generic_sim.R
#fi
#Rscript R/plot_generic.R $func
#end=`date +%s`
#runtime=$((end-start))
#echo $runtime
