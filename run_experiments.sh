start_big=`date +%s`
#for func in ackley10 levy10 rosen10 push rover pomp10log dacca lunar 
for func in ackley10 levy10 lunar rosen10 push rover pomp10log dacca
do
    echo $func
    ./parrun.sh $func
done
end_big=`date +%s`
runtime_big=$((end_big-start_big))
echo $runtime_big

