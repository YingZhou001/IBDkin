map="plink.map"
ind="ind.txt"
ibdfile="ibd.txt"

#example-1
./IBDkin -ibdfile ${ibdfile} -map ${map} -ind ${ind} -nthreads 5 -out example-1


#example-2
./IBDkin -ibdfile ${ibdfile} -map ${map} -ind ${ind} -nthreads 5 -out example-2 --outmask --outcoverage --nokinship

#example-3
for part in {1..5}
do
  ./IBDkin -ibdfile ${ibdfile} -map ${map} -ind ${ind} -nthreads 5 -out example-3.${part} -part 5 ${part}
done

