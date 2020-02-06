map="plink.map"
ind="ind.txt"
ibdfile="ibdhead.txt"

./IBDkin -ibdfile ${ibdfile} -map ${map} -ind ${ind} -nthreads 5 -out test --outmask --outcoverage
