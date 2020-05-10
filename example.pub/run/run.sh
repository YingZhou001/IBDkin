map="plink.map"
ind="ind.txt"
ibdfile="ibdhead.txt"
range="range.txt"

./IBDkin --ibdfile ${ibdfile} --map ${map} --ind ${ind} --range ${range} --nthreads 5 --out test --outmask --outcoverage
