version="v2.8.7.3"

ibdkin="/projects/browning/ukbio/joe/ibd.project.script.dev/IBDkin/src/${version}-pub/IBDkin"
cp ${ibdkin} ./

map="plink.map"
ind="ind.txt"
ibdfile="ibdhead.txt"

./IBDkin -ibdfile ${ibdfile} -map ${map} -ind ${ind} -nthreads 5 -out test --outmask --outcoverage
