newhosts='p105|p106|p107|p108|p109|p110|p111|p112|p113|p114|p115|p116'
nthreads=12
ibdcm="2"
seed=0.5
ext=0.1
gap=1000

hapibd="java -jar -Xmx120g /projects/browning/software/hap-ibd.jar"

for chr in 1 2 3 4
do
map="/projects/browning/ukbio/joe/data/maps/plink.chr${chr}.const.map"
vcf=../vcf/chr${chr}.vcf.gz
out=chr${chr}

echo "time ${hapibd} gt=${vcf} out=${out} map=${map} min-output=${ibdcm} min-seed=${seed} min-extend=${ext} nthreads=${nthreads}" |  qsub -N "o.hi.${out}" -q b-all.q -l h=${newhosts} -m n -pe local ${nthreads}

done
