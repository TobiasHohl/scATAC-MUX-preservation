#!/bin/bash

startdate=$(date +%Y-%m-%d\ %T)

task=bamcoverage

echo "Started ${task} at: $startdate"

wd=${datadir}
genomesize=2805636331
outdir=${wd}/bamcoverage
statdir=${wd}/status

mkdir -p $outdir


bamdir=${wd}/sp/filtered_bam

files=`find ${bamdir} -maxdepth 1 -name "*.bam" | sort`

for file in $files
do
	i=`basename $file .bam`
	if test -f ${outdir}/${i}.bw
	then
		echo -e "File ${i}.bw already exists."
	else
		if [ $compute = "slurm" ]; then
			SlurmEasy -t 16 -m 1G -l $logdir -n bamCoverage_${i} -k -v "bamCoverage -b $file -o ${outdir}/${i}.bw -bl $blacklist --normalizeUsing RPKM -p 16 -bs 1; touch ${statdir}/${i}_${task}.finished"
		else
			bamCoverage -b $file -o ${outdir}/${i}.bw -bl $blacklist --normalizeUsing RPKM -p 1 -bs 1
			touch ${statdir}/${i}_${task}.finished
		fi
		
	fi
done

echo "Waiting for ${task} jobs to finish..."

for file in $files
do
    fn=`basename $file .bam`

    if [ ! -f ${statdir}/${fn}_${task}.finished ]
    then
        echo "Waiting for ${fn}_${task} to finish..."
        while [ ! -f ${statdir}/${fn}_${task}.finished ]
        do
            sleep 1m
        done
    fi
	rm -rf ${statdir}/*
done

echo "done"

enddate=$(date +%Y-%m-%d\ %T)

echo "Finished ${task} at: $enddate"
