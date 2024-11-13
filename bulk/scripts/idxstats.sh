#!/bin/bash

startdate=$(date +%Y-%m-%d\ %T)

task=idxstats

echo "Started ${task} at: $startdate"

wd=${datadir}
outdir=${wd}/outs/${task}


mkdir -p $outdir

bamdir=${wd}/sp/filtered_bam

files=`find ${bamdir} -name "*.bam" | sort`

for file in $files
do
	i=`basename $file .bam`
	if test -f ${outdir}/${i}.txt
	then
		echo -e "File ${i}.txt already exists."
	else
		if [ $compute = "slurm" ]; then
			SlurmEasy -t 2 -m 1G -l $logdir -n idxstats_${i} -k -v "samtools idxstats $file > ${outdir}/${i}.txt; touch ${statdir}/${i}_${task}.finished"
		else
			samtools idxstats $file > ${outdir}/${i}.txt
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
