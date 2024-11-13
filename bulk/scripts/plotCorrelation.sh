#!/bin/bash

startdate=$(date +%Y-%m-%d\ %T)

task=plotCorrelation

echo "Started ${task} at: $startdate"

folder=${datadir}
outdir=${folder}/outs/${task}

mkdir -p $outdir

mbs=`find ${folder}/outs/multiBigWigSummary -name "*.npz" | sort`
logdir=$(pwd)/logs

echo "Submit ${task} jobs..."
for file in $mbs
do
	name=`basename $file .npz | sed s/^multiBigWigSummary_//`

	if [ $compute = "slurm" ]; then
		SlurmEasy -t 16 -m 1G -l $logdir -n plotCorrelation.${name}_pearson -k -v "plotCorrelation -in $file -c pearson -min 0 -max 1 -p heatmap --plotNumbers -o ${outdir}/plotCorrelation_${name}_pearson.pdf --outFileCorMatrix ${outdir}/plotCorrelation_${name}_pearson.tab; touch ${statdir}/${name}_${task}.finished"
	else
		plotCorrelation -in $file -c pearson -min 0 -max 1 -p heatmap --plotNumbers -o ${outdir}/plotCorrelation_${name}_pearson.pdf --outFileCorMatrix ${outdir}/plotCorrelation_${name}_pearson.tab
		touch ${statdir}/${i}_${task}.finished
	fi
done

echo "Waiting for ${task} jobs to finish..."

for file in $mbs
do
    fn=`basename $file .npz | sed s/^multiBigWigSummary_//`

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
