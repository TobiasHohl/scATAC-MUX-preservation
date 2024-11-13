#!/bin/bash

startdate=$(date +%Y-%m-%d\ %T)

task=plotPCA

echo "Started ${task} at: $startdate"

folder=${datadir}
outdir=${folder}/outs/${task}

mkdir -p ${outdir}

echo "Execute ${task}..."

mbs=`find ${folder}/multiBigWigSummary -name "*.npz" | sort`

echo "Submit ${task} jobs..."
for file in $mbs
do
	name=`basename $file .npz | sed s/^multiBigWigSummary_//`
	if [ $compute = "slurm" ]; then
		SlurmEasy -t 16 -m 1G -l $logdir -n plotPCA.${name}.rC -k -v "plotPCA -in $file -o ${outdir}/plotPCA_${name}_rowCenter.pdf --ntop 0 --rowCenter --outFileNameData ${outdir}/values_${name}.tsv; touch ${statdir}/${name}_${task}.finished"
	else
		plotPCA -in $file -o ${outdir}/plotPCA_${name}_rowCenter.pdf --ntop 0 --rowCenter --outFileNameData ${outdir}/values_${name}.tsv
		touch ${statdir}/${name}_${task}.finished
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
