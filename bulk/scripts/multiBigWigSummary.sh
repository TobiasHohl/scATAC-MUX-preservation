#!/bin/bash

startdate=$(date +%Y-%m-%d\ %T)

task=multiBigWigSummary

echo "Started ${task} at: $startdate"

folder=${datadir}
outdir=${folder}/outs/${task}

mkdir -p ${outdir}

echo "Execute ${task}..."

bws=`echo "${folder}/bamcoverage/ENCODE.filtered.bw ${folder}/bamcoverage/cryo_0.1FA.filtered.bw ${folder}/bamcoverage/cryo_0.5FA.filtered.bw ${folder}/bamcoverage/cryo_1.0FA.filtered.bw ${folder}/bamcoverage/cryo_5.0FA.filtered.bw ${folder}/bamcoverage/flash_0.1FA.filtered.bw ${folder}/bamcoverage/flash_0.5FA.filtered.bw ${folder}/bamcoverage/flash_1.0FA.filtered.bw ${folder}/bamcoverage/flash_5.0FA.filtered.bw"`
labels=`echo "ENCODE cryo_0.1 cryo_0.5 cryo_1.0 cryo_5.0 flash_0.1 flash_0.5 flash_1.0 flash_5.0"`

if [ $compute = "slurm" ]; then
	SlurmEasy -t 40 -m 2G -l $logdir -n multiBigWigSummary -k -v "multiBigwigSummary BED-file -b $bws -l $labels -o ${outdir}/${task}_all.npz --BED $bed_ref -p 40; touch ${statdir}/${task}.finished"
else
	multiBigwigSummary BED-file -b $bws -l $labels -o ${outdir}/${task}_all.npz --BED $bed_ref -p 1
	touch ${statdir}/${task}.finished
fi

echo "Waiting for ${task} jobs to finish..."

if [ ! -f ${statdir}/${task}.finished ]
then
    echo "Waiting for ${fn}_${task} to finish..."
    while [ ! -f ${statdir}/${task}.finished ]
    do
        sleep 1m
    done
fi
rm -rf ${statdir}/*

echo "done"

enddate=$(date +%Y-%m-%d\ %T)

echo "Finished ${task} at: $enddate"
