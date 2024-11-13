#!/bin/bash

startdate=$(date +%Y-%m-%d\ %T)

task=bamPEFragmentSize

echo "Started ${task} at: $startdate"

wd=${datadir}
outdir=${wd}/outs/${task}

mkdir -p $outdir

bams=`echo "${bamdir}/ENCODE.filtered.bam ${bamdir}/cryo_0.1FA.filtered.bam ${bamdir}/cryo_0.5FA.filtered.bam ${bamdir}/cryo_1.0FA.filtered.bam ${bamdir}/cryo_5.0FA.filtered.bam ${bamdir}/flash_0.1FA.filtered.bam ${bamdir}/flash_0.5FA.filtered.bam ${bamdir}/flash_1.0FA.filtered.bam ${bamdir}/flash_5.0FA.filtered.bam"`
labels=`echo "ENCODE 0.1_cryo 0.5_cryo 1.0_cryo 5.0_cryo 0.1_flash 0.5_flash 1.0_flash 5.0_flash"`

if [ $compute = "slurm" ]; then
    SlurmEasy -t 40 -m 2G -l $logdir -n bamPEFragmentSize -k -v "bamPEFragmentSize -b $bams -o ${outdir}/total.png -p 40 --samplesLabel $labels --outRawFragmentLengths ${outdir}/values.tsv; touch ${statdir}/${task}.finished"
else
    bamPEFragmentSize -b $bams -o ${outdir}/total.png -p 1 --samplesLabel $labels --outRawFragmentLengths ${outdir}/values.tsv
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
