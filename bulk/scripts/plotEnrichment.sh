#!/bin/bash

startdate=$(date +%Y-%m-%d\ %T)

task=plotEnrichment

echo "Started ${task} at: $startdate"

wd=${datadir}
outdir=${wd}/outs/${task}

mkdir -p $outdir

echo "Execute ${task}..."

bams=`echo "${wd}/sp/filtered_bam/ENCODE.filtered.bam ${wd}/sp/filtered_bam/cryo_0.1FA.filtered.bw ${wd}/sp/filtered_bam/cryo_0.5FA.filtered.bw ${wd}/sp/filtered_bam/cryo_1.0FA.filtered.bw ${wd}/sp/filtered_bam/cryo_5.0FA.filtered.bw ${wd}/sp/filtered_bam/flash_0.1FA.filtered.bw ${wd}/sp/filtered_bam/flash_0.5FA.filtered.bw ${wd}/sp/filtered_bam/flash_1.0FA.filtered.bw ${wd}/sp/filtered_bam/flash_5.0FA.filtered.bw"`

bams_labels=`echo "ref MUX-cryo_0.1FA MUX-cryo_0.5FA MUX-cryo_1.0FA MUX-cryo_5.0FA MUX-ff_0.1FA MUX-ff_0.5FA MUX-ff_1.0FA MUX-ff_5.0FA Std-ff_0.1FA Std-ff_0.5FA Std-ff_1.0FA Std-ff_5.0FA"`


beds=$bed_ref
bed_labels=`echo "ENCODE_peaks"`

if [ $compute = "slurm" ]; then
	SlurmEasy -t 10 -m 1G -l $logdir -n plotEnrichment -k -v "plotEnrichment --bamfiles $bams --BED $beds -p 1 -l $bams_labels --plotWidth 20 --plotHeight 15 --regionLabels $bed_labels -o ${outdir}/plot.png --outRawCounts ${outdir}/values.tsv; touch ${statdir}/${task}.finished"
else
	plotEnrichment --bamfiles $bams --BED $beds -p 1 -l $bams_labels --plotWidth 20 --plotHeight 15 --regionLabels $bed_labels -o ${outdir}/plot.png --outRawCounts ${outdir}/values.tsv
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
