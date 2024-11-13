#!/bin/bash

startdate=$(date +%Y-%m-%d\ %T)

task=computeMatrix

echo "Started ${task} at: $startdate"

folder=${datadir}
outdir=${folder}/outs/${task}

mkdir -p ${outdir}

echo "Execute ${task}..."

type=Fig2

bigwigs=`echo "${folder}/bamcoverage/ENCODE.filtered.bw ${folder}/bamcoverage/cryo_0.1FA.filtered.bw ${folder}/bamcoverage/cryo_0.5FA.filtered.bw ${folder}/bamcoverage/cryo_1.0FA.filtered.bw ${folder}/bamcoverage/cryo_5.0FA.filtered.bw ${folder}/bamcoverage/flash_0.1FA.filtered.bw ${folder}/bamcoverage/flash_0.5FA.filtered.bw ${folder}/bamcoverage/flash_1.0FA.filtered.bw ${folder}/bamcoverage/flash_5.0FA.filtered.bw"`
labels=`echo "ENCODE cryo_0.1 cryo_0.5 cryo_1.0 cryo_5.0 flash_0.1 flash_0.5 flash_1.0 flash_5.0"`

if [ $compute = "slurm" ]; then
	SlurmEasy -t 40 -m 2G -l $logdir -n computeMatrix_${type}_refPeaks -k -v "computeMatrix reference-point --referencePoint center -b 3000 -a 3000 -bs 1 --skipZeros -p 40 --missingDataAsZero -S $bigwigs -R $bed_ref --samplesLabel $labels -o ${outdir}/computeMatrix_${type}_refPeaks.gz; touch ${statdir}/${task}_${type}.finished"
else
	computeMatrix reference-point --referencePoint center -b 3000 -a 3000 -bs 1 --skipZeros -p 1 --missingDataAsZero -S $bigwigs -R $bed_ref --samplesLabel $labels -o ${outdir}/computeMatrix_${type}_refPeaks.gz
    touch ${statdir}/${task}_${type}.finished
fi

echo "Waiting for ${task} jobs to finish..."

if [ ! -f ${statdir}/${task}_${type}.finished ]
then
    echo "Waiting for ${fn}_${task} to finish..."
    while [ ! -f ${statdir}/${task}_${type}.finished ]
    do
        sleep 1m
    done
fi
rm -rf ${statdir}/*

echo "done"

enddate=$(date +%Y-%m-%d\ %T)

echo "Finished ${task} at: $enddate"


type=SupplFig2

beds=`echo "${folder}/outs/comparePeaksets/unique_cryo.bed ${folder}/outs/comparePeaksets/unique_ENCODE.bed ${folder}/outs/comparePeaksets/shared.bed"`

if [ $compute = "slurm" ]; then
	SlurmEasy -t 40 -m 2G -l $logdir -n computeMatrix_${type}_refPeaks -k -v "computeMatrix reference-point --referencePoint center -b 3000 -a 3000 -bs 1 --skipZeros -p 40 --missingDataAsZero -S $bigwigs -R $beds --samplesLabel $labels -o ${outdir}/computeMatrix_${type}.gz"
else
	computeMatrix reference-point --referencePoint center -b 3000 -a 3000 -bs 1 --skipZeros -p 40 --missingDataAsZero -S $bigwigs -R $beds --samplesLabel $labels -o ${outdir}/computeMatrix_${type}.gz
    touch ${statdir}/${task}_${type}.finished
fi

echo "Waiting for ${task} jobs to finish..."

if [ ! -f ${statdir}/$${task}_${type}.finished ]
then
    echo "Waiting for ${fn}_${task} to finish..."
    while [ ! -f ${statdir}/${task}_${type}.finished ]
    do
        sleep 1m
    done
fi
rm -rf ${statdir}/*

echo "done"

enddate=$(date +%Y-%m-%d\ %T)

echo "Finished ${task} at: $enddate"
