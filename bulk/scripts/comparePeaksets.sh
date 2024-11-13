#!/bin/bash

startdate=$(date +%Y-%m-%d\ %T)

task=comparePeaksets

echo "Started ${task} at: $startdate"

folder=${datadir}
outdir=${wd}/outs/${task}

mkdir -p ${outdir}

echo "Execute ${task}..."

cryo_peaks=${folder}/sp/MACS2/cryo_0.1FA.filtered.short.BAM_peaks.narrowPeak
ref_peaks=$bed_ref

# Go to $outdir

cd $outdir

# Prepare data and limit columns to chr | start | end 

cut -f 1-3 $cryo_peaks > cryo_peaks.bed
cut -f 1-3 $ref_peaks > ref_peaks.bed

# Combine

cat cryo_peaks.bed > all.bed
cat ref_peaks.bed >> all.bed

# Merge

sort -k1,1 -k2,2n all.bed > all.sorted.bed
bedtools merge -i all.sorted.bed > all.merged.bed

# Get merged peaks that fit either ref or cryo peaks

bedtools intersect -u -a all.merged.bed -b cryo_peaks.bed > cryo_all.bed
bedtools intersect -u -a all.merged.bed -b ref_peaks.bed > ENCODE_all.bed

# Compare and extract shared and unique peak sets

bedtools intersect -u -a cryo_all.bed -b ENCODE_all.bed > shared.bed
bedtools intersect -v -a cryo_all.bed -b ENCODE_all.bed > unique_cryo.bed
bedtools intersect -v -a ENCODE_all.bed -b cryo_all.bed > unique_ENCODE.bed

echo "done"

enddate=$(date +%Y-%m-%d\ %T)

echo "Finished ${task} at: $enddate"