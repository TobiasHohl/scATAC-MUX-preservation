#!/bin/bash

startdate=$(date +%Y-%m-%d\ %T)

task=pygenometracks

echo "Started ${task} at: $startdate"

wd=${datadir}
outdir=${wd}/outs/${task}

mkdir -p $outdir

prev_dir=$(pwd)
cd ${repo_dir}/bulk/scripts/pygenometracks

pyGenomeTracks --tracks 01_ANXA5.ini --region chr4:121,665,946-121,698,980 -o ${outdir}/01_ANXA5_fresh.png
pyGenomeTracks --tracks 02_TUBB2A.ini --region chr6:3,151,666-3,159,544 -o ${outdir}/02_TUBB2A_fresh.png

cd $prev_dir

echo "done"

enddate=$(date +%Y-%m-%d\ %T)

echo "Finished at: $enddate"
