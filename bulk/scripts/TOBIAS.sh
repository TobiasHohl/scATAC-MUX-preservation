#!/bin/bash

wd=${datadir}
bamdir=${wd}/sp/filtered_bams

bams=`find $bamdir -name '*bam'`

outdir=${wd}/outs/TOBIAS
statdir=${outdir}/status
logdir=${outdir}/logs

mkdir -p $statdir
mkdir -p $outdir

cd $outdir

task1=ATACcorrect

echo -e "\nTask 1: ${task1}\n"

mkdir -p ${outdir}/${task1}

for bam in $bams
do
    fn=`basename $bam | sed 's/.bam//'`
    if [ -f ${statdir}/${fn}_${task1}.finished ] ;
    then
        echo "Correction already done for ${fn}. Skipping..."
    else
        echo "Submit $fn to slurm."
        run_j=`squeue -u $user | wc -l`
		run_j=`expr $run_j - 1`
		if [[ $run_j -ge $max_j ]]
		then
			echo "Maximum job amount reached. Waiting for jobs to finish..."
		fi
		while [[ $run_j -ge $max_j ]]
		do
			sleep 60
			run_j=`squeue -u $user | wc -l`
        		run_j=`expr $run_j - 1`
		done
        SlurmEasy -t 8 -m 2G -l $logdir -n TOBIAS_${task1}_${fn} -k -v "TOBIAS ATACorrect --bam $bam --genome $genome_fa --peaks $peaks --blacklist $blacklist --outdir ${outdir}/${task1} --cores 8; touch ${statdir}/${fn}_${task1}.finished"
    fi
done

echo "Waiting for $task1 jobs to finish..."

for bam in $bams
do
    fn=`basename $bam | sed 's/.bam//'`

    if [ ! -f ${statdir}/${fn}_${task1}.finished ]
    then
        echo "Waiting for ${fn}_${task1} to finish..."
        while [ ! -f ${statdir}/${fn}_${task1}.finished ]
        do
            sleep 1m
        done
    fi
done

echo -e "All $task1 jobs finished.\n"

task2=FootprintScores

echo -e "\nTask 2: ${task2}\n"

mkdir -p ${outdir}/${task2}

for bam in $bams
do
    fn=`basename $bam | sed 's/.bam//'`
    if [ -f ${statdir}/${fn}_${task2}.finished ] ;
    then
        echo "$task2 already done for ${fn}. Skipping..."
    else
        echo "Submit ${fn}_${task2} to slurm."
        run_j=`squeue -u $user | wc -l`
		run_j=`expr $run_j - 1`
		if [[ $run_j -ge $max_j ]]
		then
			echo "Maximum job amount reached. Waiting for jobs to finish..."
		fi
		while [[ $run_j -ge $max_j ]]
		do
			sleep 60
			run_j=`squeue -u $user | wc -l`
        		run_j=`expr $run_j - 1`
		done
        SlurmEasy -t 8 -m 2G -l $logdir -n TOBIAS_${task2}_${fn} -k -v "TOBIAS FootprintScores --signal ${outdir}/${task1}/${fn}_corrected.bw --regions $peaks --output ${outdir}/${task2}/${fn}_footprints.bw --cores 8; touch ${statdir}/${fn}_${task2}.finished"
    fi
done

echo "Waiting for $task2 jobs to finish..."

for bam in $bams
do
    fn=`basename $bam | sed 's/.bam//'`

    if [ ! -f ${statdir}/${fn}_${task2}.finished ]
    then
        echo "Waiting for ${fn}_${task2} to finish..."
        while [ ! -f ${statdir}/${fn}_${task2}.finished ]
        do
            sleep 1m
        done
    fi
done

echo -e "All $task2 jobs finished.\n"

task3=BINDetect

echo -e "\nTask 3: ${task3}\n"

mkdir -p ${outdir}/${task3}

for bam in $bams
do
    fn=`basename $bam | sed 's/.bam//'`
    if [ -f ${statdir}/${fn}_${task3}.finished ] ;
    then
        echo "$task3 already done for ${fn}. Skipping..."
    else
        echo "Submit ${fn}_${task3} to slurm."
        run_j=`squeue -u $user | wc -l`
		run_j=`expr $run_j - 1`
		if [[ $run_j -ge $max_j ]]
		then
			echo "Maximum job amount reached. Waiting for jobs to finish..."
		fi
		while [[ $run_j -ge $max_j ]]
		do
			sleep 60
			run_j=`squeue -u $user | wc -l`
        		run_j=`expr $run_j - 1`
		done
        SlurmEasy -t 8 -m 2G -l $logdir -n TOBIAS_${task3}_${fn} -k -v "TOBIAS BINDetect --motifs $jaspar_motifs --signals ${outdir}/${task2}/${fn}_footprints.bw --genome $genome_fa --peaks $peaks --outdir ${outdir}/${task3}/${fn} --cond_names ${fn} --cores 8; touch ${statdir}/${fn}_${task3}.finished"
    fi
done

echo "Waiting for $task3 jobs to finish..."

for bam in $bams
do
    fn=`basename $bam | sed 's/.bam//'`

    if [ ! -f ${statdir}/${fn}_${task3}.finished ]
    then
        echo "Waiting for ${fn}_${task3} to finish..."
        while [ ! -f ${statdir}/${fn}_${task3}.finished ]
        do
            sleep 1m
        done
    fi
done

echo -e "All $task3 jobs finished.\n"


task4=PlotAggregate

echo -e "\nTask 4: ${task4}\n"

mkdir -p ${outdir}/${task4}

bws=`find ${outdir}/${task1}/ -name "*_corrected.bw" | sort | tr '\n' ' '`

motifs=`find ${outdir}/${task3}/$reference -name "*all.bed"`

for motif in $motifs
do
    fn=`echo $motif | sed 's/^.*beds\/MA[0-9]*.[0-9].//' | sed 's/_all.bed//'`
    if [ -f ${statdir}/${task4}_${fn}.finished ] ;
    then
        echo "$task4 already done for ${fn}. Skipping..."
    else
        echo "Submit ${fn}_${task4} to slurm."
        run_j=`squeue -u $user | wc -l`
		run_j=`expr $run_j - 1`
		if [[ $run_j -ge $max_j ]]
		then
			echo "Maximum job amount reached. Waiting for jobs to finish..."
		fi
		while [[ $run_j -ge $max_j ]]
		do
			sleep 60
			run_j=`squeue -u $user | wc -l`
        		run_j=`expr $run_j - 1`
		done
        SlurmEasy -t 8 -m 2G -l $logdir -n TOBIAS_${task4}_${fn} -k -v "TOBIAS PlotAggregate --TFBS $motif  --signals $bws --output ${outdir}/${task4}/${fn}.pdf --share_y both --plot_boundaries --signal-on-x --log-transform --smooth 5 --output-txt ${outdir}/${task4}/${fn}.txt; touch ${statdir}/${task4}_${fn}.finished"
    fi
done

echo "Waiting for $task4 jobs to finish..."

for motif in $motifs
do
    fn=`echo $motif | sed 's/^.*beds\/MA[0-9]*.[0-9].//' | sed 's/_all.bed//'`

    if [ ! -f ${statdir}/${task4}_${fn}.finished ]
    then
        echo "Waiting for ${fn}_${task4} to finish..."
        while [ ! -f ${statdir}/${task4}_${fn}.finished ]
        do
            sleep 1m
        done
    fi
done

echo -e "All $task4 jobs finished.\n"
