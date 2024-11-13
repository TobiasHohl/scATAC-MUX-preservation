#!/bin/bash

startdate=$(date +%Y-%m-%d\ %T)

task=plotHeatmap

echo "Started ${task} at: $startdate"

folder=${datadir}
outdir=${folder}/outs/${task}

mkdir -p ${outdir}

echo "Execute ${task}..."

matrix_files=`find ${folder}/outs/computeMatrix -maxdepth 1 -name "*.gz" | sort`

labels='ENCODE 0.1_cryo 0.5_cryo 1.0_cryo 5.0_cryo 0.1_flash 0.5_flash 1.0_flash 5.0_flash'

echo "Submit plotHeatmap jobs..."
for file in $matrix_files
do
	name=`basename $file .gz | sed s/^computeMatrix_//`

	if [ $compute = "slurm" ]; then
		SlurmEasy -t 2 -m 75G -l $logdir -n plotHeatmap_$name -k -v "plotHeatmap -m $file -o /scratch/local/plotHeatmap_${name}.pdf --dpi 600 --refPointLabel 'center' --xAxisLabel '' --colorMap Reds Blues Greens Greens Greens Greens Purples Purples Purples Purples --whatToShow 'heatmap and colorbar' --heatmapWidth 3 --heatmapHeight 25 --samplesLabel $labels ; cp /scratch/local/plotHeatmap_${name}.pdf ${outdir}/plotHeatmap_${name}.pdf; touch ${statdir}/${name}_${task}.finished"
	else
		plotHeatmap -m $file -o ${outdir}/plotHeatmap_${name}.pdf --dpi 600 --refPointLabel 'center' --xAxisLabel '' --colorMap Reds Blues Greens Greens Greens Greens Purples Purples Purples Purples --whatToShow 'heatmap and colorbar' --heatmapWidth 3 --heatmapHeight 25 --samplesLabel $labels
		touch ${statdir}/${name}_${task}.finished
	fi
done

echo "Waiting for ${task} jobs to finish..."

for file in $matrix_files
do
    fn=`basename $file .gz | sed s/^computeMatrix_//`

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


echo "done"

enddate=$(date +%Y-%m-%d\ %T)

echo "Finished at: $enddate"

