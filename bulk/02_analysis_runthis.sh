#!/bin/bash



# SET THESE PARAMETERS
datadir=/path/to/datadir
repo_dir=/path/to/git/repo


# SET TO 'slurm' IF HPC WITH SLURM AVAILABLE
compute=local

# ONLY SET IF YOU USE DIFFERENT FILES

blacklist=${repo_dir}/misc/GRCh38_ensembl_plusMito.bed
bed_ref=${repo_dir}/misc/ENCODE.bed

# Set accordingly

cryo_peaks=${folder}/sp/MACS2/cryo_0.1FA.filtered.short.BAM_peaks.narrowPeak


# DO NOT TOUCH
logdir=$(pwd)/logs
statdir=${datadir}/status

mkdir -p $statdir
mkdir -p ${datadir}/outs/Rplots


# Run tasks

. ${repo_dir}/bulk/scripts/bamcoverage.sh
. ${repo_dir}/bulk/scripts/idxstat.sh
. ${repo_dir}/bulk/scripts/bampefragmentsize.sh
. ${repo_dir}/bulk/scripts/multiBigWigSummary.sh
. ${repo_dir}/bulk/scripts/plotCorrelation.sh
. ${repo_dir}/bulk/scripts/plotEnrichment.sh
. ${repo_dir}/bulk/scripts/plotPCA.sh
. ${repo_dir}/bulk/scripts/comparePeaksets.sh
. ${repo_dir}/bulk/scripts/computeMatrix.sh
. ${repo_dir}/bulk/scripts/plotHeatmap.sh

# . ${repo_dir}/bulk/scripts/pygenometracks/pygenometracks.sh


# NECESSARY FOR RUNNING TOBIAS:

# CHANGE THIS
genome_fa=/path/to/GRCh38_ensembl/genome_fasta/genome.fa

jaspar_motifs=${repo_dir}/misc/jaspar_all-human.txt

reference=ENCODE
peaks=$bed_ref


# SLURM is required for running TOBIAS. Modify accordingly.
if [ $compute = "slurm" ]; then

	module load slurm

    read -e -p "Maximum parallel jobs (0 for infinite): " max_j
    if [[ $max_j -gt 0 ]]
    then
        read -e -p "Enter cluster user: " user
    fi

    . ${repo_dir}/bulk/scripts/TOBIAS.sh

else
	echo "SLURM is required to run TOBIAS. Please set the 'compute' parameter to 'slurm' and adjust lines 53 - 59 of script 02_analysis_runthis.sh accordingly."
fi

