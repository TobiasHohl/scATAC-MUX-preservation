# Analysis of the multiplexed single-cell ATAC-seq data

## Data preprocessing using cellranger-atac
- [Download](https://software.10xgenomics.com/single-cell-atac/software/downloads/latest?) and [install](https://software.10xgenomics.com/single-cell-atac/software/pipelines/latest/installation) cellranger-atac v2.1.0.
- [Download](https://software.10xgenomics.com/single-cell-atac/software/downloads/latest) the appropriate reference (GRCh38 Reference - 2020-A-2.0.0 (May 3, 2021))
- Download the fastq files from GEO using the accession numbers found in the paper.
- Prepare a folder with the fastq files
- Run cellranger-atac on each sample individually using the following command:
 - can be done via script and also submitted to an HPC
 - set the same output dir for all samples
```
cd /path/to/cellranger/output/dir/
cellranger-atac count --id=sample_name \
                --fastqs=/path/to/fq/folder/ \
                --sample=sample_name \
                --reference=/path/to/reference/ \
                --localcores 80 \
                --localmem 240
```

## Data analysis using Signac
- Install all necessary packages
 - Use `R_packages_singlecell.R` script in `./envs` directory
- Follow the walkthrough in the `singlecell_workflow.Rmd` markdown file to analyze the data and recreate the figures.