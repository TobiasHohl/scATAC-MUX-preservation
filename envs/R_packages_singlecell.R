## Install required packages


if (!require("tidyverse")) install.packages("tidyverse")
if (!require("Signac")) {
  setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
  install.packages("Signac")
}
if (!require("BiocManager")) install.packages("BiocManager")
if (!require("GenomeInfoDb")) BiocManager::install("GenomeInfoDb")
if (!require("EnsDb.Hsapiens.v86")) BiocManager::install("EnsDb.Hsapiens.v86")
if (!require("GenomicRanges")) BiocManager::install("GenomicRanges")
if (!require("patchwork")) install.packages("patchwork")
if (!require("png")) install.packages("png")
if (!require("DT")) install.packages("DT")
if (!require("gridExtra")) install.packages("gridExtra")
if (!require("cowplot")) install.packages("cowplot")
if (!require("scales")) install.packages("scales")
