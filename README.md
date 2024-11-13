# Enhancing single-cell ATAC sequencing with formaldehyde fixation, cryopreservation, and multiplexing for flexible analysis

Tobias Hohl<sup>1,2</sup>, Ulrike Boenisch<sup>1</sup>, Thomas Manke<sup>1</sup>, Laura Arrigoni<sup>1</sup>

- <sup>1</sup> Max-Planck-Institute of Immunobiology and Epigenetics, Freiburg

- <sup>2</sup> Faculty of Biology, Freiburg University, Freiburg

This paper has been submitted for publication in *Some Journal*.

> The paper showcases that integrating formaldehyde fixation with cryopreservation in microfluidics-based single-cell ATAC-seq experiments produces results comparable to those from fresh samples. We provide a multiplexed scATAC-seq workflow that allows for flexible sample collection, enabling time-course experiments to be processed within a single 10x scATAC-seq reaction. This method eliminates batch effects arising from different sampling and library preparation experiments due to streamlined preservation and parallel processing.


## Abstract

> The assay for transposase-accessible chromatin using sequencing (ATAC-seq) revolutionized the field of epigenetics since its emergence by providing a means to uncover chromatin dynamics and other factors affecting gene expression. The development of single-cell (sc) applications in recent years led to an even deeper understanding of cell type specific gene regulatory mechanisms. One of the major challenges while running ATAC-seq experiments, bulk or sc, is the need for freshly collected cells for successful experiments. While various freezing methods have already been tested and established for bulk and sc ATAC-seq, quality metrics for preserved cells are rather poor or dependent on sampling time when compared to fresh samples. This makes it difficult to conduct all sorts of complex experiments i.e. with multiple conditions, patients, or time course studies. Especially, accounting for batch effects can be difficult if samples need to be processed at different time points of collection. We tackled this issue by adding a fixation step prior to the freezing method. The additional fixation step improved library quality and yield data comparable to fresh samples. The workflow was also tested on multiplexed sc ATAC experiments, set-up for cost-efficient low input sample handling. Sample cross-in, typically encountered in Tn5-based multiplex approaches, were tackled with a computational procedure specifically developed for this approach.



## Software implementation

The code for generating custom Tn5 barcodes can be found in the folder `custombarcodedesign`.

The source code used to generate the results and figures in the paper are in
the `bulk` and `singlecell` folders, respectively.
The calculations and figure generation are using a variety of tools, with R being the main language used for visualization.
The data used in this study is uploaded to GEO and can be retrieved using the accession numbers mentioned in the paper.
See the `.md` files in each directory for a full description on how to reproduce the results.

## Getting the code

You can download a copy of all the files in this repository by cloning the
[git](https://git-scm.com/) repository:

    git clone https://github.com/TobiasHohl/scATAC-MUX-preservation.git

or [download a zip archive](https://github.com/TobiasHohl/scATAC-MUX-preservation/archive/master.zip).


## Dependencies

All dependencies (conda environments + R packages needed) are listed in the `envs` folder.

You'll need working conda and R (>=4.3) installation and you'll need to install the dependencies to run the code. This is marked in the respective `.md` files.


## Reproducing the results

Please refer to the `.md` files in the respective folders.


## License

TBD