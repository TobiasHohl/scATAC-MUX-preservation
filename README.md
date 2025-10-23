# Enhancing single-cell ATAC sequencing with formaldehyde fixation, cryopreservation, and multiplexing for flexible analysis

Tobias Hohl<sup>1,2</sup>, Ulrike Boenisch<sup>1</sup>, Thomas Manke<sup>1,3</sup>, Laura Arrigoni<sup>1</sup>

- <sup>1</sup> Max-Planck-Institute of Immunobiology and Epigenetics, Freiburg, Germany

- <sup>2</sup> Faculty of Biology, Albert-Ludwigs-Universität Freiburg, Freiburg, Germany
- <sup>3</sup> School of Computing, Communication and Business, HTW Berlin, Berlin, Germany

This paper has been published in [BMC Research Notes](https://link.springer.com/article/10.1186/s13104-025-07547-y).

If you use parts or the entirety of the analyses available in this repository, please cite the main publication:
> Hohl, T., Bönisch, U., Manke, T., Arrigoni, L. Enhancing single-cell ATAC sequencing with formaldehyde fixation, cryopreservation, and multiplexing for flexible analysis. BMC Res Notes 18, 437 (2025). https://doi.org/10.1186/s13104-025-07547-y

## Abstract

> ### Objective
> The need for freshly isolated cells in bulk or single cell ATAC-seq experiments creates considerable logistical barriers and increases susceptibility to batch effects. This makes it difficult to coordinate complex or longitudinal studies. Our goal was to develop a sample preservation strategy that overcomes these limitations, enabling consistent and high-quality chromatin accessibility profiling from archived samples.
> ### Results
> We established a workflow that incorporates mild formaldehyde fixation prior to cryopreservation, preserving both bulk and single-cell ATAC-seq data quality at levels comparable to fresh samples in HepG2 cells. This protocol reliably maintains key data quality metrics, including signal-to-noise ratio and fragment distributions. Furthermore, the method is fully compatible with transposase-based sample multiplexing using custom Tn5 barcodes. To address barcode hopping inherent to multiplexing, we introduced a computational demultiplexing strategy based on fragment ratios, which accurately assigns single cells to their sample of origin. Our approach streamlines experimental logistics and ensures reproducibility across diverse and temporally dispersed samples, broadening the scope for ATAC-seq–based studies, including those in clinical research settings where coordinated sample collection is challenging.



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

You'll need working conda and R (>=4.3) installations and you'll need to install the dependencies to run the code. This is marked in the respective `.md` files.


## Reproducing the results

Please refer to the `.md` files in the respective folders.


## License

TBD