# PCBS - Principal Component BiSulfite
See our [package site](https://katlande.github.io/PCBS/index.html) 

## Dependencies
PCBS is an R package with the following dependencies:
* tibble
* ggrepel
* ggplot2
* dplyr
* data.table

## Data Setup
PCBS is a tool for analyzing WGBS datasets in a fast, flexible, and accurate fashion. PCBS is designed to pipe in Bismark-aligned WGBS data. The PCBS input file can be generated from bismark .cov files with the provided Bismark2Matrix scripts (we offer scripts for R and Python with identical usage). *We recommend using the Python version of this script as it is faster, but for those more comfortable with R the final output is identical.*

Python Usage: `python Bismark2Matrix.py file_path file_tsv file_out`

R Usage: `Rscript --vanilla  Bismark2Matrix.R -p file_path -i file_tsv -o file_out`

For both Python and R, the inputs are identical:
* file_path = /path/to/cov/files
* file_tsv = A three column, tab-separated file in this format:
  * filename | sample | group
  * sam1.cov | s1_trt | trt
  * sam2.cov | s2_ctl | ctl

* file_out = output file name

#### Important note: PCBS cannot handle CpGs that are NA an any sample, and these are removed by Bismark2Matrix. For this reason, PCBS is not recommended for sparse datasets. 


## Installation

#### Install the published version
install.packages("PCBS")

#### Install the latest updates
remotes::install_github("katlande/PCBS")

## Usage
See our [vignette](https://katlande.github.io/PCBS/articles/Differential_Methylation.html)


