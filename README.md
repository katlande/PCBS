# PCBS - Principle Component BiSulfite

## Dependencies
PCBS is an R package with the following dependencies:
* tibble
* ggrepel
* ggplot2
* dplyr

## Data Setup
PCBS is a tool for analyzing WGBS datasets in a fast, flexible, and accurate fashion. PCBS is designed to pipe in Bismark-aligned WGBS data. The PCBS input file can be generated from bismark .cov files with the provided Bismark2Matrix.R script.

Usage: Rscript --vanilla  Bismark2Matrix.R file_path file_tsv file_out
* file_path = /path/to/cov/files
* file_tsv = A three column, tab-separated file in this format:
  * filename | sample | group
  * sam1.cov | s1_trt | trt
  * sam2.cov | s2_ctl | ctl

* file_out = output file name


## Installation
devtools::install_github("katlande/PCBS")


## Usage
See our [vignette](https://github.com/katlande/PCBS/blob/main/PCBS_Vignette.md)
