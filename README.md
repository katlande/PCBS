# PCBS
Principle Component BiSulfite

## Dependencies
PCBS is an R package with the following dependencies:
* tibble
* ggrepel
* ggplot2

## Data Setup
PCBS is a tool for analyzing WGBS datasets in a fast, flexible, and accurate fashion. PCBS is designed to pipe in Bismark-aligned WGBS data. The PCBS input file can be generated from bismark .cov files with the provided Bismark2Matrix.R script.

Usage: Rscript --vanilla  Bismark2Matrix.R file_path file_tsv file_out
* file_path = /path/to/cov/files
* file_tsv = A three column, tab-separated file:
# file  name  group
# s1.cov s1_trt trt
# s2.cov  s2_ctl  ctl
# etc...

* file_out = output file name

## Usage
See our vignette
