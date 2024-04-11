\name{getPCRanks}
\alias{getPCRanks}
\title{Get CpG eigenvector scores from a principle component.}
\usage{
getPCRanks(mat, IDs, PC, filter_thresh)
}
\description{
Returns eigenvector scores for input CpG sites.}
\arguments{
\item{mat}{Bismark2Matrix.R output file, or data frame object}
\item{IDs}{A character vector of IDs containing the common names for compared conditions. E.g., for samples trt1 & trt2 vs. ctl1 & ctl2, IDs=c("trt1", "ctl")}
\item{PC}{Integer, which principle component to use. Use to DefineBestPC if unsure.}
\item{filter_thresh}{Integer, a coverage threshold for filtering, where CpG coverage of all samples must be larger than this value.}
}
