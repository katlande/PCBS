\name{getPCRanks}
\alias{getPCRanks}
\title{Get CpG eigenvector scores from a principle component.}
\usage{
getPCRanks(mat, IDs, PC, filter_thresh)
}
\description{
Returns eigenvector scores for input CpG sites.}
\arguments{
\item{mat}{Input data frame.}
\item{IDs}{Character vector of treatment names.}
\item{PC}{Integer, principle component to use.}
\item{filter_thresh}{Integer, Coverage filtering threshold.}
}
\value{Bismark2Matrix.R output file, or data frame object.}
\value{A character vector of IDs containing the common names for compared conditions. E.g., for samples trt1 & trt2 vs. ctl1 & ctl2, IDs=c("trt1", "ctl")}
\value{Which principle component to use. Use to DefineBestPC if unsure.}
\value{A coverage threshold for filtering, where CpG coverage of all samples must be larger than this value.}
\examples{
getPCRanks(mat, c("trt", "ctl"), PC=1)
}
