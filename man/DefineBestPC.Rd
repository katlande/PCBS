\name{DefineBestPC}
\alias{DefineBestPC}
\title{Identify your best principle component.}
\usage{
DefineBestPC(mat, IDs, filter_thresh, return.plot)
}
\description{
Defines the best principle component to use for downstream analysis.}
\arguments{
\item{mat}{Input file.}
\item{IDs}{Character vector of treatment names.}
\item{filter_thresh}{Coverage filtering threshold.}
\item{return.plot}{T/F, return plot or numeric.}
}
\value{Bismark2Matrix.R output file, or data frame object.}
\value{A character vector of IDs containing the common names for compared conditions. E.g., for samples trt1, trt2 vs. ctl1, ctl2, IDs=c("trt", "ctl")}
\value{A coverage threshold for filtering, where CpG coverage of all samples must be larger than this value.}
\value{T/F, whether to return a PCA plot or a numeric representing the best principle component for downstream analysis.}
\examples{
DefineBestPC(mat, c("trt", "ctl"), 50, T)
}

