\name{DefineBestPC}
\alias{DefineBestPC}
\title{Identify your best principle component.}
\usage{
DefineBestPC(mat, IDs, filter_thresh, return.plot)
}
\description{
Defines the best principle component to use for downstream analysis.}
\arguments{
\item{mat}{Bismark2Matrix.R output file, or data frame object}
\item{IDs}{A character vector of IDs containing the common names for compared conditions. E.g., for samples trt1, trt2 vs. ctl1, ctl2, IDs=c("trt", "ctl")}
\item{filter_thresh}{A coverage threshold for filtering, where CpG coverage of all samples must be larger than this value}
\item{return.plot}{T/F, whether to return a PCA plot or a numeric representing the best principle component for downstream analysis}
}
\value{
If \code{return.plot=T}, a \code{grob} plotting a PCA of percent methylation of all samples is returned. Otherwise, a \code{numeric} representing the best principal component to use for PCBS analysis is returned.
}
\examples{
DefineBestPC(eigen, IDs = c("trt", "ctl"))
}

