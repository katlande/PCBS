\name{chromDictMeth}
\alias{chromDictMeth}
\title{Create a chromDict of percent methylation difference at all sites.}
\usage{
chromDictMeth(mat, IDs, filter_thresh)
}
\description{
chromDicts are lists of keyed data.tables that enable very fast computing times. Much like the primary PCBS function chromDict() which makes a chromDict of all locus ranks, the chromDictMeth() function makes is a list of chromosome-specific data.tables containing percent methylation differences at all loci..}
\arguments{
\item{mat}{\code{data.frame} object containing percent methylation and locus information for all sites, in the format of \code{eigen}}
\item{IDs}{\code{character} vector of IDs containing the common names for compared conditions. E.g., for samples trt1 & trt2 vs. ctl1 & ctl2, IDs=c("trt1", "ctl")}
\item{filter_thresh}{\code{Integer}, a coverage threshold for filtering, where CpG coverage of all samples must be larger than this value. Auto=50}
}
\value{
Returns a list of \code{data.tables} for each chromosome, for faster analysis.
}
\examples{
chromDictMethylDiff <- chromDictMeth(eigen, c("trt", "ctl"))
}
