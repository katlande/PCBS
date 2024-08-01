\name{MethylDiff}
\alias{MethylDiff}
\title{Get the mean methylation difference across a specified region.}
\usage{
MethylDiff(chromDictMeth, chrom, start, end)
}
\description{
Using a chromDictMeth() output object, quickly calculate the mean methylation difference across a user-specified region.}
\arguments{
\item{chromDictMeth}{chromDictMeth() output object}
\item{chrom}{\code{character}, chromosome}
\item{start}{\code{integer}, region start}
\item{end}{\code{integer}, region end}
}
\value{
Returns a list of \code{data.tables} for each chromosome, for faster analysis.
}
\examples{
chromDictMethylDiff <- chromDictMeth(eigen, c("trt", "ctl"))
MethylDiff(chromDictMethylDiff, "chr3", 4920450, 4923267)
}
