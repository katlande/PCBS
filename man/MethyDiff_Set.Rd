\name{MethyDiff_Set}
\alias{MethyDiff_Set}
\title{Add a mean methylation column to a data.frame of regions.}
\usage{
MethyDiff_Set(chromDictMeth, regions)
}
\description{
Using a chromDictMeth() output object, quickly calculate the mean methylation difference across set of regions.}
\arguments{
\item{chromDictMeth}{chromDictMeth() output object}
\item{regions}{\code{data.frame} of regions, where column 1 = chromosome, column 2 = region start, and column 3 = region end}
}
\value{
Returns the regions object with a mean percent methylation column.
}
\examples{
chromDictMethylDiff <- chromDictMeth(eigen, c("trt", "ctl"))
regions <- data.frame(chr=c("chr3", "chr3", "chr1"),
                      s=c(4920450, 3961576, 300000),
                      e=c(4923267, 3963805, 302900),
                      ID=c("Hypo-DMR", "partial Hyper-DMR", "random"))

MethyDiff_Set(chromDictMethylDiff, regions)
}
