\name{methylDiff_metagene}
\alias{methylDiff_metagene}
\title{Make a metagene from mean percent methylation differences.}
\usage{
methylDiff_metagene(chromDictMethObj, regions, bin, title,
xaxis, yaxis, return.data, linecol, value)
}
\description{Uses mean binned percent methylation differences across a set of regions to draw a metagene.}
\arguments{
\item{chromDictMethObj}{chromDictMeth() output object}
\item{regions}{A three-column data.frame containing a set of regions to test. Columns = chrom, start, end.}
\item{bin}{integer, number of bins to use in metagenes. Default=100.}
\item{title}{Output plot title}
\item{xaxis}{Output plot x-axis title}
\item{yaxis}{Output plot y-axis title}
\item{return.data}{T/F, whether to return a plot, or data that can be run with plot_metagene() or multiple_metagenes().}
\item{linecol}{Colour for line, auto="red"}
\item{value}{Name of the plotted metric in chromDictMethObj. Only needs to be set explicity for abnormal use cases where chromDictMethObj contains a non-rank value output by chromDict().}
}
\value{If \code{return.data=F}, returns a \code{grob} containing a metagene plot. Otherwise, returns a \code{list} of two \code{data.frames} containing metagene and metagene standard error plotting information.}
\examples{

chromDictMethylDiff <- chromDictMeth(eigen, c("trt", "ctl"))
regions <- data.frame(chr=c("chr3", "chr3", "chr1"),
                      s=c(4920450, 3961576, 300000),
                      e=c(4923267, 3963805, 302900),
                      ID=c("Hypo-DMR", "partial Hyper-DMR", "random"))
 
methylDiff_metagene(chromDictMethylDiff, regions)
                     
}
