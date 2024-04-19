\name{CheckOvercompression}
\alias{CheckOvercompression}
\title{Check if DMR calling seed number is ovcompressed.}
\usage{
CheckOvercompression(ranks, CpG_cutoff, values, max.dmr.size, return.plot)
}
\description{
Identifies if seed number to use for DMR calling causes overcompression.}
\arguments{
\item{ranks}{Rank data frame from getPCRanks.}
\item{CpG_cutoff}{NULL or numeric. If NULL, seed numbers tested will be input of the values argument. If numeric, seed numbers tested will be CpG_cutoff*values argument. Recommended to us rankDist estimate if not null}
\item{values}{Numeric vector, either seed numbers to test if CpG_cutoff=NULL or multipliers if CpG_cutoff is numeric}
\item{max.dmr.size}{Automatic=5000. Maximum DMR expansion size in downstream analysis. Note: pipeline is optimized for 5000bp max DMR size, it is not recommended to play with this value.}
\item{return.plot}{T/F, whether to return a plot or a numeric representing the best seed number for downstream analysis}
}
\value{
If \code{return.plot=T}, a \code{grob} plotting input seed number vs. compressed seed number is returned. Otherwise, a \code{numeric} is returned containing the largest tested input value without detectable overcompression.
}
\examples{
ranks <- getPCRanks(eigen, IDs = c("trt", "ctl"), PC = 1)
CheckOvercompression(ranks, 980)
}


