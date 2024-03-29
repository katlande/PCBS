\name{find_best_nSeed}
\alias{find_best_nSeed}
\title{Identify the seed number to use for DMR calling.}
\usage{
find_best_nSeed(ranks, CpG_cutoff, values, max.dmr.size, return.plot)
}
\description{
Identifies the optimal seed number to use for DMR calling.}
\arguments{
\item{ranks}{Rank data frame from getPCRanks.}
\item{CpG_cutoff}{NULL or numeric.}
\item{values}{Numeric vector.}
\item{max.dmr.size}{Automatic=5000. Maximum DMR expansion size.}
\item{return.plot}{Boolean, return plot of numeric.}
}
\value{Rank data frame from getPCRanks.}
\value{If NULL, seed numbers tested will be input of the values argument. If numeric, seed numbers tested will be CpG_cutoff*values argument. Recommended to us rankDist estimate if not null}
\value{Numeric vector, either seed numbers to test if CpG_cutoff=NULL or multipliers if CpG_cutoff is numeric.}
\value{Maximum DMR expansion size in downstream analysis. Note: pipeline is optimized for 5000bp max DMR size, it is not recommended to play with this value.}
\value{T/F, whether to return a plot or a numeric representing the best seed number for downstream analysis.}
\examples{
find_best_nSeed(ranks, 1000)
}

