\name{Get_Novel_DMRs}
\alias{Get_Novel_DMRs}
\title{Call DMRs from WGBS data.}
\usage{
Get_Novel_DMRs(ranks, nSeeds, chromDictObj, DMR_resolution, QueryLimit, minCpGs, minZ, perms)
}
\description{
DMR Calling.}
\arguments{
\item{ranks}{Rank data frame from getPCRanks.}
\item{nSeeds}{Integer, number of input seeds for DMR expansion.}
\item{chromDictObj}{chromDict output. If null, chromDict() is run internally.}
\item{DMR_resolution}{Automatic=NULL. Integer, expansion size, default=QueryLimit/15.}
\item{QueryLimit}{Automatic=5000. Maximum DMR expansion size.}
\item{minCpGs}{Minimum CpGs in a DMR region.}
\item{minZ}{Automatic=1. Absolute Z score threshold for DMR calling.}
\item{perms}{Automatic=1000. Number of permutations to use when defining the null distribution.}
}
\value{Rank data frame from getPCRanks.}
\value{Numeric. Number of starting seeds for DMR expansion.}
\value{Output of chromDict. Runs internally if NULL.}
\value{Automatic=NULL. Integer, number of bases to increase the DMR by with each expansion. If NULL, QueryLimit/15.}
\value{Maximum DMR expansion size. Not recommended to play with this value, script is optimized for 5000bp maximum size.}
\value{Removes DMRs with fewer CpGs than minCpGs.}
\value{Absolute Z score threshold for DMR calling; internal value. Not recommended to play with this setting.}
\value{Number of permutations to use when defining the null distribution. Increasing this value largely influences computational time with minimal return.}
\examples{
find_best_nSeed(ranks, 1000)
}

