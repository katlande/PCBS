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
\item{chromDictObj}{chromDict() output. If null, chromDict() is run internally.}
\item{DMR_resolution}{Automatic=NULL. Integer, number of bases to increase the DMR by with each expansion. If NULL, QueryLimit/25.}
\item{QueryLimit}{Automatic=5000. Maximum DMR expansion size (bp)}
\item{minCpGs}{Automatic=15. Minimum CpGs in a DMR region, regions with fewer CpGs will be discarded.}
\item{minZ}{Automatic=1. Absolute Z score threshold for DMR calling; internal value. Not recommended to play with this setting.}
\item{perms}{Automatic=1000. Number of permutations to use when defining the null distribution. Increasing this value largely influences computational time with minimal return}
}

