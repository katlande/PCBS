\name{oneSeed}
\alias{oneSeed}
\title{
Nested DMR calling function within within get_all_DMRs(), Get_Novel_DMRs()
}
\description{
Expands one DMR from a single point.
}
\usage{
oneSeed(chroms, seed, resolution, max.size, mincpgs, null_list, Zlim=1)
}
\arguments{
\item{chroms}{chromDict() output.}
\item{seed}{seed value input}
\item{resolution}{Get_Novel_DMRs arg DMR_resolution}
\item{max.size}{Get_Novel_DMRs arg QueryLimit}
\item{mincpgs}{Get_Novel_DMRs arg minCpGs}
\item{null_list}{get_all_DMRs arg null}
\item{Zlim}{Get_Novel_DMRs arg minZ}
}
\value{returns a \code{list} of two indices, containing a \code{data.frame} with the output from a single compressed seed expansion and a \code{data.frame} of locus information around the seed expansion, intended for use within the \code{Get_Novel_DMRs()} function.}