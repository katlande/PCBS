\name{get_all_DMRs}
\alias{get_all_DMRs}
\title{
Nested DMR calling function within Get_Novel_DMRs()
}
\description{
Expands DMRs from collapsed seeds iteratively.
}
\usage{
get_all_DMRs(chromDictObj, seeds, res=40, max.dmr.size=3000, min.dmr.cpgs=10, min.absZscore, null)
}

\arguments{
\item{chromDictObj}{chromDict() output.}
\item{seeds}{compressed seeds}
\item{res=40}{Get_Novel_DMRs arg DMR_resolution}
\item{max.dmr.size}{Get_Novel_DMRs arg QueryLimit}
\item{min.dmr.cpgs}{Get_Novel_DMRs arg minCpGs}
\item{min.absZscore}{Get_Novel_DMRs arg minZ}
\item{null}{null distribution}
}
