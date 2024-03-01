\name{checkRank}
\alias{checkRank}
\title{Check k-values manually.}
\usage{
checkRank(ranks, cutoff)
}
\description{Plots a score vs. rank plot with a manually chosen rank cut-off for manual k selection.}
\arguments{
\item{ranks}{data.frame}
\item{cutoff}{integer}
}
\value{getPCRanks output data frame.}
\value{integer, rank value to check}
\examples{
checkRank(ranks, 10000)
}
