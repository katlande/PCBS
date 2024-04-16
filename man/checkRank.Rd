\name{checkRank}
\alias{checkRank}
\title{Check k-values manually.}
\usage{
checkRank(ranks, cutoff)
}
\description{Plots a score vs. rank plot with a manually chosen rank cut-off for manual k selection.}
\arguments{
\item{ranks}{getPCRanks output data frame}
\item{cutoff}{integer, rank value to check}
}

\examples{
ranks <- getPCRanks(eigen, IDs = c("trt", "ctl"), PC = 1)
test_50 <- checkRank(ranks, 50) # set cut-off to 50
test_500 <- checkRank(ranks, 500) # set cut-off to 500
}

