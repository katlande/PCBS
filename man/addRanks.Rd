\name{addRanks}
\alias{addRanks}
\title{Add ranks to eigenvector scores.}
\usage{
addRanks(ranks)
}
\description{
Defines the best principle component to use for downstream analysis.}
\arguments{
\item{ranks}{getPCRanks output data frame.}
}
\value{
The input \code{data.frame} with rank order and absolute rank order columns.
}
\examples{
ranks <- getPCRanks(eigen, IDs = c("trt", "ctl"))
ranks <- addRanks(ranks)
}

