\name{chromDict}
\alias{chromDict}
\title{Convert a rank object into a chromDict.}
\usage{
chromDict(ranks)
}
\description{
Internal to many functions; creates a chromDict for faster computing times. chromDict can be run separately to speed up functions run iteratively.
A chromDict is a list of chromosome-specific data.tables generated from ranks.}
\arguments{
\item{ranks}{getPCRanks output data frame}
}
\value{
Returns a list of \code{data.tables} for each chromosome, for faster analysis. Used internall by many PCBS functions.
}
\examples{
ranks <- getPCRanks(eigen, IDs = c("trt", "ctl"), PC = 1)
chromDictObj <- chromDict(ranks)
}

