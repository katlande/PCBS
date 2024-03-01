\name{rankDist}
\alias{rankDist}
\title{Identify the best rank cut-off for significant CpGs.}
\usage{
rankDist(ranks, draw_intersects, noise_perc, mode, return.plot)
}
\description{
Automated rank cut-off estimator for input CpGs.}
\arguments{
\item{ranks}{Input data frame.}
\item{draw_intersects}{T/F.}
\item{noise_perc}{Automatic=0.5, numeric between 0 and 1.}
\item{mode}{"intersect" or "strict"}
\item{return.plot}{T/F.}
}
\value{getPCRanks output data frame.}
\value{T/F whether to draw intersect lines if return.plot=T.}
\value{Fraction of ranks to use to model the background noise. Not recommended to play with this value. Increasing/decreasing returns a looser/stricter threshold, respectively.}
\value{Determine cut-off with "intersect" or "strict" method. "Strict" is recommended for sets with lower variability.}
\value{T/F, whether to return a plot or a numeric.}
\examples{
rankDist(ranks, mode="intersect")
rankDist(ranks, mode="strict")
}
