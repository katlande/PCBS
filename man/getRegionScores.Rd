\name{getRegionScores}
\alias{getRegionScores}
\title{Get differential methylation scores for a set of input regions.}
\usage{
getRegionScores(ranks, regions)
}
\description{Returns p-values and Z-scores for CpGs in a set of regions, compared to a local null background distribution.}
\arguments{
\item{ranks}{data.frame}
\item{regions}{data.frame}
}
\value{addRanks output dataframe.}
\value{A four-column dataframe containing a set of regions to test. Columns = chrom, start, end, region_ID.}
\examples{
getPCRanks(ranks)
}
