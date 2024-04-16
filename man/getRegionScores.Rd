\name{getRegionScores}
\alias{getRegionScores}
\title{Calculated methylation significance in a set of regions.}
\usage{
getRegionScores(ranks, regions)
}
\description{Returns p-values and Z-scores for CpGs in a set of regions, compared to a local null background distribution.}
\arguments{
\item{ranks}{getPCRanks output data.frame}
\item{regions}{A three-column dataframe containing a set of regions to test. Columns = chrom, start, end.}
}

