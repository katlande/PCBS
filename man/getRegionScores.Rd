\name{score_metagene}
\alias{score_metagene}
\title{Make a metagene from PC Scores.}
\usage{
score_metagene(ranks, regions, bin, title, xaxis, yaxis, chromDictObj, return.data, linecol)
}
\description{Returns p-values and Z-scores for CpGs in a set of regions, compared to a local null background distribution.}
\arguments{
\item{ranks}{getPCRanks output data.frame}
\item{regions}{A three-column dataframe containing a set of regions to test. Columns = chrom, start, end.}
\item{bin}{Number of bins to use in metagenes. Default=100.}
\item{title}{Output plot title}
\item{xaxis}{Output plot x-axis title}
\item{yaxis}{Output plot y-axis title}
\item{chromDictObj}{Optional chromDictObject made from chromDict(), runs internally if set to NULL (default). Scripts that run this function multiple times will be sped up by setting this option}
\item{return.data}{Whether to return a plot, or data that can be run with plot_metagene() or multiple_metagenes().}
\item{linecol}{Line colour if return.data=F, auto="red"}
}

