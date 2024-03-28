\name{score_metagene}
\alias{score_metagene}
\title{Make a metagene from PC Scores.}
\usage{
score_metagene(ranks, regions, bin, title, xaxis, yaxis, chromDictObj, return.data, linecol)
}
\description{Uses mean binned PC scores across a set of regions to draw a metagene.}
\arguments{
\item{ranks}{getPCRanks output data.frame}
\item{regions}{A three-column data.frame containing a set of regions to test. Columns = chrom, start, end.}
\item{bin}{integer, number of bins to use in metagenes. Default=100.}
\item{title}{Output plot title}
\item{xaxis}{Output plot x-axis title}
\item{yaxis}{Output plot y-axis title}
\item{chromDictObj}{Optional chromDictObject made from chromDict(), runs internally if set to NULL (default). Scripts that run this function multiple times will be sped up by setting this option.}
\item{return.data}{T/F, whether to return a plot, or data that can be run with plot_metagene() or multiple_metagenes().}
\item{linecol}{Colour for line, auto="red"}
}

