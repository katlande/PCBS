\name{score_metagene}
\alias{score_metagene}
\title{Make a metagene from PC Scores.}
\usage{
score_metagene(ranks, regions, bin, title, xaxis, yaxis, chromDictObj, return.data, linecol)
}
\description{Uses mean binned PC scores across a set of regions to draw a metagene.}
\arguments{
\item{ranks}{data.frame}
\item{regions}{data.frame}
\item{bin}{integer}
\item{title}{character}
\item{xaxis}{character}
\item{yaxis}{character}
\item{chromDictObj}{chromDictObj}
\item{return.data}{boolean}
\item{linecol}{character}
}
\value{getPCRanks output.}
\value{A three-column dataframe containing a set of regions to test. Columns = chrom, start, end.}
\value{Number of bins to use in metagenes. Default=100.}
\value{Output plot title.}
\value{Output plot x-axis title}
\value{Output plot y-axis title.}
\value{Optional chromDictObject made from chromDict(), runs internally if set to NULL (default). Scripts that run this function multiple times will be sped up by setting this option.}
\value{Whether to return a plot, or data that can be run with plot_metagene() or multiple_metagenes().}
\value{Colour for line, auto="red".}
\examples{
score_metagene(ranks, regions, chromDictObj=chromDictObj)
}
