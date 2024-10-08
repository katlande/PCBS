\name{getRegionScores}
\alias{getRegionScores}
\title{Calculated methylation significance in a set of regions.}
\usage{
getRegionScores(ranks, regions,  chromDictObj)
}
\description{Returns p-values and Z-scores for CpGs in a set of regions, compared to a local null background distribution.}
\arguments{
\item{ranks}{getPCRanks output data.frame, only necessary if chromDictObj=NULL}
\item{regions}{A three-column dataframe containing a set of regions to test. Columns = chrom, start, end.}
\item{chromDictObj}{chromDict() output object, recommended input instead of ranks.}
}
\value{Returns a \code{data.frame} with significance scores for all input regions.}
\examples{
ranks <- getPCRanks(eigen, IDs = c("trt", "ctl"), PC = 1)
DMLs <- addRanks(ranks)

# data.frame of regions to test:
regions <- data.frame(chr=c("chr3", "chr3", "chr1"),
                      s=c(4920450, 3961576, 300000),
                      e=c(4923267, 3963805, 302900),
                      ID=c("Hypo-DMR", "partial Hyper-DMR", "random"))
                
getRegionScores(DMLs, regions)
}
