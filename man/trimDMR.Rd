\name{trimDMR}
\alias{trimDMR}
\title{
Nested DMR calling function within within Get_Novel_DMRs()
}
\description{
Trims the edges off of DMR expansions.
}
\usage{
trimDMR(df, region, min.dmr.cpgs, max.dmr.size, null_summary, null_values)
}

\arguments{
\item{df}{DMR expansion output dataframe.}
\item{region}{get_all_DMRs() region output}
\item{min.dmr.cpgs}{Get_Novel_DMRs arg minCpGs}
\item{max.dmr.size}{Get_Novel_DMRs arg QueryLimit}
\item{null_summary}{null distribution conrtainer}
\item{null_values}{null distribution conrtainer}
}
