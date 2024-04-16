\name{eigen}
\alias{eigen}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("eigen")}
\format{
  A data frame with 50000 observations on the following 13 variables.
  \describe{
    \item{\code{cpgID}}{a character vector, chrom:locus}
    \item{\code{trt1_PercMeth}}{a numeric vector, percent methylated of sample trt1}
    \item{\code{trt1_nCpG}}{a numeric vector, depth of sample trt1}
    \item{\code{trt2_PercMeth}}{a numeric vector, percent methylated of sample trt2}
    \item{\code{trt2_nCpG}}{a numeric vector, depth of sample trt2}
    \item{\code{trt3_PercMeth}}{a numeric vector, percent methylated of sample trt3}
    \item{\code{trt3_nCpG}}{a numeric vector, depth of sample trt3}
    \item{\code{ctl1_PercMeth}}{a numeric vector, percent methylated of sample ctl1}
    \item{\code{ctl1_nCpG}}{a numeric vector, depth of sample ctl1}
    \item{\code{ctl2_PercMeth}}{a numeric vector, percent methylated of sample ctl2}
    \item{\code{ctl2_nCpG}}{a numeric vector, depth of sample ctl2}
    \item{\code{ctl3_PercMeth}}{a numeric vector, percent methylated of sample ctl3}
    \item{\code{ctl3_nCpG}}{a numeric vector, depth of sample ctl3}
    
  }
}

\source{
Simulated WGBS data
}

\examples{
head(eigen)
}
\keyword{datasets}
