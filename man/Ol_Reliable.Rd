\name{Ol_Reliable}
\alias{Ol_Reliable}
\title{
PCBS ggplot theme.
}
\description{
Custom theme for ggplot used by all PCBS output figures.
}
\usage{
Ol_Reliable()
}


\examples{
df <- data.frame(A=c(1,2,3), B=c(1,2,3))

ggplot2::ggplot(df, ggplot2::aes(x=A, y=B))+
ggplot2::geom_point()+
Ol_Reliable()

}

