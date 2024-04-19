\name{se}
\alias{se}
\title{Standard error of a vector.}
\usage{
se(x)
}
\description{Takes the standard error of a vector.}
\arguments{
\item{x}{numeric vector.}
}
\value{Returns a \code{numeric}, representing the standard error of the input vector.} 
\examples{
x <- sample(1:100, 20)
se(x)
}

