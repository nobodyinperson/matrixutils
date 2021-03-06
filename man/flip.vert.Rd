% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/primitive.matrix.functions.R
\name{flip.vert}
\alias{flip.vert}
\title{Flip matrix vertically}
\usage{
flip.vert(x)
}
\arguments{
\item{x}{a 2d-matrix}
}
\value{
the given matrix flipped upside down
}
\description{
Flip given matrix on a horizontal axis, thus returning it upside down.
}
\details{
\code{apply(x,2,rev)}
}
\examples{
a<-matrix(1:9,3,3)
a
#     [,1] [,2] [,3]
#[1,]    1    4    7
#[2,]    2    5    8
#[3,]    3    6    9
flip.vert(a)
#     [,1] [,2] [,3]
#[1,]    3    6    9
#[2,]    2    5    8
#[3,]    1    4    7
}

