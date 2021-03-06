% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/primitive.matrix.functions.R
\name{rotate.right}
\alias{rotate.right}
\title{Rotate matrix 90 degrees clockwise}
\usage{
rotate.right(x)
}
\arguments{
\item{x}{a 2d-matrix}
}
\value{
the given matrix rotated 90 degrees clockwise
}
\description{
Rotate matrix 90 degrees clockwise
}
\details{
\code{t(apply(x,2,rev))}
}
\examples{
a<-matrix(1:9,3,3)
a
#     [,1] [,2] [,3]
#[1,]    1    4    7
#[2,]    2    5    8
#[3,]    3    6    9
rotate.right(a)
#     [,1] [,2] [,3]
#[1,]    3    2    1
#[2,]    6    5    4
#[3,]    9    8    7
}

