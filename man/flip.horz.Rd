% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/primitive.matrix.functions.R
\name{flip.horz}
\alias{flip.horz}
\title{Flip matrix horizontally}
\usage{
flip.horz(x)
}
\arguments{
\item{x}{a 2d-matrix}
}
\value{
the given matrix flipped from left to right
}
\description{
Flip given matrix on a vertical axis, returning it mirrored from left to right.
}
\details{
\code{t(apply(x,1,rev))}
}
\examples{
a<-matrix(1:9,3,3)
a
#     [,1] [,2] [,3]
#[1,]    1    4    7
#[2,]    2    5    8
#[3,]    3    6    9
flip.horz(a)
#     [,1] [,2] [,3]
#[1,]    7    4    1
#[2,]    8    5    2
#[3,]    9    6    3
}

