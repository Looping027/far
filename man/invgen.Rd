\name{invgen}
\alias{invgen}
\title{ Generalized inverse of a Matrix }
\description{
  Calculates the Moore-Penrose generalized inverse of a matrix X.
}
\usage{
invgen(a, tol)
}
\arguments{
  \item{a}{ Matrix for which the Moore-Penrose inverse is required. }
  \item{tol}{ A relative tolerance to detect zero singular values. }
}
\value{
  A Moore-Penrose generalized inverse matrix for X.
}
\seealso{ \code{\link{solve}},\code{\link{svd}},\code{\link{eigen}} }
\examples{
  mat1<-matrix(rnorm(100),ncol=10)
  print(invgen(mat1))
}
\keyword{ algebra }
