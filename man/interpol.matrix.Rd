\name{interpol.matrix}
\alias{interpol.matrix}
\title{Interpolation matrix}
\description{
  Calculate the matrix giving the linear interpolation of regularly
  spaced points.
}
\usage{
interpol.matrix(n = 12, m = 24, tol = sqrt(.Machine$double.eps))
}

\arguments{
  \item{n}{ Number (integer) of points in output space }
  \item{m}{ Number (integer) of points in the input function (or space) }
  \item{tol}{ A relative tolerance to detect zero singular values. }
}
\details{
  The general principle is, considering a function for which we know
  values at \code{m} equally spaced points (for instance 1/\code{m},
  2/\code{m}, ..., 1), to compute the matrix giving the linear
  approximation of \code{n} equally spaced points (for instance
  1/\code{n}, 2/\code{n}, ..., 1).

  The function works whether \code{n} or \code{m} is the largest.

  The function is vectorized, so \code{m} and \code{n} can be vectors of
  integers. In this case, they have to be of the same size and the
  resulting matrix is block diagonal.
}
\value{
  A \code{n}x\code{m} matrix if they are integer, else a
  \code{sum(n)}x\code{sum(m)} matrix.
}
\author{J. Damon}

\seealso{ \code{\link{theoretical.coef}}, \code{\link{simul.far}} or
  \code{\link{simul.farx}}. }

\examples{
  mat1 <- interpol.matrix(12,24)
  mat2 <- interpol.matrix(c(3,5),c(12,12))
  print(mat1 \%*\% base.simul.far(24,5))
  print(mat2 \%*\% base.simul.far(24,5))  
}
\keyword{algebra}
