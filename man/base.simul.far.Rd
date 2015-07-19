\name{base.simul.far}
\alias{base.simul.far}
\title{ Creating functional basis }
\description{
  Computation of a particular basis in a functional space.
}
\usage{
base.simul.far(m=24, n=5)
}
\arguments{
  \item{m}{ Number of discretization points }
  \item{n}{ Number of axis }
}
\details{
  We consider a sinusoidal basis of the functional space C[0;1] of the
  continuous functions from [0;1] to R. We compute here the values of
  the \code{n} first (functional) axis at \code{m} equi-repartited
  discretization points in [0;1] (more precisely the point
  0,\eqn{\frac{1}{\code{m}}}{1/\code{m}},...,
  \eqn{\frac{\code{m}-1}{\code{m}}}{(\code{m}-1)/\code{m}}).
}
\value{
  A matrix of size \code{m} x \code{n} containing the \code{m} values of
  the \code{n} first axis of the basis.
}
\author{ J. Damon }

\note{
  The chosen basis is orthogonal.
  
  The aim of this function is to provide an internal tool for the
  function \code{\link{simul.farx}}.
}
\seealso{ \code{\link{simul.farx}} }

\examples{
  print(temp<-base.simul.far(10,3))
  print(t(temp)\%*\%temp)
  matplot(base.simul.far(100,5),type='l')
}
\keyword{ algebra }
