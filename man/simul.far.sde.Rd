\name{simul.far.sde}
\alias{simul.far.sde}
\title{ FAR-SDE process simulation }
\encoding{utf-8}
\description{
  Simulation of a FAR process following an Stochastic
  Differential Equation
}
\usage{
simul.far.sde(coef=c(0.4, 0.8), n=80, p=32, sigma=1)
}
\arguments{
  \item{coef}{ Numerical vertor. It contains the two values of the
    coefficients (\eqn{a_1}{a(1)} and \eqn{a_2}{a(2)}, see details for
    more informations). }
  \item{n}{ Integer. The number of observations generated. }
  \item{p}{ Integer. The number of discretization points. }
  \item{sigma}{ Numeric. The standard deviation (see details for more
    informations). }
}
\details{
  This function implements the simulation proposed by Besse and Cardot
  (1996) to simulate a FAR process following the Stochastic Differential
  Equation:

  \deqn{dX^{(2)}+a_2.dX+a_1.X=\code{sigma}.dW}{%
    d^2(X)+a(2).d(X)+a(1).X=\code{sigma}.d(W)}

  Where \eqn{dX^{(2)}}{d^2(X)} and \eqn{dX}{d(X)} stand respectively for
  the second and first derivate of the process X, and W is a brownian
  process.

  The coefficients \eqn{a_1}{a(1)} and \eqn{a_2}{a(2)} are the two first
  elements of \code{coef}.

  The simulation use a order one approximation inspired by the work of
  Milstein, as described in Besse and Cardot (1996).
}
\value{
  A \code{fdata} object containing one variable ("var") which is a
  FAR(1) process of length \code{n} with \code{p} discretization
  points.
}
\references{ Besse, P. and Cardot, H. (1996).
  \emph{Approximation spline de la prévision d'un processus
  fonctionnel autorégressif d'ordre 1}.
  Revue Canadienne de Statistique/Canadian Journal of
  Statistics, \bold{24}, 467--487. }
\author{ J. Damon }

\seealso{ \code{\link{simul.far}}, \code{\link{simul.far.wiener}},
\code{\link{simul.farx}}, \code{\link{simul.wiener}}.  }

\examples{
  far1 <- simul.far.sde()
  summary(far1)
  print(far(far1,kn=2))
  par(mfrow=c(2,1))
  plot(far1,date=1)
  plot(select.fdata(far1,date=1:5),whole=TRUE,separator=TRUE)

}
\keyword{ ts }
\keyword{ misc }
