\name{simul.wiener}
\alias{simul.wiener}
\title{ Wiener process simulation }
\encoding{utf-8}
\description{
  Simulation of Wiener processes.
}
\usage{
simul.wiener(m=64, n=1, m2=NULL)
}
\arguments{
  \item{m}{ Integer. Number of discretization points. }
  \item{n}{ Integer. Number of observations. }
  \item{m2}{ Integer. Length of the Karhunen-Loève expansion (2\code{m}
    by default). }
}
\details{
  This function use the known Karhunen-Loève expansion of Wiener
  processes to simulate observations of such a process.

  The option \code{m2} is internally used to set the length of the
  expansion. This expansion need to be larger than the number of
  discretization points, but a too important value may slow down the
  generation. The default value as been chosen as a compromise.
}
\value{
  A \code{fdata} object containing one variable ("var") which is a Wiener
  process of length \code{n} with \code{m} discretization points.
}
\author{ J. Damon }
\references{ Pumo, B. (1992). \emph{Estimation et Prévision de Processus
    Autoregressifs Fonctionnels. Applications aux Processus à Temps Continu.}
    PhD Thesis, University Paris 6, Pierre et Marie Curie. }
\seealso{ \code{\link{simul.far.sde}}, \code{\link{simul.far.wiener}},
\code{\link{simul.farx}}, \code{\link{simul.far}}.  }

\examples{
  noise <- simul.wiener(m=64,n=100,m2=512)
  summary(noise)
  par(mfrow=c(2,1))
  plot(noise,date=1)
  plot(select.fdata(noise,date=1:5),whole=TRUE,separator=TRUE)
}
\keyword{ ts }
\keyword{ misc }
