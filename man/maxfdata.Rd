\name{maxfdata}
\alias{maxfdata}
\title{ Maxima of functional data }
\description{
  Extract the maxima series from a functional data object.
}
\usage{
maxfdata(data)
}
\arguments{
  \item{data}{ A \code{fdata} object }
}
\value{
  A \code{fdata} object.
}
\author{ J. Damon }

\seealso{ \code{\link{fapply}} }

\examples{
  # Simulation of a FARX process
  data1 <- simul.farx(m=10,n=400,base=base.simul.far(20,5),
                base.exo=base.simul.far(20,5),
                d.a=matrix(c(0.5,0),nrow=1,ncol=2),
                alpha.conj=matrix(c(0.2,0),nrow=1,ncol=2),
                d.rho=diag(c(0.45,0.90,0.34,0.45)),
                alpha=diag(c(0.5,0.23,0.018)),
                d.rho.exo=diag(c(0.45,0.90,0.34,0.45)),
                cst1=0.0)

  print(data2 <- maxfdata(data1))
  print(unclass(data2))
}
\keyword{ ts }
\keyword{ univar }
