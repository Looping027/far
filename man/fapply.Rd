\name{fapply}
\alias{fapply}
\title{ Apply functions over a fdata object }
\description{
  \code{fapply} returns a \code{fdata} object of the same length as
  data. Each element of which is the result of applying \code{FUN} to
  the corresponding element of data.
}
\usage{
fapply(data, FUN, row.names, ...)
}
\arguments{
  \item{data}{ A \code{fdata} object }
  \item{FUN}{ the function to be applied. In the case of functions like
  +, \%*\%, etc., the function name must be quoted. }
  \item{row.names}{ a vector giving the names describing the results of
  \code{FUN} }
  \item{\dots}{ optional arguments to \code{FUN}. }
}
\details{
  This function has to be used only with \code{fdata} objects, unless it
  stop, returning no value.
}
\value{
  The returned value is a \code{fdata} object too.
}
\author{ J. Damon }


\seealso{ \code{\link{apply}}, \code{\link{lapply}}.}

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

  fapply(data1,sum)
  multplot(fapply(fapply(data1,abs),cumsum))
}

\keyword{ univar }
