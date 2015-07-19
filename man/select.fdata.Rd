\name{select.fdata}
\alias{select.fdata}
\title{ Subscript of fdata }
\description{
  Use this function to subscript some functional observations of a
  functional data.
}
\usage{
select.fdata(data, date, name)
}
\arguments{
  \item{data}{ A \code{fdata} object. }
  \item{date}{ A vector of character containing the chosen dates (could
    be \code{NULL}). }
  \item{name}{ A vector giving the chosen name (could be \code{NULL}). }
}
\details{
  This function select one or several variables from \code{data}
  and can also subset the dates. This is useful in order to
  study the endogenous variables of a FARX process.
}
\value{
  A \code{fdata} object.
}
\author{ J. Damon }

\seealso{ \code{\link{fdata}} }

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

  print(data1)
  print(data1.X <- select.fdata(data1,name="X"))
  print(data2 <- select.fdata(data1,date=paste((1:5)*5)))
  date.fdata(data2)
}
\keyword{ manip }
