\name{fdata}
\alias{as.fdata}
\alias{as.fdata.default}
\alias{as.fdata.matrix}
\alias{as.fdata.array}
\alias{as.fdata.list}
\alias{fdata}
\alias{print.fdata}
\alias{summary.fdata}
\alias{print.summary.fdata}
\title{ Functional Data class }
\description{
  Object of class 'fdata' and its methods.
}
\usage{
as.fdata(object,...)
as.fdata.matrix(object,..., col, p, dates, name)
as.fdata.list(object,..., dates, name)
}
\arguments{
  \item{object}{ A matrix or a list. }
  \item{col}{ A vector giving the names of the variables to include in
  the 'fdata' object. }
  \item{p}{ A real value giving the number of discretization point chosen.}
  \item{dates}{ A vector of character containing the dates of the observations.}
  \item{name}{ A vector of character containing the names of the
  variables (generated if not provided). }
  \item{\dots}{ Additional arguments. }
}
\details{
  Fdata objects are mainly used to modelize functional data in the purpose
  of computing functional autoregressive model by the
  \code{\link{far}} and \code{\link{kerfon}} functions.

  An fdata is composed of one or several variables.
  Each ones is a functional time series.

  To be more precise, every variable got a functional data by
  element of the \code{dates} (explicitly given or implicitly
  deduced). So the number of functional observations is a common data.
  
  In the contrary, each variable can be expressed in a different
  functional space. For example, if you got two variables,
  Temperature and Wind, measured during 30 days. Choosing a daily
  representation, the \code{fdata} will contain a 30 elements long
  \code{dates} vector. Nevertheless, the variables measurement can be
  different. If Temperature is measured every hour and Wind every two
  hours, the \code{fdata} object can handle such a representation.
  The only constraint is to get a regular measurement: no changes in the
  methodology.
  
  Basically, the \code{fdata} objects are discrete measurements but the
  modelization which can be used on it will make it functional.
  Indeed, The first methods implemented as \code{far} and \code{kerfon}
  use a linear approximation, but more sophisticate modelization, as
  splines or wavelets approximations may come.
}
\value{
  An object of class fdata.
}
\author{ J. Damon }

\seealso{ \code{\link{far}}, \code{\link{multplot}},
  \code{\link{maxfdata}}, \code{\link{kerfon}}.}

\examples{
# Reading of the data
library(stats)
data(UKDriverDeaths)

# Making the data of class 'fdata'
fUKDriverDeaths <- as.fdata(UKDriverDeaths,col=1,p=12,dates=1969:1984,
                            name="UK Driver Deaths")
summary(fUKDriverDeaths)

# ploting of the data : whole and 1 year
par(mfrow=c(2,1))
plot(fUKDriverDeaths,xval=1969+(1:192)/12,whole=TRUE,
     name="Whole Evolution : ")
plot(fUKDriverDeaths,date="1984",xval=1:12,
     name="Evolution during year 1984 : ")

# Matrix conversion
print(as.fdata(matrix(rnorm(50),10,5)))
print(as.fdata(matrix(rnorm(500),100,5),col=1:2,p=5))

# List Conversions
print(as.fdata(list("X"=matrix(rnorm(100),10,10),
"Z"=matrix(rnorm(50),5,10))))
}
\keyword{ ts }
\keyword{ misc }
