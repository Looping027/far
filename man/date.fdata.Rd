\name{date.fdata}
\alias{date.fdata}
\title{ Extract the date of fdata }
\description{
  Extract the date(s) of \code{fdata} objects
}
\usage{
date.fdata(data)
}
\arguments{
  \item{data}{ A \code{fdata} object }
}
\details{
  The dates are the labels of the functionals observations of the
  \code{fdata} object.

  \code{fdata} are not constructed as \code{ts} object so a specific
  function to obtain the date is useful.
}
\value{
  A vector giving the dates (as character).
}
\author{ J. Damon }

\seealso{ \code{\link{fdata}} }

\examples{
# Reading the data
library(stats)
data(UKDriverDeaths)

# Conversion of the data
fUKDriverDeaths <- as.fdata(UKDriverDeaths,col=1,p=12,dates=1969:1984,
                            name="UK Driver Deaths")
date.fdata(fUKDriverDeaths)
}
\keyword{ ts }
\keyword{ misc }
