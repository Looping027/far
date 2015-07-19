\name{is.na.fdata}
\alias{is.na.fdata}
\title{ Not Available / ``Missing'' Values }
\description{
  The generic function is.na returns a logical vector of the same
  ``form'' as its argument x, containing TRUE for those elements marked
  NA or NaN (!) and FALSE otherwise. dim, dimnames and names attributes
  are preserved.
}
\usage{
\method{is.na}{fdata}(x)
}
\arguments{
  \item{x}{ A \code{fdata} object }
}
\details{
  An observation is considered as NA if any of its values is NA.
}
\value{
  A matrix of Logical values giving as rows the variables of \code{x}
  and as columns the observation.
}
\author{ J. Damon }

\seealso{ \code{\link{NA}} }

\examples{
  # Reading of the data
  library(stats)
  data(UKDriverDeaths)
  UKDriverDeaths[20]<-NA

  # Making the data of class 'fdata'
  fUKDriverDeaths <- as.fdata(UKDriverDeaths,col=1,p=12,dates=1969:1984,
                            name="UK Driver Deaths")
  summary(fUKDriverDeaths)
  is.na(fUKDriverDeaths)
}
\keyword{ NA }
