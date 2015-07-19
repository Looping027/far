\name{plot.fdata}
\alias{plot.fdata}
\title{ Plot Functional Data }
\description{
  Plot Functional Data (more precisely \code{fdata} objects).
}
\usage{
\method{plot}{fdata}(x,...,date, xval, name, main, whole, separator)
}
\arguments{
  \item{x}{ A \code{fdata} object. }
  \item{date}{ A vector of character giving the chosen dates. }
  \item{xval}{ A numerical vector giving the values to appear on the x
  axis. }
  \item{name}{ A vector of character giving the plotted variables. }
  \item{main}{ an overall title for the plot. }
  \item{whole}{ Logical. If \code{TRUE} all the observations are plot on
    the same graphic. }
  \item{separator}{ Logical. It will be used when
    \code{whole}=\code{TRUE}. If \code{TRUE} then dashed lines are
    plotted to separated the observations. }
  \item{\dots}{ Additional arguments to the plot. }
}
\details{
  This function facilitate the plotting of \code{fdata} objects.
  It is dedicated to univariate plots, please take a look at
  \code{\link{multplot}} if you need multivariate plots in one graphic.

  The default behaviour is to plot the observation called "1" of all the
  variables available in \code{x} (so it will produce as many plots as
  the number of variables).

  If you want less variables, use the \code{name} argument. If you need
  more observations, use the \code{date} argument. When provided, the
  \code{xval} argument allow you to change the labels of the x-axis.

  It is also possible to plot the complete series on the same plot using
  the \code{whole} argument. In this case, the \code{separator} allow
  you to draw line to distinguish the different observations of the
  functional data.
}
\author{ J. Damon }

\seealso{ \code{\link{fdata}}, \code{\link{multplot}}.}

\examples{
  # Reading of the data
  library(stats)
  data(UKDriverDeaths)

  # Making the data of class 'fdata'
  fUKDriverDeaths <- as.fdata(UKDriverDeaths, col=1, p=12,
                            dates=1969:1984,
                            name="UK Driver Deaths")
  summary(fUKDriverDeaths)

  # plotting of the data : whole and 1 year
  par(mfrow=c(2,1))
  plot(fUKDriverDeaths, xval=1969+(1:192)/12,
       whole=TRUE, name="Whole Evolution : ", separator=TRUE)
  plot(fUKDriverDeaths, date="1984", xval=1:12,
       name="Evolution during year 1984 : ")
}
\keyword{ ts }
\keyword{ aplot }
