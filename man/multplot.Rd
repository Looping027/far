\name{multplot}
\alias{multplot}
\alias{multplot.default}
\alias{multplot.fdata}
\title{ Multivariate plots }
\description{
  Multivariate plots of Functional Data (more precisely \code{fdata}
  objects).
}
\usage{
multplot(object, ...)

\method{multplot}{fdata} (object, date = 1, xval = NULL, name = NULL, legend = FALSE,
    yleg, xlab = NULL, ylab = NULL, main = NULL, whole = FALSE, ...)
}
\arguments{
  \item{object}{ An \code{fdata} object for which a multplot is
    desired. }
  \item{date}{ String vector. List of the dates to work with. }
  \item{xval}{ Numerical vector. Values of the axis x. }
  \item{name}{ String vector. The set of variables to plot. }
  \item{legend}{ Boolean. Plot a legend ? }
  \item{yleg}{ Numeric. Where to put the legend box (y value). }
  \item{xlab}{ String. Title of the axis x. }
  \item{ylab}{ String. Title of the axis y. }
  \item{main}{ String. Title of the plot. }
  \item{whole}{ Boolean. A global plot (TRUE) or a plot by day (FALSE) }
  \item{\dots}{ Additional arguments. }
}
\details{
  This function facilitate the plotting of \code{fdata} objects.
  It is dedicated to multivariate plots, please take a look at
  \code{\link{plot.fdata}} if you need univariate plots in one graphic.

  The default behaviour is to produce one plot containing all the
  variables of the observation called "1".

  If you want less variables, use the \code{name} argument. If you need
  more observations, use the \code{date} argument. When provided, the
  \code{xval} argument allow you to change the labels of the x-axis.

  It is also possible to plot the complete series on the same plot using
  the \code{whole} argument.

  Moreover a legend facility is provided using the \code{legend} and
  \code{yleg} arguments.
}
\author{ J. Damon }

\seealso{ \code{\link{fdata}}, \code{\link{plot.fdata}}.}
\examples{
  # Simulation of a FARX process
  data1 <- simul.farx(m=10,n=100,base=base.simul.far(20,5),
                base.exo=base.simul.far(20,5),
                d.a=matrix(c(0.5,0),nrow=1,ncol=2),
                alpha.conj=matrix(c(0.2,0),nrow=1,ncol=2),
                d.rho=diag(c(0.45,0.90,0.34,0.45)),
                alpha=diag(c(0.5,0.23,0.018)),
                d.rho.exo=diag(c(0.45,0.90,0.34,0.45)),
                cst1=0.0)

 # 2 variables : X et Z
 # number of points per curve : 10
 # number of curves : 100
 # corresponding dates
 date.fdata(data1)

 multplot(data1) # plot the date "1" of the variables "X" and "Z"
 multplot(data1,legend=TRUE) # Same thing with a legend
 multplot(data1,legend=TRUE,yleg=-0.5) # same thing with a legend misplaced
 multplot(data1,main="day 1",legend=TRUE,xlab="hour",
          ylab="object of study")

 par(mfrow=c(1,3))
 multplot(data1,date=c("3","4","5")) # days "3", "4" and "5" are plotted
 par(mfrow=c(1,1))

 # to plot the whole series, we used whole = TRUE
 # but we have to give the x values

 multplot(data1,xval=seq(from=0,to=99.9,by=0.1),whole=TRUE)

 # to plot a subset of the series,
 # it is recommended to create a subset object with select.fdata
 data2 <- select.fdata(data1,date=c("4","5","6"))
 multplot(data2,xval=seq(from=4,to=6.9,by=0.1),whole=TRUE)
}
\keyword{ hplot }
