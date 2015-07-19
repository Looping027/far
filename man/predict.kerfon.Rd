\name{predict.kerfon}
\alias{predict.kerfon}
\title{ Forecasting of functional kernel model }
\description{
  Computation of the prediction based on a functional kernel model
}
\usage{
\method{predict}{kerfon}(object, ..., newdata=NULL, label, na.rm=TRUE, positive=FALSE)
}
\arguments{
  \item{object}{ A \code{kerfon} object result of the \code{kerfon} function. }
  \item{newdata}{ A \code{\link{fdata}} object used in the
    \code{kerfon} model to compute the prediction, or \code{NULL}
    to predict one step forward with the data in \code{object}. }
  \item{label}{ A vector of character giving the dates to associate to
  the predicted observations. }
  \item{na.rm}{ Logical. Does the \code{n.a.} need to be removed. }
  \item{positive}{ Logical. Does the result must be forced to positive values. }
  \item{\dots}{ Additional arguments. }
}
\details{
  This function computes one step forward prediction for a
  \code{\link{kerfon}} model.

  Use the \code{newdata} option to input the past values,
  and the \code{label} option value to define the labels for the new
  observations. Notices that the output as the same length as
  \code{newdata}.

  In some special context, the user may need to suppress the
  \code{na.rm} observations with the \code{na.rm} option, or force the
  prediction to be positive with the \code{positive} option (in this
  case the result will be maximum of 0 and the predicted value).
}
\value{
  A \code{fdata} object.
}
\author{ J. Damon }

\seealso{ \code{\link{kerfon}} }

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

  # Cross validation
  model1 <- kerfon(data=data1, x="X", r=10, na.rm=TRUE)
	print(model1)

  # Predicting values
  pred1 <- predict(model1,newdata=select.fdata(data1,date=1:399))
  # Persistence
  persist1 <- pred.persist(select.fdata(data1,date=1:399),x="X")
  # Real values
  real1 <- select.fdata(data1,date=2:400)

  errors0 <- persist1[[1]]-real1[[1]]
  errors1 <- pred1[[1]]-real1[[1]]

  # Norm of observations
  summary(real1)
  # Persistence
  summary(as.fdata(errors0))
  # kerfon model
  summary(as.fdata(errors1))
}
\keyword{ models }
\keyword{ ts }
