\name{pred.persist}
\alias{pred.persist}
\title{ Forecasting using functional persistence }
\description{
  Compute prediction of functional data using the persistence.
}
\usage{
pred.persist(data, x, na.rm=TRUE, label, positive=FALSE)
}
\arguments{
  \item{data}{ A \code{fdata} object. }
  \item{x}{ A vector of character giving the names of the variables
    predicted. }
  \item{na.rm}{ Logical. Does the \code{n.a.} need to be removed. }
  \item{label}{ A vector of character giving the dates to associate to
  the predicted observations. }
  \item{positive}{ Logical. Does the result must be forced to positive
    values. }
}
\details{
  The persistence model is a beautiful way to name the simplest model
  ever. This model just suppose that the next observation will be equal
  to the previous one, that is to say, noting \eqn{\hat{X}_{n}}{Xp(n)}
  the prediction for \eqn{X_{n}}{X(n)} that we "compute" :
  
  \deqn{\hat{X}_{n+1}=X_{n}}{Xp(n+1)=X(n)}

  Of course, the intrinsic purpose of this model is to be a comparison
  for more complicated  models.

  The \code{x} option is provided to select the variable to predict,
  using the \code{label} option value as the labels for the new
  observations. Notices that the output as the same length as the input
  as it is only a shift in time.

  In some special context, the user may need to suppress the
  \code{na.rm} observations with the \code{na.rm} option, or force the
  prediction to be positive with the \code{positive} option (in this
  case the maximum of 0 and the past value is computed).
}
\value{
  A \code{fdata} object.
}
\author{ J. Damon }
\note{
  This has been more instinctive to call this function predict.persist
  but, due to the naming mechanism introduced by the object
  oriented programming, this would have reefer to the predict method for
  the persist objects. As it isn't the meaning of this function, we
  preferred the name pred.persist. 
}

\seealso{ \code{\link{predict.far}},\code{\link{predict.kerfon}}. }

\examples{
  # Simulation of a FARX process
  data1 <- simul.farx(m=10,n=40,base=base.simul.far(20,5),
                base.exo=base.simul.far(20,5),
                d.a=matrix(c(0.5,0),nrow=1,ncol=2),
                alpha.conj=matrix(c(0.2,0),nrow=1,ncol=2),
                d.rho=diag(c(0.45,0.90,0.34,0.45)),
                alpha=diag(c(0.5,0.23,0.018)),
                d.rho.exo=diag(c(0.45,0.90,0.34,0.45)),
                cst1=0.0)
  print(data2 <- pred.persist(data1,x="X",label="41"))
  print(unclass(select.fdata(data1,date=paste(38:40)))$X)
  print(unclass(select.fdata(data2,date=paste(39:41))))
}
\keyword{ ts }
\keyword{ misc }
