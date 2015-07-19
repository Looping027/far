\name{predict.far}
\alias{predict.far}
\title{ Forecasting of FARX(1) model }
\description{
  Forecasting using FAR(1) or FARX(1) model
}
\usage{
\method{predict}{far}(object, ..., newdata=NULL, label, na.rm=TRUE, positive=FALSE)
}
\arguments{
  \item{object}{ A \code{far} object result of the \code{far} function. }
  \item{newdata}{ A data matrix (one column for each observation) used to
    predict the FAR(1) model from the values in newdata, or \code{NULL}
    to predict one step forward with the data in \code{object}. }
  \item{label}{ A vector of character giving the dates to associate to
  the predicted observations. }
  \item{na.rm}{ Logical. Does the \code{n.a.} need to be removed. }
  \item{positive}{ Logical. Does the result must be forced to positive values. }
  \item{\dots}{ Additional arguments. }
}
\details{
  This function computes one step forward prediction for a
  \code{\link{far}} model.

  Use the \code{newdata} option to input the past values,
  and the \code{label} option value to define the labels for the new
  observations. Notices that the output as the same length as
  \code{newdata} in the case of a FAR model, and the length of
  \code{newdata} minus one in the case of a FARX model. This is due to
  the time shift of the exogeneous variable: \eqn{X_{t+1}}{Xt+1} and
  \eqn{Y_{t}}{Yt} are used in the computation of
  \eqn{\hat{Y}_{t+1}}{Ypt+1}.

  In some special context, the user may need to suppress the
  \code{na.rm} observations with the \code{na.rm} option, or force the
  prediction to be positive with the \code{positive} option (in this
  case the result will be maximum of 0 and the predicted value).
}
\value{
  A \code{fdata} object.
}
\author{ J. Damon }

\seealso{
  \code{\link{far}}, \code{\link{pred.persist}},
  \code{\link{predict.kerfon}}. 
}

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

  # Cross validation (joined and separate)
  model1.cv <- far.cv(data=data1, y="X", x="Z", kn=8, ncv=10, cvcrit="X",
                center=FALSE, na.rm=FALSE, joined=TRUE)
  model2.cv <- far.cv(data=data1, y="X", x="Z", kn=c(4,4), ncv=10, cvcrit="X",
                center=FALSE, na.rm=FALSE, joined=FALSE)
  print(model1.cv)
  print(model2.cv)
  k1 <- model1.cv$minL2[1]
  k2 <- model2.cv$minL2[1:2]

  # Modelization of the FARX process (joined and separate)
  model1 <- far(data=data1, y="X", x="Z", kn=k1,
                center=FALSE, na.rm=FALSE, joined=TRUE)
  model2 <- far(data=data1, y="X", x="Z", kn=k2,
                center=FALSE, na.rm=FALSE, joined=FALSE)

  # Predicting values
  pred1 <- predict(model1,newdata=data1)
  pred2 <- predict(model2,newdata=data1)
  # Persistence
  persist1 <- pred.persist(select.fdata(data1,date=1:399),x="X")
  # Real values
  real1 <- select.fdata(data1,date=2:400)

  errors0 <- persist1[[1]]-real1[[1]]
  errors1 <- pred1[[1]]-real1[[1]]
  errors2 <- pred2[[1]]-real1[[1]]

  # Norm of observations
  summary(real1)
  # Persistence
  summary(as.fdata(errors0))
  # FARX models
  summary(as.fdata(errors1))
  summary(as.fdata(errors2))
}
\keyword{ ts }
\keyword{ methods }

