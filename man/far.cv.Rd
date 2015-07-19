\name{far.cv}
\alias{far.cv}
\title{ Cross Validation for FARX(1) model }
\description{
  Cross Validation for FAR(1) and FARX(1) models
}
\usage{
far.cv(data, y, x, kn, ncv, cvcrit, center=TRUE, na.rm=TRUE, joined=FALSE)
}
\arguments{
  \item{data}{ A \code{fdata} object. }
  \item{y}{ A vector giving the name(s) of the endogenous variable(s) of
  the model. }
  \item{x}{ A vector giving the name(s) of the exogenous variable(s) of
  the model. }
  \item{kn}{ A vector giving the maximum values of the various
  \code{kn} (dimension of plug-in in the algorithm). If it not supplied,
  the number of discretization point is used. }
  \item{ncv}{ Number of observations used to the cross validation }
  \item{cvcrit}{ A vector of characters. Name of the variable used to
    measure the errors (\code{y} by default). }
  \item{center}{ Logical. Does the observation need to be centered. }
  \item{na.rm}{ Logical. Does the \code{n.a.} need to be removed. }
  \item{joined}{ Logical. If \code{TRUE}, the joined (whole) far model
  is computed, otherwise the model work with the separated variables. }
}
\details{
  In order to perform good forecasting with a FAR or FARX model, you need
  to determine the dimensions \code{kn} of the subspace in which the linear
  operator is estimated (see \code{\link{far}} for more details).

  This function helps the user to do this choice by performing a cross
  validation on a test sample. The usage is close of the
  \code{\link{far}} function, so we will discuss about the options which
  differ.

  First, the \code{kn} option is used to restrict the values searched:
  this is a vector containing the maxima values. As in
  \code{\link{far}}, the dimension of this vector is function of the
  number of variables involved in the model and the type of estimation
  done (joined or not).

  \code{ncv} is the number of observation used to test the models. If it
  is not provided, the function use the last fifth of the observations in
  \code{data}. In such a case, the four first fifth are used to
  estimates the models. This is in general a good compromise.

  Finally, \code{cvcrit} list the variables used to test the models. If
  more than one variable is provided, the test is calculated as a mean
  of the errors over all the variables.

  The criteria used to test the (functional) errors are the norms L1,
  L2, L infinite, L1 on the maxima, L2 on the maxima, and L infinite on the
  maxima.
}
\value{
  It is a LIST with the following elements
  \item{cv }{Matrix giving the various errors (L1, L2, L infinite, L1 on the
  maxima, L2 on the maxima, L infinite on the maxima) for the tested
  values of \code{kn}}
  \item{minL1 }{A vector corresponding to the row of \code{cv} where the
  L1 error minima is obtained}
  \item{minL2 }{A vector corresponding to the row of \code{cv} where the
  L2 error minima is obtained}
  \item{minLinf }{A vector corresponding to the row of \code{cv} where the
  L infinite error minima is obtained}
  \item{minL1max }{A vector corresponding to the row of \code{cv} where the
  L1 maxima's error minima is obtained}
  \item{minL2max }{A vector corresponding to the row of \code{cv} where the
  L2 maxima's error minima is obtained}
  \item{minLinfmax }{A vector corresponding to the row of \code{cv} where the
  L infinite maxima's error minima is obtained}
}

\author{ J. Damon }

\seealso{ \code{\link{far}}, \code{\link{fdata}} }

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
  print(model1)
  print(model2)
}
\keyword{ ts }
\keyword{ models }
