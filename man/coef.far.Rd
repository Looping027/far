\name{coef.far}
\alias{coef.far}
\title{ Extract Model Coefficients }
\description{
  'coef' method to extract the linear operator of a FAR model.
}
\usage{
  \method{coef}{far}(object, ...)
}
\arguments{
  \item{object}{ An object of type \code{far}. }
  \item{\dots}{ Other arguments (not used in this case). }
}
\details{
  Give the matricial representation of the linear operator express in
  the canonical basis. See \code{\link{far}} for more details about the
  meaning of this operator.

  If the \code{far} model is used on a one dimensional variable or with
  the \code{joined=TRUE} option, then the matrix has a dimension equal
  to the subspace dimension.

  In the other case, the dimension of the matrix is equal to the sum of
  the dimensions of the various subspaces. In such a case, the order of
  the variables in the matrix is the same as in the vector
  \code{c(y,x)}. For instance, if \code{kn=c(3,2)} with \code{y="Var1"}
  and \code{x="Var3"} then:
  \itemize{
    \item The first 3x3 first bloc of the matrix is the autocorrelation
    of ``Var1''.
    \item The 3x2 up right bloc of the matrix is the correlation of
    ``Var3'' on ``Var1''.
    \item The 2x3 down left bloc of the matrix is the correlation of
    ``Var1'' on ``Var3''.
    \item The 2x2 down right bloc of the matrix is the autocorrelation of
    ``Var3''.
  }
}
\value{
  A square matrix of size (raw and column) equal to the sum of the
  element of \code{kn}.
}
\author{ J. Damon, S. Guillas }

\seealso{ \code{\link{far}},\code{\link{coef}} }

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

  # Modelization of the FARX process (joined and separate)
  model1 <- far(data1,kn=4,joined=TRUE)
  model2 <- far(data1,kn=c(3,1),joined=FALSE)

  # Calculation of the theoretical coefficients
  coef.theo <- theoretical.coef(m=10,base=base.simul.far(20,5),
                base.exo=base.simul.far(20,5),
                d.a=matrix(c(0.5,0),nrow=1,ncol=2),
                alpha.conj=matrix(c(0.2,0),nrow=1,ncol=2),
                d.rho=diag(c(0.45,0.90,0.34,0.45)),
                alpha=diag(c(0.5,0.23,0.018)),
                d.rho.exo=diag(c(0.45,0.90,0.34,0.45)),
                cst1=0.0)

  # Joined coefficient
  round(coef(model1),2)
  coef.theo$rho.T

  # Separate coefficient
  round(coef(model2),2)
  coef.theo$rho.X.Z
}
\keyword{ algebra }
