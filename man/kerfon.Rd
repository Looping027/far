\name{kerfon}
\alias{kerfon}
\alias{print.kerfon}
\title{ Functional Kernel estimation }
\description{
  Modelization of \code{fdata} using functional kernel.
}
\usage{
kerfon(data, x, r, hmin, hmax, na.rm=TRUE)
}
\arguments{
  \item{data}{ A \code{fdata} object. }
  \item{x}{ The name of the studied variable. }
  \item{r}{ Number of observations used to cross validate the model. }
  \item{hmin}{ Minimal value of the bandwidth. }
  \item{hmax}{ Maximal value of the bandwidth. }
  \item{na.rm}{ Logical. Does the \code{n.a.} need to be removed. }
}
\details{
  This function constructs a functional kernel model and performs the
  estimation of it's bandwidth.

  One nonparametric way to deal with the conditional expectation
  \eqn{\rho(x)=\mbox{I$\!$E}\left[  X_{i}\left|  X_{i-1}=x\right.  \right]}{
    r(x)=E[X(i)|X(i-1)=x]}, where \eqn{\left(X_{i}\right)}{X(i)} is a
  $H$-valued process, is to consider a predictor inspired by the
  classical kernel regression, as in Nadaraja and Watson. This estimator
  is defined by :

  \deqn{\hat{\rho}_{h_{n}}(x)=\frac{\sum\limits_{i=1}^{n-1}X_{i+1}\cdot%
      K\left(\frac{\left\|  X_{i}-x\right\|  _{H}}{h_{n}}\right)%
    }{\sum\limits_{i=1}^{n-1}K\left(  \frac{\left\|  X_{i}-x\right\|%
	_{H}}{h_{n}}\right)},x\in H}{r*(x,hn)=sum%
    of(X(i+1)*K(||X(i)-x||/hn)))/((n-1) * sum%
    of(K(||X(i)-x||/hn))))}

  Where K is a kernel, \eqn{\left\|.\right\|  _{H}}{||.||} is the
  norm in H, and \eqn{h_n}{hn} is the bandwidth (\eqn{\in
    \mbox{I$\!$R}^+_*}{in R+*}).
  
  The function \code{kerfon} use the cross validation to determinate a
  value for \eqn{h_n}{hn}. This method have been chosen because of the
  lack of theoretical results about this model. The parameters
  \code{hmin} and \code{hmax} are used, when provided, to control the
  permissible values of \eqn{h_n}{hn}. By default, those parameters are
  respectively equals to \eqn{\sigma/8} and \eqn{4*\sigma}, where
  \eqn{\sigma} is the estimated squared root of the variance operator of
  X. To choose the value of \eqn{h_n}{hn}, you need to provide the same
  value for both \code{hmin} and \code{hmax}.

  During the cross-validation, considering that the fdata object
  \code{x} contains \eqn{n} observations, the function use the first
  \eqn{(n-r)} observations as the past values, and compute the mean
  square norm of the errors on the last \eqn{r} observations.

  Of course, if the model created is then used to compute prediction
  through \code{\link{predict.kerfon}}, the whole set of observations (the
  \eqn{n} observations) are used as the past values.

  As \code{fdata} object may contains several variables, a way is
  provided to select the studied variable (the function only works
  with one variable for the moment).
}
\value{
  A kerfon object. A method for the \code{\link{print}} function is
  provided.

  For information, the object is a list with the following elements :

  \item{call}{the call of the function.}
  \item{h}{the bandwidth (three values : optimal, minimum, maximum)}
  \item{x}{the name of the chosen variable}
  \item{xdata}{the past values for \code{x}}
  \item{ydata}{the associated values for \code{x}}
}
\author{ J. Damon }

\seealso{ \code{\link{predict.kerfon}} }

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
}
\keyword{ ts }
\keyword{ nonlinear }
