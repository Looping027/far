\name{simul.far.wiener}
\alias{simul.far.wiener}
\title{ FAR(1) process simulation with Wiener noise }
\encoding{utf-8}
\description{
  Simulation of a FAR(1) process using a Wiener noise.
}
\usage{
simul.far.wiener(m=64, n=128,
d.rho=diag(c(0.45, 0.9, 0.34, 0.45)), cst1=0.05, m2=NULL)
}
\arguments{
  \item{m}{ Integer. Number of discretization points. }
  \item{n}{ Integer. Number of observations. }
  \item{d.rho}{ Numerical matrix. Expression of the first bloc of the
    linear operator in the Karhunen-Loève basis. }
  \item{cst1}{ Numeric. Perturbation coefficient on the linear
    operator. }
  \item{m2}{ Integer. Length of the Karhunen-Loève expansion (2\code{m}
    by default). }
}
\details{
  This function simulate a FAR(1) process with a Wiener noise. As for
  the \code{\link{simul.wiener}}, the function use the Karhunen-Loève
  expansion of the noise. The FAR(1) process, defined by its linear
  operator (see \code{\link{far}} for more details), is computed in the
  Karhunen-Loève basis then projected in the natural basis. The
  parameters given in input (\code{d.rho} and \code{cst1}) are expressed
  in the Karhunen-Loève basis.

  The linear operator, expressed in the Karhunen-Loève basis, is of the
  form:

  \deqn{
    \left(\begin{array}[c]{lc}%
    \code{d.rho} & 0 \cr%
    0 & eps.rho
    \end{array}\right)
  }{\code{d.rho} 0 \cr
    0 eps.rho}

  Where \code{d.rho} is the matrix provided in ths call, the two 0 are
  in fact two blocks of 0, and eps.rho is a diagonal matrix having on
  his diagonal the terms:

  \deqn{\left(\varepsilon_{k+1}, \varepsilon_{k+2}, \ldots, %
    \varepsilon_{\code{m2}}\right)}{(eps(k+1), eps(k+2), %
    \ldots, eps(\code{m2}))}

  where

  \deqn{\varepsilon_{i}=\frac{\code{cst1}}{i^2}+%
    \frac{1-\code{cst1}}{e^i}}{eps(i)=%
    \code{cst1}/(i^2)+(1-\code{cst1})/exp(i)}

  and k is the length of the \code{d.rho} diagonal.

  The \code{d.rho} matrix can be viewed as the information and the
  eps.rho matrix as a perturbation. In this logic, the norm of eps.rho
  need to be smaller than the one of \code{d.rho}.
}
\value{
  A \code{fdata} object containing one variable ("var") which is a
  FAR(1) process of length \code{n} with \code{m} discretization
  points.
}
\references{ Pumo, B. (1992). \emph{Estimation et Prévision de
  Processus Autoregressifs Fonctionnels. Applications aux
  Processus à Temps Continu.}
  PhD Thesis, University Paris 6, Pierre et Marie Curie. }
\author{ J. Damon }

\seealso{ \code{\link{fdata}}, \code{\link{far}} ,
    \code{\link{simul.far.wiener}}.  }

\examples{
  far1 <- simul.far.wiener(m=64,n=100)
  summary(far1)
  print(far(far1,kn=4))
  par(mfrow=c(2,1))
  plot(far1,date=1)
  plot(select.fdata(far1,date=1:5),whole=TRUE,separator=TRUE)
}
\keyword{ ts }
\keyword{ misc }
