\name{simul.far}
\alias{simul.far}
\title{ FAR(1) process simulation }
\description{
  Simulation of a FAR process using a Gram-Schmidt basis.
}
\usage{
simul.far(m=12,
  n=100,
  base=base.simul.far(24, 5),
  d.rho=diag(c(0.45, 0.9, 0.34, 0.45)),
  alpha=diag(c(0.5, 0.23, 0.018)),
  cst1=0.05)
}
\arguments{
  \item{m}{ Integer. Number of discretization points. }
  \item{n}{ Integer. Number of observations. }
  \item{base}{ A functional basis expressed as a matrix, as the matrix
    created by \code{\link{base.simul.far}} or with
    \code{\link{orthonormalization}}.  }
  \item{d.rho}{ Numerical matrix. Expression of the first bloc of the
    linear operator in the Gram-Schmidt basis. }
  \item{alpha}{ Numerical matrix. Expression of the first bloc of the
    covariance operator in the Gram-Schmidt basis. }
  \item{cst1}{ Numeric. Perturbation coefficient on the linear
    operator. }
}
\details{
  This function simulate a FAR(1) process with a strong white noise.

  The simulation is realized in two steps.

  First step, the function compute a FAR(1) process \eqn{T_n}{Tn} in a
  functional space (that we call in the sequel H) using a simple
  equation and the \code{d.rho}, \code{alpha} and \code{cst} parameters.

  Second step, the process \eqn{T_n}{Tn} is projected in the canonical
  basis using the \code{base} linear projector. 
  
  The \code{base} basis need to be a orthonormal basis wide enought. In the 
  contrary, the function use the \code{\link{orthonormalization}} function
  to make it so. Notice that the size of this matrix corresponds to the
  dimension of the "modelization space" H (let's call it
  \eqn{m_2}{m2}). Of course, the larger \code{m2} the better the
  functionnal approximation is. Whatever, keep in mind that
  \code{m2}=2\code{m} is a good compromise, in order to avoid the memory
  limits.

  In H, the linear operator \eqn{\rho}{rho} is expressed as:

  \deqn{
    \left(\begin{array}[c]{cc}%
    \code{d.rho} & 0 \cr%
    0 & eps.rho%
    \end{array}\right)%
  }{\code{d.rho} 0 \cr%
    0 eps.rho}

  Where \code{d.rho} is the matrix provided in the call, the two 0 are
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

  In H, \eqn{C^T}{CT}, the covariance operator of \eqn{T_n}{Tn}, is
  defined by:

    \deqn{
    \left(\begin{array}[c]{cc}%
    m_2 * \code{alpha} & 0 \cr%
    0 & eps.alpha%
    \end{array}\right)%
  }{m2 * \code{alpha} 0 \cr%
    0 eps.alpha}

  Where \code{alpha} is the matrix provided in the call, the two 0 are
  in fact two blocks of 0, and eps.alpha is a diagonal matrix having on
  his diagonal the terms:

  \deqn{\left(\epsilon_{k+1}, \epsilon_{k+2}, \ldots, %
    \epsilon_{\code{m2}}\right)}{(epsilon(k+1), epsilon(k+2), %
    \ldots, epsilon(\code{m2}))}

  where

  \deqn{\epsilon_{i}=\frac{\code{cst1}}{i}}{epsilon(i)=%
    \code{cst1}/i}
}
\value{
  A \code{fdata} object containing one variable ("var") which is a
  FAR(1) process of length \code{n} with \code{p} discretization
  points.
}
\note{
  To simulate \eqn{T_n}{Tn}, the function creates a white noise
  \eqn{E_n}{En} having the following covariance operator:

  \deqn{C^T - \rho * C^T * t({\rho})}{CT - rho * CT t(CT)}

  where t(.) is the transposition operator.
  \eqn{T_n}{Tn} is the computed using the equation:

  \deqn{T_{n+1} = \rho * T_n + E_n}{Tn+1 = rho * Tn + En}
}
\author{ J. Damon, S. Guillas }

\seealso{ \code{\link{simul.far.sde}}, \code{\link{simul.far.wiener}},
\code{\link{simul.farx}}, \code{\link{simul.wiener}},
\code{\link{base.simul.far}}.  }

\examples{
  far1 <- simul.far(m=64,n=100)
  summary(far1)
  print(far(far1,kn=4))
  par(mfrow=c(2,1))
  plot(far1,date=1)
  plot(select.fdata(far1,date=1:5),whole=TRUE,separator=TRUE)
}
\keyword{ ts }
\keyword{ misc }
