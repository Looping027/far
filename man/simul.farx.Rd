\name{simul.farx}
\alias{simul.farx}
\alias{theoretical.coef}
\title{ FARX(1) process simulation}
\description{
  Simulation of functional data with exogenous variables using a
  Gram-Schmidt basis.
	}
\usage{
simul.farx(m=12,n=100,base=base.simul.far(24,5),
  base.exo=base.simul.far(24,5),
  d.a=matrix(c(0.5,0),nrow=1,ncol=2),
  alpha.conj=matrix(c(0.2,0),nrow=1,ncol=2),
  d.rho=diag(c(0.45,0.90,0.34,0.45)),
  alpha=diag(c(0.5,0.23,0.018)),
  d.rho.exo=diag(c(0.45,0.90,0.34,0.45)),
  cst1=0.05)
theoretical.coef(m=12,base=base.simul.far(24,5),
  base.exo=NULL,
  d.rho=diag(c(0.45,0.90,0.34,0.45)),
  d.a=NULL,
  d.rho.exo=NULL,
  alpha=diag(c(0.5,0.23,0.018)),
  alpha.conj=NULL,
  cst1=0.05)
}
\arguments{
  \item{m}{ Integer. Number of discretization points. }
  \item{n}{ Integer. Number of observations. }
  \item{base}{ A functional basis expressed as a matrix, as the matrix
    created by \code{\link{base.simul.far}} or with
    \code{\link{orthonormalization}}.  }
  \item{base.exo}{ A functional basis expressed as a matrix, as the matrix
    created by \code{\link{base.simul.far}} or with
    \code{\link{orthonormalization}}.  }
  \item{d.rho}{ Numerical matrix. Part of the linear operator in the
    Gram-Schmidt basis (see details for more informations). }
  \item{d.a}{ Numerical matrix. Part of the linear operator in the
    Gram-Schmidt basis (see details for more informations). }
  \item{d.rho.exo}{ Numerical matrix. Part of the linear operator in the
    Gram-Schmidt basis (see details for more informations). }
  \item{alpha}{ Numerical matrix. Part of the linear operator in the
    Gram-Schmidt basis (see details for more informations). }
  \item{alpha.conj}{ Numerical matrix.  Part of the linear operator in the
    Gram-Schmidt basis (see details for more informations). }
  \item{cst1}{ Numeric. Perturbation coefficient on the linear
    operator. }
}
\details{
  The \code{simul.farx} function simulates a FARX(1) process with one
  endogeneous variable, one exogeneous variable and a strong white
  noise. To do so, the function uses the fact that a FARX(1) model can
  be seen as a FAR(1) model in a wider space. Therefore, the method is
  very similar to the one used by the function \code{\link{simul.far}}.

  The simulation is realized in two steps.

  First step, the function compute a FAR(1) process \eqn{T_n}{Tn} in a
  functional space (that we call in the sequel H) using a simple
  equation and the given parameters. \eqn{T_n}{Tn} is of the form
  \eqn{(T_{1n},T_{2n})}{(T1n,T2n)} where \eqn{T_{1n}}{T1n} and
  \eqn{T_{2n}}{T2n} are respectively the endogeneous and the exogeneous
  parts of the process.

  Second step, the process \eqn{T_n}{Tn} is projected in the canonical
  basis using the \code{base} and \code{base.exo} linear projectors to
  give the endogeneous (\eqn{X_n}{Xn}) and the exogeneous
  (\eqn{Z_n}{Zn}) variables respectively. 
  
  Those two basis need to be orthonormal and wide enought. In the
  contrary, the function use the \code{\link{orthonormalization}}
  function to make it so. Notice that the size of this matrix
  corresponds to the dimension of the "modelization space" H (let's call
  it \eqn{m_2=m1_2+m2_2}{m2=m12+m22}). Of course, the larger \code{m2}
  the better the functionnal approximation is. Whatever, keep in mind
  that \code{m2}=2\code{m} is a good compromise, in order to avoid the
  memory limits.

  In H, the linear operator \eqn{\rho}{rho} is expressed as:

  \deqn{
    \left(\begin{array}[c]{lc}%
    d.rho.mod & \code{d.a} \cr%
    0 & d.rho.exo.mod%
    \end{array}\right)%
  }{d.rho.mod \code{d.a} \cr%
    0 d.rho.exo.mod}

  Where d.rho.mod and d.rho.exo.mod are modified version of the provided
  \code{d.rho} and \code{d.rho.exo} respectively to avoid 0 on their
  diagonal. More precisely, the 0 on their diaginals are replaced by:

  \deqn{\left(\varepsilon_{k+1}, \varepsilon_{k+2}, \ldots, %
    \varepsilon_{\code{m2}}\right)}{(eps(k+1), eps(k+2), %
    \ldots, eps(\code{m2}))}

  where

  \deqn{\varepsilon_{i}=\frac{\code{cst1}}{i^2}+%
    \frac{1-\code{cst1}}{e^i}}{eps(i)=%
    \code{cst1}/(i^2)+(1-\code{cst1})/exp(i)}

  and k is the position in the \code{d.rho} or \code{d.r.ho.exo}
  diagonal.

  In H, \eqn{C^T}{CT}, the covariance operator of \eqn{T_n}{Tn}, is
  defined by:

    \deqn{
    \left(\begin{array}[c]{lc}%
    alpha.mod & alpha.conj.mod \cr%
    t(alpha.conj.mod) & alpha.exo%
    \end{array}\right)%
  }{alpha.mod alpha.conj.mod \cr%
    t(alpha.conj.mod) alpha.exo}

  Where alpha.mod and alpha.exo.mod are modified versions of
  \eqn{m1_2 * }{m12 * }\code{alpha} and \eqn{m2_2 * }{m22 *
  }\code{alpha.conj} respectively to avoid 0 on their diagonal. More
  precisely, the 0 on their diaginals are replaced by:

  \deqn{\left(\epsilon_{k+1}, \epsilon_{k+2}, \ldots, %
    \epsilon_{\code{m2b}}\right)}{(epsilon(k+1), epsilon(k+2), %
    \ldots, epsilon(\code{m2b}))}

  where

  \deqn{\epsilon_{i}=\frac{\code{cst1}}{i}}{epsilon(i)=%
    \code{cst1}/i}

  alpha.exo is a matrix representation of the covariance operator of
  \eqn{T_{2n}}{T2n} and is obtained by inverting the following relation:

  \deqn{alpha.conj.mod = d.rho.exo.mod * alpha.conj.mod * t(d.rho.mod) +%
    d.rho.exo.mod * mod.alpha * t(\code{d.a})}

  The \code{theoretical.coef} function is provided to help the user
  making comparison. Calling this function with the same parameters that
  where used in a simulation (realized with \code{simul.farx} or
  \code{simul.far}), we obtain the parameters used internaly by the
  function to make the simulation. Those values can therefore be
  compared to those obtained with the estimation function \code{far}
  (examples are provided below).
}
\value{
  A \code{fdata} object containing two variables ("X" the endogeous
  variable, and "Z" the exogeneous variable) which is a FARX(1) process
  of length \code{n} with \code{p} discretization points.
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
\code{\link{simul.far}}, \code{\link{simul.wiener}}.  }

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

  # Modelisation of the FARX process (joined and separate)
  model1 <- far(data1,k=4,joined=TRUE)
  model2 <- far(data1,k=c(3,1),joined=FALSE)

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
\keyword{ ts }
