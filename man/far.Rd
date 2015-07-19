\name{far}
\alias{far}
\alias{plot.far}
\alias{print.far}
\title{ FARX(1) model estimation }
\encoding{utf-8}
\description{
  Estimates the parameters of FAR(1) and FARX(1) processes
  (mean and autocorrelation operator)
}
\usage{
far(data, y, x, kn, center=TRUE, na.rm=TRUE, joined=FALSE)
}
\arguments{
  \item{data}{ A \code{fdata} object. }
  \item{y}{ A vector giving the name(s) of the endogenous variable(s) of
  the model. }
  \item{x}{ A vector giving the name(s) of the exogenous variable(s) of
  the model. }
  \item{kn}{ A vector giving the values of the various
  \code{kn} (dimension of plug-in in the algorithm). If it not supplied,
  the default value is one. }
  \item{center}{ Logical. Does the observation need to be centered. }
  \item{na.rm}{ Logical. Does the \code{n.a.} need to be removed. }
  \item{joined}{ Logical. If \code{TRUE}, the joined (whole) far model
  is computed, otherwise the model work with the separated variables. }
}
\details{
  \bold{The models}

  A Functional AutoRegressive of order 1 (\bold{FAR(1)}) process is, in a
  general way,  defined by the following equation:

  \deqn{T_{n}=\rho\left(T_{n-1}\right)+\epsilon_{n}, n \in
    Z}{Tn=r(Tn-1)+En, n in Z}

  where \eqn{T_{n}}{Tn} and \eqn{\epsilon_{n}}{En} take their values in a
  functional space (for instance an Hilbertian one), and \eqn{\rho}{r}
  is a linear operator. \eqn{\epsilon_{n}}{En} is a strong white noise.

  Now, let us consider a vector of observations, for instance:

  \deqn{\left(T_{1,n},...,T_{i,n},...,T_{m,n}\right)}{%
    (T1,n,...,Ti,n,...,Tm,n)}

  where each  \eqn{T_{i,n}}{Ti,n} lives in a one dimension functional
  space (not necessary the same). In the following, we will cut this
  list into two parts: the endogeneous variables \eqn{Y_{n}}{Yn} (the
  ones we are interested in), and the exogeneous variables
  \eqn{X_{n}}{Xn} (which influence the endogeneous ones).

  Then an order 1 Functional AutoRegressive process with eXogeneous
  variables (\bold{FARX(1)}) is defined by the equation:

  \deqn{Y_{n}=\rho\left(Y_{n-1}\right)+a\left(X_{n}\right)+\epsilon_{n}
    , n \in Z}{Yn=r(Yn-1)+a(Xn-1)+En, n in Z}

  where \eqn{\rho}{r} and \eqn{a} are linear operators in the
  adequate spaces.

  \bold{Estimation}

  This function estimates the parameters of FAR and FARX models.

  First, if the mean of the \code{data} is not zero (which is required
  by the model), you can substance this mean using the \code{center}
  option. Moreover, if the \code{data} contains \code{NA} values, you
  can work with it using the \code{na.rm} option.

  \bold{FAR Estimation}

  The estimation is mainly about estimating the \eqn{\rho}{r}
  operator. This estimation is done in a appropriate subspace (computed
  from the variance of the observations). What is important to know is
  that the best dimension \code{kn} for this subspace is not determined
  by this function. So the user have to supply this dimension using the
  \code{kn} option. A way to chose this dimension is to first use the
  \code{\link{far.cv}} function on the history.

  \bold{FARX Estimation}

  The FARX estimation can be realized by two methods: joined or not.

  The \bold{joined} estimation is done by ``joining'' the variables
  into one and estimating a FAR model on the resulting variable. For
  instance, with the previous notations, the transformation is:

  \deqn{T_{n}=\left(Y_{n},X_{n+1}\right)}{Tn=(Yn,Xn+1)}

  and \eqn{T_{n}}{Tn} is then a peculiar FAR(1) process. In such a case,
  you have to use the \code{joined=TRUE} oto the interpretation of
  this operatorption and specify \bold{one}
  value for \code{kn} (corresponding to the \eqn{T_{n}}{Tn} variable).

  Alternatively, you can choose not to estimate the FARX model by the
  joined procedure, then \code{kn} need to be a vector with a length
  equal to the number of variables involved in the FARX model
  (endogeneous and exogeneous).

  In both procedures, the endogeneous and exogeneous variables are
  provided through the \code{y} and \code{x} options respectively.


  \bold{Results}

  The function returns a \code{far} object. Use the \code{print},
  \code{coef} and \code{predict} functions to get more informations
  about the model.
}
\note{
  This function could be used to estimate FAR and FARX with order higher
  than 1 as a change of variables can transform the process to an
  order 1 FAR or FARX. For instance, if \eqn{T_{n}}{Tn} is a FAR(2)
  process then \eqn{Y_{n}=\left(T_{n},T_{n-1}\right)}{Yn=(Tn,Tn-1)} is a
  FAR(1) process.

  However, this is not a basic use of this function and may require a
  hard work of the user to get the result.
}
\value{
  A \code{far} object, see details for more informations.
}
\references{ Besse, P. and Cardot, H. (1996). Approximation spline de la
   prévision d'un processus fonctionnel autorégressif d'ordre 1.
   \emph{Revue Canadienne de Statistique/Canadian Journal of Statistics},
   \bold{24}, 467--487.

   Bosq, D. (2000) \emph{Linear Processes in Function Spaces: Theory and
   Applications}, (Lecture Notes in Statistics,
   Vol. 149). New York: Springer-Verlag.
}
\author{ J. Damon }

\seealso{ \code{\link{predict.far}}, \code{\link{far.cv}} }

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
\keyword{ models }
\keyword{ ts }
