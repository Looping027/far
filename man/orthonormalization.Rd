\name{orthonormalization}
\alias{orthonormalization}
\title{ Orthonormalization of a set of a matrix }
\description{
  Gram-Schmidt orthogonalization of a matrix considering
  its columns as vectors. Normalization is provided as will.
}
\usage{
orthonormalization(u, basis=TRUE, norm=TRUE)
}
\arguments{
  \item{u}{ a matrix (n x p) representing n different vectors 
  in a n dimensional space  }
  \item{basis}{ does the returned matrix have to be a basis }
  \item{norm}{ does the returned vectors have to be normed }
}
\details{
  This is a simple application of the Gram-Schmidt algorithm of 
  orthogonalization (please note that this process was
  presented first by Laplace).
  
  The user provides a set of vector (structured in a matrix)
  and the function calculate a orthogonal basis of the same
  space. If desired, the returned basis can be normed,
  or/and completed to cover the hole space.
  
  If the number of vectors in \code{u} is greater than the
  dimension of the space (that is if n > p), only the first
  p columns are taken into account to computed the result.
  A warning is also provided.
  
  The only assumption made on \code{u} is that the span space
  is of size min(n,p). In other words, there must be no
  colinearities in the initial set of vector.
}
\value{
  The orthogonalized matrix obtained from \code{u} where the
  vector are arranged in columns.
  
  If \code{basis} is set to \code{TRUE}, the returned matrix
  is squared.
}
\author{ J. Damon }

\examples{
  mat1 <- matrix(c(1,0,1,1,1,0),nrow=3,ncol=2)
  orth1 <- orthonormalization(mat1, basis=FALSE, norm=FALSE)
  orth2 <- orthonormalization(mat1, basis=FALSE, norm=TRUE)
  orth3 <- orthonormalization(mat1, basis=TRUE, norm=TRUE)
  crossprod(orth1)
  crossprod(orth2)
  crossprod(orth3)
}
\keyword{ algebra }

