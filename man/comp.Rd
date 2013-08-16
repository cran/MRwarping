\name{comp}
\alias{comp}
\title{Constructs and evaluates a single warplet. 
}
\description{
A quartic warplet kernel is used to construct
a warplet with parameters a,lambda,r1,r2,x.
This function is used within MRwarp.
}
\usage{
comp(a, lambda, r1, r2, x)
}
\arguments{
  \item{a}{
Center of the warplet kernel.
}
  \item{lambda}{
Intensity of the warp, should be between (-1,1).
}
  \item{r1}{
Radius on the left-hand side of the center a, such that a-r1 
is the lower bound of the warping domain.
}
  \item{r2}{
Radius on the right-hand side of the center a, such that a+r2
is the upper bound of the warping domain.
}
  \item{x}{
Time point where to evaluate the warplet.
}
}

\value{
The vector of warped time points.
}
\references{
See the papers: Slaets, Claeskens and Silverman (2010). Warping functional data in R and C via a Bayesian Multiresolution approach.       Journal of Statistical Software, 55(3), 1-22,\cr
URL http://www.jstatsoft.org/v55/i03/.

Claeskens, Silverman and Slaets (2010). A multiresolution approach to time warping 
achievec by a Bayesian prior-posterior transfer fitting strategy. Journal of the Royal
Statistical Society, Series B, 72(5), 673-694.
}
\author{
L. Slaets, G. Claeskens, B.W. Silverman.
}
