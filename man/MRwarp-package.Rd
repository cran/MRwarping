\name{MRwarping-package}
\alias{MRwarping-package}
\alias{MRwarping}
\docType{package}
\title{
This package performs Bayesian multiresolution time warping for functional data.
}
\description{
Time warping is performed via a composition of warplets. The Bayesian model starts with one
warplet and adds warplets one at a time until the warping action becomes negligible in the sense
of having almost zero intensity or too narrow domains. The posterior distributions are used
as prior distributions for the extended model in the next step.
Warplets have an immediate interpretation as warping functions and the inverse warplet is trivial
to obtain. 
}
\details{
\tabular{ll}{
Package: \tab MRwarping\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-10-22\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}
}
\author{
L. Slaets, G. Claeskens, B.W. Silverman

Maintainer:  <Gerda.Claeskens@kuleuven.be>
}
\references{
Slaets, Claeskens and Silverman (2010). Warping functional data in R and C 
via a Bayesian Multiresolution approach. Journal of Statistical Software, 55(3), 1-22,\cr
URL http://www.jstatsoft.org/v55/i03/.


Claeskens, Silverman and Slaets (2010). A multiresolution approach to time warping 
achievec by a Bayesian prior-posterior transfer fitting strategy. Journal of the Royal
Statistical Society, Series B, 72(5), 673-694.
}
\keyword{ package }
