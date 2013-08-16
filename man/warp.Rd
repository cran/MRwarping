\name{warp}
\alias{warp}
\title{Evaluates a composition of warplets.
}
\description{The function warp evaluates a composition of warplets that 
are constructed by the function comp.
}
\usage{
warp(A, Lambda, R1, R2, x)
}
\arguments{
  \item{A}{Vector of centers of the warplets.
}
  \item{Lambda}{Vector of intensities.
}
  \item{R1}{Vector of radii on the left-hand side of the centers.
}
  \item{R2}{Vector of radii on the right-hand side of the centers.
}
  \item{x}{Vector of time points at which to evaluate the warping function.
}
}
\value{warp contains the warping function evaluated at x, see also comp.
}
\references{
See the papers: Slaets, Claeskens and Silverman (2010). Warping functional data in R and C 
via a Bayesian Multiresolution approach. Journal of Statistical Software, 55(3), 1-22,\cr
URL http://www.jstatsoft.org/v55/i03/.


Claeskens, Silverman and Slaets (2010). A multiresolution approach to time warping 
achievec by a Bayesian prior-posterior transfer fitting strategy. Journal of the Royal
Statistical Society, Series B, 72(5), 673-694.
}
\author{L. Slaets, G. Claeskens, B.W. Silverman.
}

\examples{
t = seq(0,10,length.out=1000)
tau.t = warp(c(5,2),c(0.6,0.4),c(2,1.5),c(3,2),t)

## The function is currently defined as
warp =
function(A,Lambda,R1,R2,x)
{
Wx = x
for (i in 1:length(A))
    	{
	warp = comp(A[i],Lambda[i],R1[i],R2[i],Wx)
	Wx = warp
	}
return(Wx)
}
}


