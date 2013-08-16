\name{MRwarp}
\alias{MRwarp}
\title{
Main function to perform multiresolution warping for functional
data in a Bayesian way. 
}
\description{
The Bayesian procedure starts with one warplet in the model and
uses the posterior distributions as priors for a more extended model 
with one more warplet. The model is built with adding one warplet 
at a time and allows for amplitude variations.
}
\usage{
MRwarp(Xdata, Ydata, chain = 400, thin = 10, burnin = 200, kernel.s,
 	components = 1, selection = "FIXED", shr = 0.3, outputfit = 1, 
	alpha = 0.1)
}
\arguments{
  \item{Xdata}{
N by T matrix containing the T
x-coordinates or time points of the N curve observations.
Each row corresponds to a particular subject. No default.
}
  \item{Ydata}{
N by T matrix containing the T y-coordinates
or response values of the N curve observations. 
Each row corresponds to a particular subject. No default.
}
  \item{chain}{
The total number of MCMC iterations (default=400).
}
  \item{thin}{
The thinning factor of the MCMC algorithm (default=10).
}
  \item{burnin}{
The number of MCMC iterations that are discarded (default=200).}
  \item{kernel.s}{
Vector containing the starting values for the kernel parameters. No default.}
  \item{components}{
The number of warping components in the final model (default=1). The value is
ignored when selection="STEP".}
  \item{selection}{
Whether we want to estimate a  fixed number of warplets ("FIXED"),
or evaluate the warping procedure after each component ("STEP") (default=
"FIXED").}
  \item{shr}{
Determines the variance of the prior on the warplet intensities and shifts
(default = 0.3).}
  \item{outputfit}{
1 if the warped curves should be plotted after each estimated model, 0 otherwise (default=1).}
  \item{alpha}{
The significance level to be used in the model selection procedure (default=0.1).}
}
\value{
\item{last }{List of output values for the last fitted model.}
\item{previous }{List of output values for the one but last fitted model.}
\item{shift }{Component of the output list. 
A MCMC chain of the estimated horizontal shift for curves 1
to \eqn{N}.}
\item{warping}{Component of the output list that is itself a list 
containing four quantities: lower, A, upper and Intensities.}
\item{lower }{Component of warping. Adaptive MCMC chains of the estimated warping lower bounds 
\eqn{(w_{l,1}, \ldots, w_{l,Q-1})}.} 
\item{A }{Component of warping.
Adaptive MCMC chains of the estimated warping centers  \eqn{(a_1, \ldots, a_{Q-1})}.} 
\item{upper }{Component of warping.
 AMCMC chains of the estimated warping upper bounds  
\eqn{(w_{u,1}, \ldots, w_{u,Q-1})}.}
\item{Intensities }{Component of warping. 
Adaptive MCMC chains of the estimated warping intensities, first all \eqn{N}
values for warplet 1,\ldots, all \eqn{N} values for warplet \eqn{Q}.
}
\item{kernels }{Component of the output list. 
Adaptive MCMC chains of the estimated kernel lower 
bounds, centers and upper bounds.}
\item{error.variance }{Component of the output list.
The estimated value of the error variance.} 
\item{max.post.dens }{Component of the output list. 
The row in the parameter chain vectors/matrices
corresponding to the highest posterior pseudo-log-likelihood.} 
}
\references{
See the papers: Slaets, Claeskens and Silverman (2010). Warping functional data in R and C 
via a Bayesian Multiresolution approach. Journal of Statistical Software, 55(3), 1-22,\cr
URL http://www.jstatsoft.org/v55/i03/.


Claeskens, Silverman and Slaets (2010). A multiresolution approach to time warping 
achievec by a Bayesian prior-posterior transfer fitting strategy. Journal of the Royal
Statistical Society, Series B, {72}(5), 673-694.
}
\author{
L. Slaets, G. Claeskens, B.W. Silverman.
}


\examples{
data(TICdata)
TIC = as.matrix(TICdata)

index = 1:200*2-1
TICy = t(matrix(index,200,11))
x = 1:400
for (i in 1:11)
	{
	TIC.sm = spm(TIC[i,]~f(x))
	TICy[i,] = TIC.sm$fit$fitted[index]
	}
TICx = t(matrix(index,200,11))
kernel.s = c(70,100,130,270,285,300)

\dontrun{
output = MRwarp(Xdata=TICx,Ydata=TICy,chain=100,thin=5,burnin=50,kernel.s,
components=1,selection="FIXED",shr=0.3,outputfit=1,alpha=0.1)
}
}
