\name{TICdata}
\alias{TICdata}
\docType{data}
\title{TIC data.
}
\description{A sample of 11 curves with TIC counts of a Liquid Chromatography-Mass Spectrometry (LS-MS) data set.
}
\usage{data(TICdata)}
\format{
  A data frame with 400 observations for each of 11 curves. The different rows
correspond to the different curves.
}
\source{
 Listgarten, J., Neal, R.M., Roweis, S.T. and Emili, A. (2005).
Multiple Alignment of Continuous Time Series, in 
Advances in Neural Information Processing Systems 17,
 Eds Saul, L.K., Weiss Y. and Bottou, L.,
 MIT Press, Cambridge, MA, 817--824.
}
\examples{
data(TICdata)
TIC=as.matrix(TICdata)

## Preparation of the TIC data for use in warping.

# for smoothing the LC-MS data TIC
library("SemiPar")

index = 1:200*2-1
TICy = t(matrix(index,200,11))
TIC = as.matrix(TICdata)
x = 1:400
for (i in 1:11)
	{
	TIC.sm = spm(TIC[i,]~f(x))
	TICy[i,] = TIC.sm$fit$fitted[index]
	}
TICx = t(matrix(index,200,11))
}
\keyword{datasets}
