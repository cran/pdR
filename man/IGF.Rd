\name{IGF}

\alias{IGF}

\title{Unit root test based on Change(2002)
}

\encoding{latin1}

\description{This function estimates the unit root regression based on instrument generating function of Change(2002) and returns useful outputs.
}

\usage{
IGF(y, maxp, ic, spec)
}

\arguments{
  \item{y}{A univariate time series data

}
  \item{maxp}{the max number of lags

}
  \item{ic}{Information criteria, either "AIC" or "BIC"

}
  \item{spec}{regression model specification.\cr
=0, no intercept and trend.\cr
=1, intercept only.\cr
=2, intercept and trend.

}
}
\details{Estimate univariate unit root test of Chang(2002).}

\value{
\item{tstat.IGF}{IGF unit root test}
\item{beta}{regression coefficients. The first one is the AR(1) coefficient of unit root, and the last one is the intercept or trend}
\item{sdev}{The IGF standard error for unit root coefficient}
\item{cV}{The scalar C in IGF equation}
\item{p}{The optimal number of lag}

}
\references{
Chang, Y. (2002) Nonlinear IV Unit Root Tests in Panels with Cross-Sectional Dependency. Journal of Econometrics, 110, 261-292.

}
\author{
Ho Tsungwu <tsungwu@cc.shu.edu.tw>
}

\examples{
data(inf19)
y <- inf19[,1]
IGF(y,maxp=35,ic="BIC",spec=2)$tstat.IGF

}

\keyword{ instrument generating functions }
\keyword{ unit root }