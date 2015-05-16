\name{pIGF}
\alias{pIGF}

\title{Panel unit root test of Chang(2002)}

\encoding{latin1}


\description{Compute the panel unit root test statistic of Chang(2002).

}


\usage{
pIGF(datamat, maxp, ic, spec)
}

\arguments{
  \item{datamat}{T by N panel data.T is the time length,N is the number of cross-section units.

}
  \item{maxp}{the max number of lags

}
  \item{ic}{Information criteria, either "AIC" or "BIC".

}
  \item{spec}{model specification.\cr
=0, no intercept and trend.\cr
=1, intercept only.\cr
=2, intercept and trend.

}
}

\details{This function estimates the panel unit root test based on univariate instrument generating function of (Chang,2002).

}
\value{
\item{ panel.tstat }{ panel IGF test statistics }
\item{ pvalue }{ P-value of the panel.tstat }

}

\references{
Chang,Y.(2002) Nonlinear IV Unit Root Tests in Panels with Cross-Sectional Dependency. Journal of Econometrics, 110, 261-292.
}

\author{Ho Tsungwu <tsungwu@cc.shu.edu.tw>}

\examples{
data(inf19)
datam <- inf19
pIGF(datam,maxp=25,ic="BIC",spec=2)
}