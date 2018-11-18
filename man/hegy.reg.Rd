\name{hegy.reg}
\alias{hegy.reg}
\title{
Generate the HEGY regressors.
}
\encoding{latin1}

\description{
This function generates the level regresors in HEGY regression, without differenced lag terms.
}
\usage{
hegy.reg(wts)
}

\arguments{
  \item{wts}{
Univariate time series, with a possibly seasonal stochastic trend
}
}
\details{
This function automatically identifies the frequency of time series data, and generate necessary level components as described in  Eq.(3.7) of Hylleberg et. al (1990).
}

\references{
Hylleberg, S., Engle, R.F., Granger, C.W.J., and Yoo, B.S.(1990) Seasonal integration and cointegration. Journal of Econometrics,44, 215-238.\cr
Javier Lopez-de-Lacalle in https://github.com/cran/uroot/blob/master/R/hegy.R
}
\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, modifed from Javier Lopez-de-Lacalle
}

\examples{
data(inf_Q)
y=inf_Q[,1]
hegy.reg(y)
}

