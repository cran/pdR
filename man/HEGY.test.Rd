\name{HEGY.test}
\alias{HEGY.test}
\title{
Seasonal unit root test based on Hylleberg et al. (1990)
}
\encoding{latin1}

\description{
The function performs seasonal unit root test based on Eq.(3.6) of Hylleberg et al. (1990), univariate time series.
}
\usage{
HEGY.test(wts, itsd, regvar = 0, selectlags = list(mode = "signf", Pmax = NULL))
}
\arguments{
  \item{wts}{Univariate time series}
  \item{itsd}{Options for c(i,t,sd)\cr 
  i=1, intercept;=0 no intercept\cr
	t=1, trend;=0 no deterministic trend\cr
	sd=1, season dummy 1:(s-1);=0 no.
}
  \item{regvar}{Additional regressors}
  \item{selectlags}{Selection of lags\cr
mode, Criteria for selection, having three options: "signf","bic","aic". \cr
Pmax, maximum number of lags.
}
}
\details{
Mode for selectlags has three options, AIC and BIC use R built-in functions for linear model and their meanings are popular and straightforward. They include only lags that meet specific criterion, others are dropped from regressors. That is, lag orders of your model may not be a regular sequence. See also selPsignf() and selPabic().

}
\value{

\item{stats}{Tests statistics for HEGY regression coefficients.}
\item{hegycoefs}{HEGY regression coefficients.}
\item{lagsorder}{Lags order. "aic" or "bic" returns a scalar; "signf" returns a sequence of numbers}
\item{lagcoefs}{Coefficients of lag terms.}
\item{regvarcoefs}{Coefficient(s) of additional regressor(s).}

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
y<-inf_Q[,1]
hegy.out<-HEGY.test(wts=y, itsd=c(1,0,c(1:3)),regvar=0, selectlags=list(mode="aic", Pmax=12))

hegy.out$stats #HEGY test statistics
names(hegy.out) # HEGY objects, which can be called by using $, see below.
hegy.out$hegycoefs
hegy.out$regvarcoefs
}

\keyword{ Seasonal unit root test }
