\name{SMPLSplit_est}
\alias{SMPLSplit_est}

\title{Estimation of sub-sampled data

}

\encoding{latin1}

\description{

A function for estimating the subsampled data.

}
\usage{
SMPLSplit_est(data,dep,indep,th,plot,h=1,nonpar=2)
}
\arguments{
  \item{data}{
  the data in either data.frame or matrix.  

}

  \item{dep}{
  the name of dependent variable.

}
  \item{indep}{
  the name(s) of independent variable(s). 
}

  \item{th}{

  the name of threshold variable.

}

  \item{plot}{
  =1, plot; =0, do not plot.
}

\item{h}{
h=1, heteroskedasticity-consistent covariance; h=0, homoskedastic case. 
}

\item{nonpar}{
  Indicator for non-parametric method used to estimate nuisance scale in the presence of
  heteroskedasticity (only relevant if h=1).Set nonpar=1 to estimate regressions using
  a quadratic.Set nonpar=2 to estimate regressions using an Epanechnikov kernel with automatic bandwidth.
}

}


\value{

\item{threshold}{values of threshold estimates.}
\item{est0}{coefficient estimates of global data.}
\item{est.low}{coefficient estimates of low regime.}
\item{est.high}{coefficient estimates of high regime.}
\item{est0.info}{additional information of global data.}
\item{est.joint.info}{additional information of joint threshods.}
\item{est.low.info}{additional information of est.low.}
\item{est.high.info}{additional information of est.high.}
}


\details{This code estimates the parameters of sub-sampled data. It generalizes the simple code of Dr. Hansen, allowing White Corrected Heteroskedastic Errors.
}

\note{Original code offered by Dr. B. E.Hansen (http://www.ssc.wisc.edu/~bhansen/). 
}

\references{Hanse B. E. (2000) Sample Splitting and Threshold Estimation. Econometrica, 68, 575-603.
}

\examples{
## Not run, becasue of bootstrap replicaiton takes time. Users may unmark # and run. 
data("dur_john")
rep <- 500
trim_per <- 0.15
dep <- "gdpGrowth"
indep <- colnames(dur_john)[c(2,3,4,5)]

SMPLSplit_est(data=dur_john,dep,indep,th="GDP60",plot=0,h=1,nonpar=2)


}
