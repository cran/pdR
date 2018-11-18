\name{SMPLSplit_est}
\alias{SMPLSplit_est}

\title{Estimation of sub-sampled data

}

\encoding{latin1}

\description{

A function for estimating the subsampled data.

}
\usage{
SMPLSplit_est(dat,dep,indep,th,plot)
}
\arguments{
  \item{dat}{
  The data in either data.frame or matrix  

}

  \item{dep}{
  The number of culumn of dependent variable

}
  \item{indep}{
  The number of culumns of  regime dependent independent variables: 
}

  \item{th}{

  The threshold variable

}

  \item{plot}{
  =1, plot; =0, do not plot
}

}
\details{This code estimates the parameters of sub-sampled data. It generalizes the simple code of Dr. Hansen, allowing Heteroskedastic Errors (White Corrected).
}

\note{Original code offered by Dr. B. E.Hansen (http://www.ssc.wisc.edu/~bhansen/). 
}

\references{Hanse B. E. (2000) Sample Splitting and Threshold Estimation. Econometrica, 68, 575-603.
}
