\name{SMPLSplit_het}
\alias{SMPLSplit_het}

\title{Testing for sample splitting

}

\encoding{latin1}

\description{

A function for testing sample split given subsampled data.

}
\usage{
SMPLSplit_het(data,dep,indep,th,trim_per,rep,plot)
}
\arguments{
  \item{data}{
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
  \item{trim_per}{
  Trimmed percentage

}
  \item{rep}{
  Number of bootstrap repetition

}
  \item{plot}{
  =1, plot; =0, do not plot
}

}
\details{This code tests for the presence of threshold. It generalizes the simple code of Dr. Hansen, allowing Heteroskedastic Errors (White Corrected). }

\note{Original code offered by Dr. B. E.Hansen (http://www.ssc.wisc.edu/~bhansen/). 
}

\references{Hanse B. E. (2000) Sample Splitting and Threshold Estimation. Econometrica, 68, 575-603.
}
