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
  the data in either data.frame or matrix  

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
  \item{trim_per}{
  trimmed percentage.

}
  \item{rep}{
  number of bootstrap repetition.

}
  \item{plot}{
  =1, plot; =0, do not plot.
}

}

\value{

\item{fstat}{LM-test for no threshold.}
\item{pvalue}{bootstrap P-Value.}

}




\details{This code tests for the presence of threshold. It generalizes the simple code of Dr. Hansen, allowing Heteroskedastic Errors (White Corrected). }

\note{Original code offered by Dr. B. E.Hansen (http://www.ssc.wisc.edu/~bhansen/). 
}

\references{Hanse B. E. (2000) Sample Splitting and Threshold Estimation. Econometrica, 68, 575-603.
}
