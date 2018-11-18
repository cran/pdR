\name{SMPLSplit_example}
\alias{SMPLSplit_example}

\title{Example code for sample splitting
}
\encoding{latin1}

\description{

A sample code for learning sample splitting.

}
\usage{
SMPLSplit_example(data,dep,indep,th1,th2,trim_per,rep,plot)
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

  \item{th1}{
  The first threshold variable
}
  \item{th2}{
  
  The second threshold variable

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

\details{This code is the learning example for learning Hansen's econometric sample splitting. I detailed the description of each threshold stage.  
}

\note{Original code offered by Dr. B. E.Hansen (http://www.ssc.wisc.edu/~bhansen/). 
}

\references{Hanse B. E. (2000) Sample Splitting and Threshold Estimation. Econometrica, 68, 575-603.
}

\examples{
## Not run, becasue of bootstrap replicaiton takes time. Users may unmark # and run. 
data("dur_john")
#rep <- 500
#trim_per <- 0.15
#dep <- 1
#indep <- c(2,3,4,5)
#th1 <- 6
#th2 <- 7
#OUT=SMPLSplit_example(data=dur_john,dep,indep,th1,th2,trim_per,rep,plot=1)
#OUT$TEST
#OUT$Hypothesis
#OUT$Threshold
#stat=matrix(as.numeric(OUT$TEST),byrow = TRUE,8,2)
#colnames(stat)=c("F-Stat","P-value")
#rownames(stat)=OUT$Hypothesis
#stat

}