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
  the data in either data.frame or matrix. 

}

  \item{dep}{
  the name of dependent variable.

}
  \item{indep}{
  the name(s) of independent variable(s) 
}

  \item{th1}{
  the first threshold variable.
}
  \item{th2}{
  
  the second threshold variable.

}
  \item{trim_per}{
  trimmed percentage.

}
  \item{rep}{
  nNumber of bootstrap repetitions.

}
  \item{plot}{
  =1, plot; =0, do not plot.
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
#dep <- "gdpGrowth"
#indep <- colnames(dur_john)[c(2,3,4,5)]
#th1 <- "GDP60"
#th2 <- "Literacy"
#OUT=SMPLSplit_est(data=dur_john,dep,indep,th=th1,plot=0,h=1,nonpar=2)
#OUT$TEST
#OUT$Hypothesis
#OUT$Threshold
#stat=matrix(as.numeric(OUT$TEST),byrow = TRUE,8,2)
#colnames(stat)=c("F-Stat","P-value")
#rownames(stat)=OUT$Hypothesis
#stat

}