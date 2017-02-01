\name{ptm}
\alias{ptm}

\title{Threshold Specification of panel data

}
\description{

A generalized specification for estimating panel threshold model.

}
\usage{
ptm(dep, ind1, ind2, d, bootn, trimn, qn, conf_lev, t, n)
}
\arguments{
  \item{dep}{
  Dependent variable

}
  \item{ind1}{
  Independent variables: regime dependent 
}
  \item{ind2}{
 Independent variables:regime independent

}
  \item{d}{
  Threshold variable
}
  \item{bootn}{
  
  Vector of bootstrap repetition

}
  \item{trimn}{
  Vector of trimmed percentage

}
  \item{qn}{
  Number of quantiles to examine

}
  \item{conf_lev}{
  Confidence level

}

  \item{t}{
  Length of time period

}
  \item{n}{
  Number of cross-section units

}
}
\details{This code fits only balanced panel data. It generalizes the simple code of Dr. Hansen (http://www.ssc.wisc.edu/~bhansen/), allowing multiple (more-than-one) regime-dependent (ind1) variables. We generalize the original code to better fit general need of threshold modeling in panel data. \cr
bootn and trimn are vector of 3 by 1, indicating numbers of three corresponding regimes.\cr
This version corrects a slight error incurred by argument max_lag, which is used by Hansen to arrange investment data via lags. In this package, users manipulate data to fit personal research to ptm(), hence this argument is omitted lest should degree of freedom will suffer a loss of N.
}

\references{Hanse B. E. (1999) Threshold effects in non-dynamic panels: Estimation, testing and inference. Journal of Econometrics,93, 345-368.
}
\author{Ho Tsungwu <tsungwu@cc.shu.edu.tw>
}


\examples{
# library(pdR)
#data(invest)
#dat<-invest[1:1500,]    # subsetting the first 1500 obs., #for simplicity
#t <- 15            #Length of time period
#nt <- nrow(dat)
#n <- nt/t           # number of cross-section units

#dep<- as.matrix(dat[,1])       # investment/assets
#th1<- as.matrix(dat[,2])  #Tobin's Q
#th2<- as.matrix(dat[,3])  # cash-flow/assets
#ind1<- cbind(th1,th2)  #regime-dep covariates 
#d <- as.matrix(dat[,4])      # Threshold variable      
#ind2 <- cbind((th1^2),(th1^3),(th1*d)) # regime-indep covariates:
#bootn<-c(100,200,300)  # bootstrapping replications for each  threshold esitmation
#trimn<-c(0.05,0.05,0.05)   #trimmed percentage for each threshold  esitmation

#qn<-400
#conf_lev<-0.95

#Output=ptm(dep,ind1,ind2,d,bootn,trimn,qn,conf_lev,t,n)
#Output[[1]] #Formatted output of 1st threshold, 2 regimes
#Output[[2]] #Formatted output of 2nd threshold, 3 regimes
#Output[[3]] #Formatted output of 3rd threshold, 4 regimes

# In the output, the Regime-dependent Coefficients matrix
# is, from top to bottom, regime-wise.
}