\name{ptm}
\alias{ptm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Threshold Specification of panel data
%%  ~~function to do ... ~~ 
}
\description{panel threshold specification
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ptm(dep, ind1, ind2, d, bootn, trimn, qn, conf_lev, max_lag, t, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dep}{dependent variable
%%     ~~Describe \code{qn} here ~~
}
  \item{ind1}{regime dependent independent variables
%%     ~~Describe \code{ind1} here ~~
}
  \item{ind2}{ regime independent independent variables
%%     ~~Describe \code{ind2} here~~
}
  \item{d}{threshold variable
%%     ~~Describe \code{d} here~~
}
  \item{bootn}{vector of bootstrap repetition
%%     ~~Describe \code{bootn} here~~
}
  \item{trimn}{vector of trimmed percentage
%%     ~~Describe \code{trimn} here~~
}
  \item{qn}{number of quantiles to examine
%%     ~~Describe \code{qn} here~~
}
  \item{conf_lev}{confidence level
%%     ~~Describe \code{conf_lev} here~~
}
  \item{max_lag}{maximum number of lags
%%     ~~Describe \code{max_lag} here~~
}
  \item{t}{length of time period
%%     ~~Describe \code{t} here~~
}
  \item{n}{number of cross-section units
%%     ~~Describe \code{n} here~~
}
}
\details{This code fits only balanced panel data. It generalizes the simple code of Dr. Hansen (http://www.ssc.wisc.edu/~bhansen/), allowing multiple (more-than-one) regime-dependent (ind1) variables. We generalize the original code to better fit general need of threshold modeling in panel data. 
bootn and trimn are vector of 3 by 1, indicating numbers of three corresponding regimes. 
%%  ~~ If necessary, more details than the description above ~~
}

\references{Hanse B. E. (1999) Threshold effects in non-dynamic panels: Estimation, testing and inference. Journal of Econometrics,93, 345-368.
%% ~put references to the literature/web site here ~
}
\author{Prof. Tsung-wu Ho, Dept of Finance, Shih Hsin University, Taiwan, Taipei
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
#max_lag <- 1

#ptm(dep,ind1,ind2,d,bootn,trimn,qn,conf_lev,max_lag,t,n)

# In the output, the Regime-dependent Coefficients matrix
# is, from top to bottom, regime-wise.
}