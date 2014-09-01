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
%%  ~~who you are~~
}
\note{In addition to Dr. Hanse's contributions to threshold modeling, in the future, this package plans to include some special panel unit root tests based upon instrument generating functions(IGF), proposed by Chang (2002, J. of Econometrics), and seasonal panel unit root test, and regression quantile panel linear model(rqplm).. 
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# not run. This example is merely an illustration, to execute it,
# please unmark the symbol # of each line below

#library(plm)
#data("Produc")
#dat=Produc
#attach(dat)
#nt <- nrow(dat)
#n <- 48;t=nt/n

#dep<- log(gsp)       # dependent variable             
#ind1<- cbind(log(pcap),log(pc))      # regime-dep:  
#ind2 <- cbind(log(emp),unemp)        # regime-indep:
#d <- log(emp)                        # Threshold variable  

#bootn=c(100,100,100)
#trimn=c(0.01,0.01,0.05)

#qn=400
#conf_lev=0.95
#max_lag <- 1

# ptm(dep,ind1,ind2,d,bootn,trimn,qn,conf_lev,max_lag,t,n)
}