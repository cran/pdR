\name{model}
\alias{model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate specified panel threshold model
%%  ~~function to do ... ~~
}
\description{This function is the main function estimating threshold regression for function ptm()
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
model(r, trim, rep, it, qq1, cf, xt, ct, thresh, tt, qn1, n, qn, cc, yt, ty, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{vector of threshold estimate(s)
%%     ~~Describe \code{r} here ~~
}
  \item{trim}{value of trimmed percentage
%%     ~~Describe \code{trim} here~~
}
  \item{rep}{number bootstrap repetition
%%     ~~Describe \code{rep} here ~~
}
  \item{it}{number of regime during computation, used in a for loop
%%     ~~Describe \code{it} here~~
}
  \item{qq1}{ defined parameter
%%     ~~Describe \code{qq1} here~~
}
  \item{cf}{special declaration, e.g. lag()
%%     ~~Describe \code{cf} here~~
}

  \item{xt}{regime independent variables
%%     ~~Describe \code{xt} here~~
}
  \item{ct}{trace of regime dependent variables
%%     ~~Describe \code{ct} here~~
}
  \item{thresh}{threshold variable
%%     ~~Describe \code{threshold} here~~
}
  \item{tt}{length of time period
%%     ~~Describe \code{tt} here~~
}
  \item{qn1}{as defined by nrow(qq1)
%%     ~~Describe \code{qn1} here~~
}
  \item{n}{number of cross-section units
%%     ~~Describe \code{n} here~~
}
  \item{qn}{number of quantiles to examine
%%     ~~Describe \code{qn} here~~
}
  \item{cc}{as defined by 2*log(1-sqrt(conf_lev))
%%     ~~Describe \code{cc} here~~
}
  \item{yt}{vectorized dependent variable
%%     ~~Describe \code{yt} here~~
}
  \item{ty}{trace of yt
%%     ~~Describe \code{ty} here~~
}
  \item{k}{number of regime-independent independent variables
%%     ~~Describe \code{k} here~~
}
}
\note{Original code offered by Dr. B. E.Hansen (http://www.ssc.wisc.edu/~bhansen/). 
%%  ~~further notes~~
}

\references{Hanse B. E. (1999) Threshold effects in non-dynamic panels: Estimation, testing and inference. Journal of Econometrics,93, 345-368.
}

