\name{r_est}
\alias{r_est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A subroutine for model()
%%  ~~function to do ... ~~
}
\description{This function is a subroutine for model(), estimation procedure
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
r_est(y, r, trim, tt, qq1, qn1, qn, n, cf, xt, ct, thresh)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{vector of dependent variable
%%     ~~Describe \code{y} here ~~
}

  \item{r}{numer of regime 
%%     ~~Describe \code{r} here ~~
}

  \item{trim}{value of trimmed percentage
%%     ~~Describe \code{trim} here ~~
}

  \item{tt}{length of time period
%%     ~~Describe \code{tt} here ~~
}

  \item{qq1}{parameter defined by as.matrix(unique(thresh)[floor(sq*nrow(as.matrix(sort(unique(thresh)))))])
%%     ~~Describe \code{qq1} here ~~
}
  \item{qn1}{as defined by nrow(qq1)
%%     ~~Describe \code{qn1} here ~~
}
  \item{qn}{number of quantiles to examine
%%     ~~Describe \code{qn} here ~~
}
  \item{n}{parameter of cross-section units
%%     ~~Describe \code{n} here ~~
}
  \item{cf}{special declaration, e.g. lag()
%%     ~~Describe \code{cf} here ~~
}
  \item{xt}{regime independent variables
%%     ~~Describe \code{xt} here ~~
}
  \item{ct}{trace of regime dependent variables
%%     ~~Describe \code{ct} here ~~
}
  \item{thresh}{threshold variable
%%     ~~Describe \code{thresh} here ~~
}


}

\references{Hanse B. E. (1999) Threshold effects in non-dynamic panels: Estimation, testing and inference. Journal of Econometrics,93, 345-368.

Original code from Dr. Hansen (http://www.ssc.wisc.edu/~bhansen/).
%% ~put references to the literature/web site here ~
}
