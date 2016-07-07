\name{model}
\alias{model}

\title{Estimate specified panel threshold model

}
\description{This function is the main function estimating threshold regression for function ptm()

}
\usage{
model(r, trim, rep, it, qq1, cf, xt, ct, thresh, tt, qn1, n, qn, cc, yt, ty, k)
}

\arguments{
  \item{r}{vector of threshold estimate(s).

}
  \item{trim}{value of trimmed percentage.

}
  \item{rep}{number bootstrap repetition.

}
  \item{it}{number of regime during computation, used in a for loop.

}
  \item{qq1}{ defined parameter.

}
  \item{cf}{special declaration, e.g. lag().

}

  \item{xt}{regime independent variables.

}
  \item{ct}{trace of regime dependent variables.

}
  \item{thresh}{threshold variable.

}
  \item{tt}{length of time period.

}
  \item{qn1}{as defined by nrow(qq1).

}
  \item{n}{number of cross-section units.
}
  \item{qn}{number of quantiles to examine.
}
  \item{cc}{as defined by 2*log(1-sqrt(conf_lev)).
}
  \item{yt}{vectorized dependent variable.
}
  \item{ty}{trace of yt.
}
  \item{k}{number of regime-independent independent variables.
}
}
\note{Original code offered by Dr. B. E.Hansen (http://www.ssc.wisc.edu/~bhansen/). 

}

\references{Hanse B. E. (1999) Threshold effects in non-dynamic panels: Estimation, testing and inference. Journal of Econometrics,93, 345-368.
}

