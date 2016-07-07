\name{r_est}
\alias{r_est}

\title{A subroutine for model()

}
\description{This function is a subroutine for model(), estimation procedure.

}
\usage{
r_est(y, r, trim, tt, qq1, qn1, qn, n, cf, xt, ct, thresh)
}

\arguments{
  \item{y}{vector of dependent variable.

}

  \item{r}{numer of regime. 

}

  \item{trim}{value of trimmed percentage.

}

  \item{tt}{length of time period.

}

  \item{qq1}{parameter defined by as.matrix(unique(thresh)[floor(sq*nrow(as.matrix(sort(unique(thresh)))))]).

}
  \item{qn1}{as defined by nrow(qq1).
}
  \item{qn}{number of quantiles to examine.
}
  \item{n}{parameter of cross-section units.
}
  \item{cf}{special declaration, e.g. lag().
}
  \item{xt}{regime independent variables.
}
  \item{ct}{trace of regime dependent variables.
}
  \item{thresh}{threshold variable.
}


}

\references{Hanse B. E. (1999) Threshold effects in non-dynamic panels: Estimation, testing and inference. Journal of Econometrics,93, 345-368.

Original code from Dr. Hansen (http://www.ssc.wisc.edu/~bhansen/).
}
