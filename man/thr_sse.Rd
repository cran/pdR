\name{thr_sse}
\alias{thr_sse}

\title{a subroutine calculating SSE

}
\description{This function is a sub-routine for model(), calculating SSE of each threshold regression.

}
\usage{
thr_sse(y, q, r, cf, xt, ct, thresh, tt, n)
}

\arguments{
  \item{y}{parameter.
}
  \item{q}{qq1 in model().
}
  \item{r}{parameter.
}
  \item{cf}{as defined in model().
}
  \item{xt}{as defined in model().
}
  \item{ct}{as defined in model().
}
  \item{thresh}{as defined in model().
}
  \item{tt}{as defined in model().
}
  \item{n}{as defined in model().
}

}


\references{Hanse B. E. (1999) Threshold effects in non-dynamic panels: Estimation, testing and inference. Journal of Econometrics,93, 345-368.

Original code from Dr. Hansen (http://www.ssc.wisc.edu/~bhansen/).
}

