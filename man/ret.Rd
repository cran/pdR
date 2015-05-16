\name{ret}
\alias{ret}
\title{
Returns a data.frame of sequential lag matrix.
}
\description{
ret() is similar to embed(), but returns a data.frame specified with colnames, not matrix.
}
\usage{
ret(wts, k)
}

\arguments{
  \item{wts}{
Univariate time series.
}
  \item{k}{
k-1 lagged terms. 
}
}
\details{
ret() is similar to embed(), but returns a data.frame with colnames, not matrix.  Moreover, unlike embed(),ret() fills lagged cells with NA, instead of trimming them.
}
\value{
A T by k dataframe returns. If you need 2 lags, you have to specify k=3, then it returns a dataframe with T by 3 dataframe, the first column is lag0.

}
\references{
Javier Lopez-de-Lacalle in https://github.com/cran/uroot/blob/master/R/hegy.R 
}
\author{
Ho Tsungwu <tsungwu@cc.shu.edu.tw>, modifed from Javier Lopez-de-Lacalle
}
\examples{
data(inf_Q)
y=inf_Q[,2]
ret(y,3)
}
