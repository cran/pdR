\name{lookupCVtable}
\alias{lookupCVtable}
\title{
Function for looking up tabulated critical values and associated p-values of HEGY test.
}
\description{
Function for looking up tabulated critical values and associated p-values,Hylleberg et. al (1990, Table 1a and Table 1b).
}
\usage{
lookupCVtable(code)
}
\arguments{
  \item{code}{
Type of HEGY model, this will be automatically identified.
}
}
\value{

\item{table }{Table for critical value and p-value.}

}

\references{
Hylleberg, S., Engle, R.F., Granger, C.W.J., and Yoo, B.S.(1990) Seasonal integration and cointegration. Journal of Econometrics,44, 215-238.\cr
Javier Lopez-de-Lacalle in https://github.com/cran/uroot/blob/master/R/hegy.R
}
\author{
Ho Tsungwu <tsungwu@cc.shu.edu.tw>, modifed from Javier Lopez-de-Lacalle
}
