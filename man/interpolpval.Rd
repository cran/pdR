\name{interpolpval}
\alias{interpolpval}

\title{
Extracting critical value and p-value from Table 1 of Hylleberg et. al (1990)
}
\description{
Hylleberg et. al (1990,pp.226-227) offer simulated critical values for seasonal unitr  to test. interpolpval() is an internal call and should not be used independently.
}
\usage{
interpolpval(code, stat, N, swarn = TRUE)
}

\arguments{
  \item{code}{
Type of HEGY model, this will be automatically identified.
}
  \item{stat}{
Empirical test statistics.
}
  \item{N}{
Sample size calculating stat above.
}
  \item{swarn}{
Logical. Whether the warning message for negative p-value will be returned? The default is TRUE.
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
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, modifed from Javier Lopez-de-Lacalle
}
