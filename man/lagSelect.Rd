\name{lagSelect}
\alias{lagSelect}
\title{Select the optimal number of lags, given criteria}
\description{Determine the optimal number of lags for dynamic regression
}
\usage{
lagSelect(y, maxp, ic)
}

\arguments{
  \item{ y }{A univariate time series data}
  \item{ maxp }{the max number of lags}
  \item{ ic }{Information criteria, either "AIC" or "BIC"}
}
\details{
Information criteria "AIC" and "BIC" use the R built-in functions.
}
\value{
It returns an integer, indicating the optimal lags
}

\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\examples{
#library(pdR)
#data(inf19)
#y<-inf19[,1]
#lagSelect(y,maxp=25,ic="BIC")
}

