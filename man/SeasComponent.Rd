\name{SeasComponent}
\alias{SeasComponent}

\title{
Generate a data matrix of seasonal components
}
\description{
Generate a data matrix of seasonal components, having two pattern cycles.
}
\usage{
SeasComponent(wts, type)
}

\arguments{
  \item{wts}{
A univeriate time series with monthly or quarterly frequency.
}
  \item{type}{Types of patterns of seasonal cycle.\cr
="dummyCycle", generating dummy variables for the pattern of seasonal cycle, Barsky & Miron (1989)\cr
="trgCycle", generating trigonometric variables for the pattern of seasonal cycle, Granger & Newbold (1986).
}
}
\details{
This function generates data matrix for controlling the pattern of seasonal cycles. type="dummyCycle" generates DUMMY variables with season frequency. type="trgCycle" generates trigonometric pattern.  
}
\value{
A dataframe returns. Number of columns is determined by frequency.
}
\references{
Barsky, Robert B. and Jeffrey A. Miron (1989) The Seasonal Cycle and the Business Cycle. Journal of Political Economy, 97 (3): 503-32.\cr
Granger, Clive William John and Newbold, Paul (1986) Forecasting Economic Time Series. 2nd edition. Published by New Milford, Connecticut, U.S.A.: Emerald Group Pub Ltd.\cr
Javier Lopez-de-Lacalle in https://github.com/cran/uroot/blob/master/R/hegy.R 
}
\author{
Ho Tsungwu <tsungwu@cc.shu.edu.tw>, modifed from Javier Lopez-de-Lacalle
}
\examples{
data(inf_Q)
y=inf_Q[,2]
SeasComponent(y,type="dummyCycle")
SeasComponent(y,type="trgCycle")
}

