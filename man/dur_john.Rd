\name{dur_john}
\alias{dur_john}
\docType{data}
\title{
The cross-country growth data in Durlauf and Johnson(1995)
}

\description{
The Durlauf-Johnson data manupipulated by Hansen(2000),excluding missing variables and oil states 
}
\usage{data(dur_john)
}


\format{
  A data frame with 19 countries
  \describe{
    \item{\code{gdpGrowth}}{Economic growth measured by GDP of 1960 and 1985}
    \item{\code{logGDP60}}{log Per capita GDP in 1960}
    \item{\code{Inv_GDP}}{Average ratio of investment (including Government Investment) to GDP      from 1960 to 1985}
    \item{\code{popGrowth}}{Average growth rate of working-age population 1960 to 1985}
    \item{\code{School}}{Average fraction of working-age population enrolled in secondary school
       from 1960 to 1985}
    \item{\code{GDP60}}{Per capita GDP in 1960}
    \item{\code{Literacy}}{fraction of the population over 15 years old that is able to read 
      and write in 1960}
  }
}

\details{
Steven N. Durlauf and Paul A. Johnson, "Multiple Regimes and Cross-Country
Growth Behavior," Journal of Applied Econometrics, Vol. 10, No. 4, 1995,
365-384. 
}

\examples{
data(dur_john)
head(dur_john)
}
