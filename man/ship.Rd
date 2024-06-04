\name{ship}
\alias{ship}
\docType{data}
\title{
Panel data on the number of ship accidents
}
\encoding{latin1}
\description{
 Panel data on the number of ship accidents, McCullagh and
Nelder(1989)
}
\usage{data("ship")}

\format{


  \describe{
    \item{\code{accident}}{the number of ship accidents}
    \item{\code{ship}}{Ship iD}
    \item{\code{service}}{the number of months in service}
    \item{\code{op_75_79}}{the operating period between 1975 and 1979}
    \item{\code{co_65_69}}{consecutive construction periods of 5 years}
    \item{\code{co_70_74}}{consecutive construction periods of 5 years}
    \item{\code{co_75_79}}{consecutive construction periods of 5 years}
    \item{\code{yr_con}}{years of construction}
    \item{\code{yr_op}}{years of operation}
  }
}

\references{McCullagh, P., and J. A. Nelder (1989) Generalized Linear Models. 2nd ed. London: Chapman and Hall/CRC.}


\examples{
data(ship)
}

