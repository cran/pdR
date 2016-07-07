\name{cigaretts}
\alias{cigaretts}

\title{
Cigaretts consumption of US states 
}
\encoding{latin1}

\description{Cigaretts consumption of US states }
\usage{
data(cigaretts)
}

\format{
  A data frame of 48 US states' cigaretts consumption

  \describe{

    \item{\code{State}}{State abbreviation, N}
    \item{\code{Year}}{Year, t}
    \item{\code{Y_SALES}}{Cigarette sales in packs per capita, deflated by population}
    \item{\code{X1_PRICE}}{P=Real price per pack of cigarettes, deflated by 1983 CPI.
}
    \item{\code{X2_PIMIN}}{Real minimum price in adjoining states per pack of cigarettes, deflated by CPI}
    \item{\code{X3_NDI}}{Per capita disposable income}



  }
}


\references{
Baltagi Badi H. (2005) Econometric Analysis of Panel Data.  John Wiley.
}

\examples{
data(cigaretts)
head(cigaretts)
}