\name{crime}
\alias{crime}

\title{
Annual crime dataset of US counties 
}

\usage{
data(crime)
}
\encoding{latin1}
\description{Annual crime dataset of US counties }

\format{
  A data frame of US counties
  \describe{
    \item{\code{county}}{counties index, N}
    \item{\code{year}}{Year, t}
    \item{\code{crmrte}}{crime rate(crime/population)}
    \item{\code{prbarr}}{probability of arrest (arrests/offenses)}
    \item{\code{prbconv}}{probability of conviction, given arrest}
    \item{\code{prbpris}}{probability of a prison, given conviction }
    \item{\code{avgsen}}{sanction severity(average prison sentence in days )}
    \item{\code{polpc}}{ability of police force to detect crime(# of police per capita)}
    \item{\code{density}}{population density(POP/area)
    }
    \item{\code{taxpc}}{Taxpayment per capita}
    \item{\code{region}}{region index of county}
    \item{\code{smsa}}{=1 if SAMA, POP>50000; =0 else }
    \item{\code{pctmin}}{See Baltagi(2006) for details}
    \item{\code{wcon}}{See Baltagi(2006) for details }
    \item{\code{wtuc}}{See Baltagi(2006) for details}
    \item{\code{wtrd}}{See Baltagi(2006) for details}
    \item{\code{wfir}}{See Baltagi(2006) for details}
    \item{\code{wser}}{See Baltagi(2006) for details}
    \item{\code{wmfg}}{See Baltagi(2006) for details}
    \item{\code{wfed}}{See Baltagi(2006) for details}
    \item{\code{wsta}}{See Baltagi(2006) for details}
    \item{\code{wloc}}{See Baltagi(2006) for details}
    \item{\code{mix}}{See Baltagi(2006) for details}
    \item{\code{pctymle}}{See Baltagi(2006) for details}
   
  }
}


\references{
Baltagi Badi H. (2005) Econometric Analysis of Panel Data.  John Wiley.
Baltagi Badi H. (2006) Estimating an Economic Model of Crime Using Panel Data from North Carolina. J.of Applied Econometrics 21: 543¡V547.
}

