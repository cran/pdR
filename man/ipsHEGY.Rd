\name{ipsHEGY}
\alias{ipsHEGY}

\title{
IPS-HEGY seasonal unit root test in panel data, Otero et al.(2007).

}
\encoding{latin1}

\description{
This function performs panel data-based HEGY seasonal unit root test, the asymptotics is based upon Otero et al.(2007).

}

\usage{
ipsHEGY(data,deterministic =c(1,0,0),lag.method="AIC",maxlag=12,CIPS = TRUE)

}

\arguments{
  \item{data}{
Panel data, T by N

}
  \item{deterministic}{Options for deterministic components.\cr 
  c(1,0,0), intercept only.\cr
  c(0,1,0), trend only.\cr
	c(0,0,1),, season dummy only.\cr
  

}
  \item{lag.method}{
Selection of optimal lags, having for options: "fixed","AIC", "BIC", "AICc".

}  \item{maxlag}{
Maximum number of lags for searching optimal criteria.

}
  \item{CIPS}{
Logical. If TRUE, using Pesaran(2007) to account for cross-section correlation. The default is TRUE.

}
}
\details{
Mode for selectlags has four options, AIC,AICc and BIC use R built-in functions for linear model and their meanings are popular and straightforward. "fixed" fixes maxlag as the number of lags. \cr
The critical values for panel HEGY are standard normal for individual t-ratios, however,you need to perform simulation for the critical values of F joint test, at pdR 1.3. To this end, you are encouraged to work this out for yourself: using arima.sim() to sample seasonal time series with unit root (1-order difference) and obtain their statistics under the null using ipsHEGY(), then it is streightforward to obtain critical values.\cr
Otero et al. (2007) provide critical values for quarterly frequency. \cr
The univariate HEGY is based on R package \code{uroot}, simply modified to include cross-secionaal average.
}
\value{

\item{P_HEGY}{Panel HEGY statistics.}
\item{U_HEGY}{Idividual HEGY statistics of N units.}
}

\references{
Otero, J., Smith, J., and Giulietti, M. (2005) Testing for seasonal unit roots in heterogeneous panels. Economics Letters,  86, 229-235.\cr
Otero, J., Smith, J., and Giulietti, M. (2007) Testing for seasonal unit roots in heterogeneous panels in the presence of cross section dependence. Economics Letters, 86, 179-184.\cr
Pesaran M. Hashem (2007) A simple panel unit root test in the presence of cross-section dependence.  Journal of Applied Econometrics, 22, 265-312.
}

\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}


\examples{
data(inf_M)

#Seasonal dummy only takes quarters 1:3, 
#becasue of the presence of common intercept.

OUT<-ipsHEGY(inf_M, deterministic = c(1,1,1), lag.method = "AIC", maxlag = 12,CIPS=TRUE)
OUT$P_HEGY
OUT$U_HEGY


# Simulation of critical values



}