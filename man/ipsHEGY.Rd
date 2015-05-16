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
ipsHEGY(data, itsd, Sel, pmax, CIPS = TRUE)

}

\arguments{
  \item{data}{
Panel data, T by N

}
  \item{itsd}{Options for c(i,t,sd).\cr 
  i=1, intercept;=0 no intercept.\cr
  t=1, trend;=0 no deterministic trend.\cr
	sd=1, season dummy 1:(s-1);=0 no.

}
  \item{Sel}{
Selection of lags, having three options: "signf","bic","aic".

}  \item{pmax}{
Maximum number of lags for searching optimal criteria.

}
  \item{CIPS}{
Logical. If TRUE, using Pesaran(2007) to account for cross-section correlation. The default is TRUE.

}
}
\details{
Mode for selectlags has three options, AIC and BIC use R built-in functions for linear model and their meanings are popular and straightforward. "signf" includes only statistically significant lags, and statistically insignificant lags are dropped from regressors. That is, once you select this option, lags of your model may not be continuous.\cr
The critical values for panel HEGY are standard normal for individual t-ratios, however,you need to perform simulation for the critical values of F joint test, at pdR 1.3. To this end, you are encouraged to work this out for yourself: using arima.sim() to sample seasonal time series with unit root (1-order difference) and obtain their statistics under the null using ipsHEGY(), then it is streightforward to obtain critical values.\cr
Otero et al. (2007) provide critical values for quarterly frequency.
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

\author{Ho Tsungwu <tsungwu@cc.shu.edu.tw>
}


\examples{
data(inf_Q)
dataz<-inf_Q
itsd<-c(1,0,c(1:3))
#Seasonal dummy only takes quarters 1:3, 
#becasue of the presence of common intercept.
Sel<-"bic" # "aic","bic", "signf".
pmax<-12

OUT<-ipsHEGY(dataz,itsd,Sel,pmax,CIPS=FALSE)
OUT$P_HEGY
OUT$U_HEGY

# Simulation of critical values





}