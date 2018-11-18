\name{selPabic}
\alias{selPabic}
\title{
Selection of lags.
}
\description{
 Lagged coefficient estimates are kept if they meet the inequality condition of AIC or BIC.
}
\usage{
selPabic(lmdet, type, Pmax = NULL)
}

\arguments{
  \item{lmdet}{
Object of lm()
}
  \item{type}{Take the value of "aic" or "bic".
}
  \item{Pmax}{
The maximum number of lag orders.
}
}
\details{
This is an internal function used for HEGY.test().  Beginning with pamx, the lag order will be drop  if its inclusion worsens the minimum condition. Hence, they may not be a regular sequence.  For example, for pmax=10, the selected lags may look like (1,4,5,8,9), rather than 1,2,3,...10.
}
\value{
This function returns the lag orders.
}
\references{
Javier Lopez-de-Lacalle in https://github.com/cran/uroot/blob/master/R/hegy.R 
}
\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, modifed from Javier Lopez-de-Lacalle
}

\examples{
data(inf_Q)
y=inf_Q[,1]
hegy.out<-HEGY.test(wts=y, itsd=c(1,0,c(1:3)),regvar=0, selectlags=list(mode="aic", Pmax=12))
hegy.out$lagsorder
hegy.out$lagcoefs
}