\name{contts}
\alias{contts}
\title{
Function for extracting components from a lm object
}
\description{
Extract the standard error and t-stat of the a-th parameter estimate of a lm object
}
\usage{
contts(lm, a)
}
\arguments{
  \item{lm}{
lm object
}
  \item{a}{
The a-th parameter estimate of a linear model regression
}
}

\value{

\item{se.coef}{The standard error of the selected coefficient}
\item{t.stat}{The t-stat of the selected coefficient}
}
\references{
Javier Lopez-de-Lacalle in https://github.com/cran/uroot/blob/master/R/hegy.R 
}
\author{
Ho Tsungwu <tsungwu@cc.shu.edu.tw>, modifed from Javier Lopez-de-Lacalle
}

\examples{
x=rnorm(100)
y=1+0.2*x+rnorm(100)
LMout=lm(y~x)
contts(LMout,1)

#$se.coef
#[1] 0.1081023

#$t.stat
#(Intercept) 
#   10.60401 

}
