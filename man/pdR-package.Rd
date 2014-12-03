\name{pdR-package}
\alias{pdR-package}
\alias{pdR}
\docType{package}

\title{
Panel Data Regression for threshold model and unit root tests
}
\encoding{latin1}

\description{
Panel threshold regression (Hansen,1999) and IGF-based panel unit root test (Chang,2002)

}

\details{ 
This package designs a specification function ptm() to estimate the panel threshold model of Hansen(1999). The key feature of ptm() is to generalize Hansen's original code to allow multiple (more-than-one) regime-dependent independent variables; Dr. Hansen's original code admits only 1 regime-dependent independent variable. This version also includes panel unit root test based on the instrument generating functions(IGF), proposed by Chang (2002, J. of Econometrics)

\tabular{ll}{
Package: \tab pdR \cr
Type: \tab Package \cr
Version: \tab 1.2 \cr
Date: \tab 2014-12-03 \cr
License: \tab GPL-2 \cr
}

}

\references{ 
Chang, Y. (2002) Nonlinear IV Unit Root Tests in Panels with Cross-Sectional Dependency.  \emph{Journal of Econometrics}, 110, 261¡V292.\cr
Hanse B. E. (1999) Threshold effects in non-dynamic panels: Estimation, testing and inference. \emph{Journal of Econometrics},93, 345-368.
}

\author{
Ho Tsung-wu <tsungwu@cc.shu.edu.tw>\cr
Maintainer: Ho Tsung-wu <tsungwu@cc.shu.edu.tw>
}