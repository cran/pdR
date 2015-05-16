\name{pdR-package}
\alias{pdR-package}
\alias{pdR}
\docType{package}
\title{
Panel Data Regression: Threshold Model and Unit Root Tests

}
\encoding{latin1}

\description{Functions for analysis of panel data, including the panel threshold model of Hansen (1999,JE), panel unit root test of Chang(2002,JE) based upon instuments generating functions (IGF), and panel seasonal unit root test based upon Hylleberg et al.(1990,JE) .

}

\details{ This package designs a specification function ptm() to estimate the panel threshold model of Hansen(1999). The key feature of ptm() is to generalize Hansen's original code to allow multiple (more-than-one) regime-dependent right-hand-side independent variables; Dr. Hansen's original code admits only 1 regime-dependent right-hand-side independent variable. This version also includes panel unit root tests based on the instrument generating functions(IGF), proposed by Chang (2002, J. of Econometrics), and the panel version of Hylleberg et al.(1990) seasonal unit root test, proposed by Otero, et al. (2005, 2007).
\tabular{ll}{
Package: \tab pdR\cr
Type: \tab Package\cr
Version: \tab 1.3\cr
Date: \tab 2015-5-16\cr
License: \tab GPL-2\cr
}

}
\references{ 
Chang, Y. (2002) Nonlinear IV Unit Root Tests in Panels with Cross-Sectional Dependency. Journal of Econometrics, 110, 261-292.\cr
Hansen B. E. (1999) Threshold effects in non-dynamic panels: Estimation, testing and inference. Journal of Econometrics,93, 345-368.\cr
Hylleberg, S., Engle, R.F., Granger, C.W.J., and Yoo, B.S.(1990) Seasonal integration and cointegration. Journal of Econometrics,44, 215-238.\cr
Otero, J., Smith, J., and Giulietti, M. (2005) Testing for seasonal unit roots in heterogeneous panels. Economics Letters,  86, 229-235.\cr
Otero, J., Smith, J., and Giulietti, M. (2007) Testing for seasonal unit roots in heterogeneous panels in the presence of cross section dependence. Economics Letters, 86, 179-184.\cr
Pesaran M. Hashem (2007) A simple panel unit root test in the presence of cross-section dependence.  Journal of Applied Econometrics, 22, 265-312.
}

\author{
Ho Tsung-wu <tsungwu@cc.shu.edu.tw>\cr
Maintainer: Ho Tsung-wu <tsungwu@cc.shu.edu.tw>
}