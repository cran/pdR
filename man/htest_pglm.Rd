\name{htest_pglm}
\alias{htest_pglm}

\title{Specification test for panel glm models

}

\encoding{latin1}


\description{This function performs Hausman specification test for panel glm. 

}
\usage{
htest_pglm(RE, FE, re.method)
}

\arguments{
  \item{RE}{Random effect objects. Support pglm, glmer, glmmTMB}
  \item{FE}{Fixed effect objects.}
  \item{re.method}{Method that used to estimate the random effect estimation, in addition to "pglm", it also supports "glmmTMB" of package glmmTMB, and "glmer" of package lme4.
}

}


\references{
Hausman J.A. (1978). Specification Tests in Econometrics. Econometrica, 46, 1251-1271.
}
\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, modifed from phtest() of plm
}

\examples{
data(ship)
library(pglm)
Eq1="accident  ~  op_75_79+co_65_69+co_70_74+co_75_79"
FE.pois <- pglm(Eq1,data=ship,family = "poisson",model = "within",index = 'ship',R=10)

RE.pois <- pglm::pglm(Eq1,data=ship,family = "poisson", model = "random", index = 'ship')

## Hausman test
htest_pglm(RE=RE.pois,  FE=FE.pois, re.method="pglm")

Eq2=accident  ~  op_75_79+co_65_69+co_70_74+co_75_79 + (1 | ship)
re.glmmTMB=glmmTMB::glmmTMB(Eq2,data=ship, family="poisson")

## Hausman test
htest_pglm(RE=re.glmmTMB,  FE=FE.pois, re.method="glmmTMB")


}