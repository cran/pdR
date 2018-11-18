\name{iClick.plm1way}
\alias{iClick.plm1way}

\title{
iClick GUI for one-way panel data analysis, based on package plm.
}

\encoding{latin1}

\description{
This function generates analysis of panel data by iClick.plm. 
}


\usage{
iClick.plm1way(dep,indep,Formula,data,bootrep=99,ENDOG,IV,inst.method)
}

\arguments{
  \item{dep}{Column number of dependent variable; e.g., dep=data[,2]. Default is NULL

}

  \item{indep}{Column number of Independent variables; e.g., indep=data[,c(3,5,8)]. Default is NULL

}

  \item{Formula}{Equation input by explicit formula; e.g., y=x1+x2+x3.Default is NULL

}

  \item{data}{ A panel data class declared by plm.

}

  \item{bootrep}{ Replication number of bootstrapping for fixed effect, the default number is 99 to avoid unnecessary computation.

}

  \item{ENDOG}{For 2SLS, declare endogeneous variables here; otherwise, keep it as default by NULL.

}

  \item{IV}{For 2SLS, declare IV variables here; otherwise, keep it as default by NULL.

}

  \item{inst.method}{For 2SLS, select estimation method.

}


}

\value{
GUI output button.
}

\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}


\seealso{
Package plm.
}
\examples{
#library("pdR")
##data("bank_income")
#data1.plm=plm.data(bank_income,index="ID");
#head(data1.plm,2)
#iClick.plm1way(dep=4,indep=c(5,7,8),data=data1.plm)


#data("productivity")
#data2.plm=plm.data(productivity,index="state")
#head(data2.plm,2)
#formula2="log(y_gsp)~log(x1_hwy)+log(x2_water)"
#iClick.plm2way(Formula=formula2,data=data2.plm)


#data("crime")
#data3.plm= plm.data(crime, index = c("county"));head(data3.plm,2)
#formula3="log(crmrte)~log(prbarr)+log(polpc)+log(prbconv)+log(prbpris)
#+log(avgsen)+log(density)+log(wcon)+log(wtuc)+log(wtrd)+log(wfir)
#+log(wser) +log(wmfg)+log(wfed)+log(wsta)+log(wloc)+log(pctymle)
#+log(pctmin)+smsa+region"
#endo=c("log(prbarr)","log(polpc)")
#iv=c("log(taxpc)","log(mix)")
#iClick.plm1way(Formula=formula3,data=data3.plm,ENDOG=endo,IV=iv)

}

