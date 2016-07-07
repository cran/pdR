\name{productivity}
\alias{productivity}

\title{
Productivity data of 48 US state,1970-1986 
}
\description{
Gross state production data
}
\usage{
data(productivity)
}
							
\format{
  A data frame with US production
\describe{
  \item{state}{US state index, 1-48
}
  \item{year}{Year index
}

  \item{y_gsp}{ Gross state product
}
  \item{x1_hwy}{Expediture of public utility- highway construction
}
  \item{x2_water}{Expediture of public utility- water
}
  \item{x3_other}{Expediture of others}

  \item{x4_private}{Private consumption of each state}

  \item{x5_emp}{Employment rate of each state}

  \item{x6_unemp}{Unemployment rate of each state}
}
}

\examples{
data(productivity)
head(productivity)
}
