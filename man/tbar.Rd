\name{tbar}
\alias{tbar}
\title{
Compute the resursive mean
}
\description{Compute the resursive mean of each series
}

\usage{ tbar(x)
}

\arguments{
  \item{ x }{A univariate time series data}

}

\details{
This function computes the resursive mean
}

\author{Ho Tsung-wu <tsungwu@ntnu.edu.tw>
}


\examples{
data(inf19)
y <- inf19[,1]
tbar(y)

}
