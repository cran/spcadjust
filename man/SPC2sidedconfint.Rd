\name{SPC2sidedconfint}
\alias{SPC2sidedconfint}
\title{
Computes a two-sided confidence interval for properties of a CUSUM chart.
}
\description{
Computes a two-sided confidence interval for properties of a CUSUM chart.
}
\usage{
SPC2sidedconfint(data, nrep = 500, property, covprob = 0.9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{data}{ The observed data.
}
  \item{nrep}{The number of bootstrap repetitions.
}
  \item{property}{The property to be computed. Must be a descendant of SPCproperty.
}
\item{covprob}{The coverage probability of the adjustment.
}
}
\value{The desired confidence interval, a vector of length 2.
}
\author{
Axel Gandy
}

\seealso{\code{\link{SPCproperty}}
}
\examples{
# compute 2-sided CI for the ARL of a CUSUM control chart assuming normality
 X <-  rnorm(100) #observed data
 chart <- new("SPCCUSUMNormal",Delta=1) # CUSUM chart with normal observations
 SPC2sidedconfint(data=X,nrep=100,covprob=0.95,
            property=new("ARLCUSUM",chart=chart,threshold=4))
}

