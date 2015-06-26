\name{SPCproperty}
\alias{SPCproperty}
\title{
Computes bootstrap adjusted properties for control charts.
}
\description{
Computes bootstrap adjusted properties for control charts.
}
\usage{
SPCproperty(data, nrep = 500, property, covprob = 0.9)
}
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
\value{An object of type SPCpropertyres.
}
\author{
Axel Gandy
}


\seealso{
 \code{\link{SPC2sidedconfint}}
}
\examples{
 # calibrate CUSUM chart to an in-control ARL of 100
 # run with a larger number of replications in real examples!

 X <-  rnorm(100) #observed data
 chart <- new("SPCCUSUMNormal",Delta=1) # CUSUM chart with normal observations
 SPCproperty(data=X,nrep=15,
            property=new("calARLCUSUM",chart=chart,target=100))
}
