\name{SPCShewNormalCenterScale-class}
\Rdversion{1.1}
\docType{class}
\alias{SPCShewNormalCenterScale-class}

\title{Class \code{"SPCShewNormalCenterScale"}}
\description{
Implements a Shewart chart, where the observations are centered and
scaled and past data is resampled parametrically from an estimated
normal distribution.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SPCShewNormalCenterScale", ...)}.
}
\section{Extends}{
Class \code{"\linkS4class{SPCShew}"}, directly.
Class \code{"\linkS4class{SPCchart}"}, by class "SPCShew", distance 2.
}
\section{Methods}{
  \describe{
    \item{getcdfupdates}{\code{signature(chart = "SPCShewNormalCenterScale")}: ... }
    \item{Pofdata}{\code{signature(chart = "SPCShewNormalCenterScale")}: ... }
    \item{resample}{\code{signature(chart = "SPCShewNormalCenterScale")}: ... }
    \item{updates}{\code{signature(chart = "SPCShewNormalCenterScale")}: ... }
    \item{xiofP}{\code{signature(chart = "SPCShewNormalCenterScale")}: ... }
	 }
}
\author{
Axel Gandy
}
\seealso{
\code{\linkS4class{SPCShewNonparCenterScale}},
\code{\linkS4class{ARLShew}},
\code{\linkS4class{hitprobShew}},
\code{\linkS4class{calARLShew}},
\code{\linkS4class{calhitprobShew}}.
}
\examples{
X<-rnorm(100);

##calibrate to ARL 100
chartShew <- new("SPCShewNormalCenterScale",twosided=TRUE)
SPCproperty(data=X,nrep=500,
            property=new("calARLShew",chart=chartShew,target=100),covprob=c(0.7,0.9))

chartShewOneSided <- new("SPCShewNormalCenterScale",twosided=FALSE)
SPCproperty(data=X,nrep=500,
            property=new("calARLShew",chart=chartShewOneSided,target=100),covprob=c(0.7,0.9))


##calibrate to a hitting probability of 0.01 in 100 steps
SPCproperty(data=X,nrep=500,
            property=new("calhitprobShew",chart=chartShew,target=0.01,nsteps=100))
SPCproperty(data=X,nrep=500,
            property=new("calhitprobShew",chart=chartShewOneSided,target=0.01,nsteps=100))

## work out  for ARL for a fixed threshold of 4
SPCproperty(data=X,nrep=500,
            property=new("ARLShew",chart=chartShew,threshold=4))
SPCproperty(data=X,nrep=500,
            property=new("ARLShew",chart=chartShewOneSided,threshold=4))


SPCproperty(data=X,nrep=500,
            property=new("hitprobShew",chart=chartShew,nsteps=100,threshold=4))

SPCproperty(data=X,nrep=500,
            property=new("hitprobShew",chart=chartShewOneSided,nsteps=100,threshold=4))

}
\keyword{classes}
