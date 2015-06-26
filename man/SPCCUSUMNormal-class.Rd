\name{SPCCUSUMNormal-class}
\Rdversion{1.1}
\docType{class}
\alias{SPCCUSUMNormal-class}

\title{Class \code{"SPCCUSUMNormal"}}
\description{
CUSUM chart with centering and scaling assuming normally distributed observations.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SPCCUSUMNormal", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{Delta}:}{The amount subtracted at each step of the CUSUM
      chart}
  }
}
\section{Extends}{
Class \code{"\linkS4class{SPCCUSUM}"}, directly.
Class \code{"\linkS4class{SPCchart}"}, by class "SPCCUSUM", distance 2.
}
\section{Methods}{
  \describe{
    \item{getcdfupdates}{\code{signature(chart = "SPCCUSUMNormal")}: ... }
    \item{Pofdata}{\code{signature(chart = "SPCCUSUMNormal")}: ... }
    \item{resample}{\code{signature(chart = "SPCCUSUMNormal")}: ... }
    \item{updates}{\code{signature(chart = "SPCCUSUMNormal")}: ... }
    \item{xiofP}{\code{signature(chart = "SPCCUSUMNormal")}: ... }
	 }
}
\author{
Axel Gandy
}


\seealso{
  \code{\linkS4class{SPCCUSUMNonparCenterScale}},
  \code{\linkS4class{calARLCUSUM}},
  \code{\linkS4class{calhitprobCUSUM}},
  \code{\linkS4class{hitprobCUSUM}},
  \code{\linkS4class{ARLCUSUM}}.
}
\examples{
##increase nrep for real applications
X <-  rnorm(1000)
chart <- new("SPCCUSUMNormal",Delta=1)
SPCproperty(data=X,nrep=10,
            property=new("calARLCUSUM",chart=chart,target=100))

SPCproperty(data=X,nrep=10,
            property=new("calhitprobCUSUM",chart=chart,target=0.05,nsteps=1000))

SPCproperty(dat=X,nrep=10,
            property=new("hitprobCUSUM",chart=chart,threshold=5,nsteps=100),
            covprob=c(0.7,0.8,0.9))
SPCproperty(dat=X,nrep=10,
            property=new("ARLCUSUM",chart=chart,threshold=2),
            covprob=c(0.7,0.8,0.9))
}
\keyword{classes}
