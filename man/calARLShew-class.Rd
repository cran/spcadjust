\name{calARLShew-class}
\Rdversion{1.1}
\docType{class}
\alias{calARLShew-class}

\title{Class \code{"calARLShew"}}
\description{
Calibrates a Shewhart chart to a given in-control ARL.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("calARLShew", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{target}:}{Object of class \code{"numeric"}  }
    \item{\code{chart}:}{Object of class \code{"SPCchart"}  }
    \item{\code{lowerconf}:}{Object of class \code{"logical"}  }
  }
}
\section{Extends}{
Class \code{"\linkS4class{SPCcalARL}"}, directly.
Class \code{"\linkS4class{SPCproperty}"}, by class "SPCcalARL", distance 2.
}
\section{Methods}{
  \describe{
    \item{SPCq}{\code{signature(property = "calARLShew")}: ... }
	 }
}
\author{
Axel Gandy
}

\seealso{
 \code{\linkS4class{ARLShew}}.
}
\examples{
X<-rnorm(100)
chartShew <- new("SPCShewNormalCenterScale")
SPCproperty(data=X,nrep=500,
            property=new("calARLShew",chart=chartShew,target=1000))
}
\keyword{classes}
