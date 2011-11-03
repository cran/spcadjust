\name{calARLCUSUM-class}
\Rdversion{1.1}
\docType{class}
\alias{calARLCUSUM-class}

\title{Class \code{"calARLCUSUM"}}
\description{
Computes for a CUSUM chart the threshold needed to give a desired ARL.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("calARLCUSUM", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{gridpoints}:}{Object of class \code{"numeric"}  }
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
    \item{SPCq}{\code{signature(property = "calARLCUSUM")}: ... }
	 }
}
\author{
Axel Gandy
}


\seealso{
 \code{\linkS4class{ARLCUSUM}} gives the inverse.
}
\examples{
X <-  rnorm(1000)
chart <- new("SPCCUSUMNormal",Delta=1)
SPCproperty(data=X,nrep=10,
            property=new("calARLCUSUM",chart=chart,target=100))
#increase the number of repetitions nrep for real applications.
}
\keyword{classes}
