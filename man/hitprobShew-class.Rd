\name{hitprobShew-class}
\Rdversion{1.1}
\docType{class}
\alias{hitprobShew-class}

\title{Class \code{"hitprobShew"}}
\description{
Computes the hitting probability for Shewhart charts.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("hitprobShew", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{threshold}:}{Object of class \code{"numeric"}  }
    \item{\code{nsteps}:}{Object of class \code{"numeric"}  }
    \item{\code{chart}:}{Object of class \code{"SPCchart"}  }
    \item{\code{lowerconf}:}{Object of class \code{"logical"}  }
  }
}
\section{Extends}{
Class \code{"\linkS4class{SPChitprob}"}, directly.
Class \code{"\linkS4class{SPCproperty}"}, by class "SPChitprob", distance 2.
}
\section{Methods}{
  \describe{
    \item{SPCq}{\code{signature(property = "hitprobShew")}: ... }
	 }
}
\author{
  Axel Gandy
}
\seealso{
 \code{\linkS4class{calhitprobShew}}.
}
\examples{
X <- rnorm(100)
chartShew <- new("SPCShewNormalCenterScale")
SPCproperty(data=X,nrep=500,
            property=new("hitprobShew",chart=chartShew,nsteps=100,threshold=4))
}
\keyword{classes}
