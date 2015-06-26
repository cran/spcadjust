\name{calhitprobShew-class}
\Rdversion{1.1}
\docType{class}
\alias{calhitprobShew-class}

\title{Class \code{"calhitprobShew"}}
\description{
Derives the threshold needed to give a certain false alarm property for
Shewhart charts.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("calhitprobShew", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{target}:}{Object of class \code{"numeric"}  }
    \item{\code{nsteps}:}{Object of class \code{"numeric"}  }
    \item{\code{chart}:}{Object of class \code{"SPCchart"}  }
    \item{\code{lowerconf}:}{Object of class \code{"logical"}  }
  }
}
\section{Extends}{
Class \code{"\linkS4class{SPCcalhitprob}"}, directly.
Class \code{"\linkS4class{SPCproperty}"}, by class "SPCcalhitprob", distance 2.
}
\section{Methods}{
  \describe{
    \item{SPCq}{\code{signature(property = "calhitprobShew")}: ... }
	 }
}
\author{
Axel Gandy
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\linkS4class{hitprobShew}},
 \code{\linkS4class{calhitprobCUSUM}}.
}
\examples{
chartShew <- new("SPCShewNormalCenterScale")
X <- rnorm(100)
SPCproperty(data=X,nrep=500,
            property=new("calhitprobShew",chart=chartShew,target=0.01,nsteps=100))
}
\keyword{classes}
