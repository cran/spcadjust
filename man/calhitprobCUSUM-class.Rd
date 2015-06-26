\name{calhitprobCUSUM-class}
\Rdversion{1.1}
\docType{class}
\alias{calhitprobCUSUM-class}

\title{Class \code{"calhitprobCUSUM"}}
\description{
Derives the threshold needed to give a certain false alarm property for
CUSUM charts.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("calhitprobCUSUM", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{gridpoints}:}{Object of class \code{"numeric"}  }
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
    \item{SPCq}{\code{signature(property = "calhitprobCUSUM")}: ... }
	 }
}
\author{
Axel Gandy
}

\seealso{
 \code{\linkS4class{hitprobCUSUM}}.
}
\examples{
X<-rnorm(100);
chart <- new("SPCCUSUMNormal",Delta=1)
SPCproperty(data=X,nrep=10,
            property=new("calhitprobCUSUM",chart=chart,target=0.05,nsteps=1000))
#increase the number of repetitions for real applications!
}
\keyword{classes}
