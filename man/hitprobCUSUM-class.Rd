\name{hitprobCUSUM-class}
\Rdversion{1.1}
\docType{class}
\alias{hitprobCUSUM-class}

\title{Class \code{"hitprobCUSUM"}}
\description{
Computes hitting probabilities for CUSUM charts.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("hitprobCUSUM", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{gridpoints}:}{Object of class \code{"numeric"}  }
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
    \item{SPCq}{\code{signature(property = "hitprobCUSUM")}: ... }
	 }
}
\author{
Axel Gandy
}


\seealso{
 \code{\linkS4class{calhitprobCUSUM}}.
}
\examples{
X <-  rnorm(1000)
chart <- new("SPCCUSUMNormal",Delta=1)
SPCproperty(data=X,nrep=10,
            property=new("calhitprobCUSUM",chart=chart,target=0.05,nsteps=1000))
#for real applications increase nrep, the number of bootstrap replications.
}
\keyword{classes}
