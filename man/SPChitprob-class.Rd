\name{SPChitprob-class}
\Rdversion{1.1}
\docType{class}
\alias{SPChitprob-class}


\title{Class \code{"SPChitprob"}}
\description{
Virtual base class for hitting probability computations.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{threshold}:}{Object of class \code{"numeric"}  }
    \item{\code{nsteps}:}{Object of class \code{"numeric"}  }
    \item{\code{chart}:}{Object of class \code{"SPCchart"}  }
    \item{\code{lowerconf}:}{Object of class \code{"logical"}  }
  }
}
\section{Extends}{
Class \code{"\linkS4class{SPCproperty}"}, directly.
}
\section{Methods}{
  \describe{
    \item{qtrafo}{\code{signature(property = "SPChitprob")}: ... }
    \item{SPCoutput}{\code{signature(property = "SPChitprob")}: ... }
	 }
}
\author{ Axel Gandy
}

\seealso{
  \code{\linkS4class{hitprobCUSUM}} ,
  \code{\linkS4class{hitprobShew}}.
}
\examples{
showClass("SPChitprob")
}
\keyword{classes}
