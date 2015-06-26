\name{SPCcalhitprob-class}
\Rdversion{1.1}
\docType{class}
\alias{SPCcalhitprob-class}

\title{Class \code{"SPCcalhitprob"}}
\description{
Virtual base class for calibrating control charts to a certain hitting probability.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{target}:}{Object of class \code{"numeric"}  }
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
    \item{qtrafo}{\code{signature(property = "SPCcalhitprob")}: ... }
    \item{SPCoutput}{\code{signature(property = "SPCcalhitprob")}: ... }
	 }
}
\author{
Axel Gandy
}
\seealso{
  \code{\linkS4class{calhitprobCUSUM}},
  \code{\linkS4class{calhitprobShew}}.
}
\examples{
showClass("SPCcalhitprob")
}
\keyword{classes}
