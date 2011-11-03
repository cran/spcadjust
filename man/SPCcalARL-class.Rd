\name{SPCcalARL-class}
\Rdversion{1.1}
\docType{class}
\alias{SPCcalARL-class}

\title{Class \code{"SPCcalARL"}}
\description{
Virtual class for calibrations to a given ARL.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{target}:}{Object of class \code{"numeric"}  }
    \item{\code{chart}:}{Object of class \code{"SPCchart"}  }
    \item{\code{lowerconf}:}{Object of class \code{"logical"}  }
  }
}
\section{Extends}{
Class \code{"\linkS4class{SPCproperty}"}, directly.
}
\section{Methods}{
  \describe{
    \item{qtrafo}{\code{signature(property = "SPCcalARL")}: ... }
    \item{SPCoutput}{\code{signature(property = "SPCcalARL")}: ... }
	 }
}
\author{
Axel Gandy
}
\seealso{
  \code{\linkS4class{calARLCUSUM}},
  \code{\linkS4class{calARLShew}}.
}
\examples{
showClass("SPCcalARL")
}
\keyword{classes}
