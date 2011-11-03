\name{ARLShew-class}
\Rdversion{1.1}
\docType{class}
\alias{ARLShew-class}

\title{Class \code{"ARLShew"}}
\description{
Computes average run lengths (ARLs) for Shewhart charts.
}
\section{Slots}{
  \describe{
    \item{\code{threshold}:}{Object of class \code{"numeric"} -
      threshold of the chart. }
    \item{\code{chart}:}{Object of class \code{"SPCchart"} - chart type
      for which properties are to be computed. }
    \item{\code{lowerconf}:}{Object of class \code{"logical"}  defaults
      to FALSE - \code{\link{SPCproperty}} will report an upper
      confidence interval. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{SPCARL}"}, directly.
Class \code{"\linkS4class{SPCproperty}"}, by class "SPCARL", distance 2.
}
\section{Methods}{
  \describe{
    \item{SPCq}{\code{signature(property = "ARLShew")}: ... }
	 }
}
\author{
Axel Gandy
}

\seealso{
  \code{\linkS4class{ARLCUSUM}},  \code{\linkS4class{SPCARL}}.
}
\examples{
showClass("ARLShew")
}
\keyword{classes}
