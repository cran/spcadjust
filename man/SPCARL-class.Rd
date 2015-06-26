\name{SPCARL-class}
\Rdversion{1.1}
\docType{class}
\alias{SPCARL-class}

\title{Class \code{"SPCARL"}}
\description{
Virtual class for average runlength (ARL) computations of CUSUM control charts.
}
\section{Objects from the Class}{
Virtual class - no objects can be created from it.
}
\section{Slots}{
  \describe{
    \item{\code{threshold}:}{ threshold of the chart, a positive double. }
    \item{\code{chart}:}{chart for which the ARL should be
      computed. Should be a descendant of \code{\linkS4class{SPCchart}}.}
    \item{\code{lowerconf}:}{Set to \code{FALSE}. An upper confidence
      interval will be reported by  \code{\link{SPCproperty}}. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{SPCproperty}"}, directly.
}
\section{Methods}{
No methods defined with class "SPCARL" in the signature.
}
\author{
Axel Gandy
}

\seealso{
 \code{\linkS4class{ARLCUSUM}}, \code{\linkS4class{ARLShew}} for implementations.
}
\examples{
showClass("SPCARL")
}
\keyword{classes}
