\name{ARLCUSUM-class}
\Rdversion{1.1}
\docType{class}
\alias{ARLCUSUM-class}

\title{Class \code{"ARLCUSUM"}}
\description{
Class for computing the average in-control length  of CUSUM charts.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ARLCUSUM",
    chart, threshold,...)}.
}
\section{Slots}{
  \describe{
    \item{\code{gridpoints}:}{Object of class \code{"numeric"}: number
      of gridpoints to be used in the approximating Markov chain. }
    \item{\code{threshold}:}{Object of class \code{"numeric"} -
      threshold for running the chart. }
    \item{\code{chart}:}{Object of class \code{"SPCchart"} - chart on
      which the property should be computed. }
    \item{\code{lowerconf}:}{Object of class \code{"logical"} - defaults
      to FALSE - \code{\link{SPCproperty}} will report an upper
      confidence interval.}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{SPCARL}"}, directly.
  Class \code{"\linkS4class{SPCproperty}"}, by class "SPCARL", distance 2.
}
\section{Methods}{
  \describe{
    \item{SPCq}{\code{signature(property = "ARLCUSUM")}: ... }
  }
}
\author{
  Axel Gandy
}


\seealso{
  \code{\linkS4class{ARLShew}},  \code{\linkS4class{SPCARL}}.
}
\examples{
showClass("ARLCUSUM")
}
\keyword{classes}
