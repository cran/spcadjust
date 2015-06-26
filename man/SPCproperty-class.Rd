\name{SPCproperty-class}
\Rdversion{1.1}
\docType{class}
\alias{SPCproperty-class}

\title{Class \code{"SPCproperty"}}
\description{
Virtual class from which all properties of control charts should be derived.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{chart}:}{Object of class \code{"SPCchart"}; the chart
      for which the property is to be computed. }
    \item{\code{lowerconf}:}{TRUE if a lower confindence interval should
      be reported by default, false otherwise. }
  }
}
\section{Methods}{
No methods defined with class "SPCproperty" in the signature.
}
\author{
Axel Gandy
}
\examples{
#list all descendants.
showClass("SPCproperty")
}
\keyword{classes}
