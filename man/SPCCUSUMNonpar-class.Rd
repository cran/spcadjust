\name{SPCCUSUMNonpar-class}
\Rdversion{1.1}
\docType{class}
\alias{SPCCUSUMNonpar-class}

\title{Class \code{"SPCCUSUMNonpar"}}
\description{
Virtual class that provides basic functionality for nonparametric CUSUM
charts. To use this, create a descendant of the class and implement
\code{\link{xiofP-methods}} and \code{\link{updates-methods}}.
}
\section{Objects from the Class}{
Virtual class - no objects can be created.
}
\section{Extends}{
Class \code{"\linkS4class{SPCCUSUM}"}, directly.
Class \code{"\linkS4class{SPCchart}"}, by class "SPCCUSUM", distance 2.
}
\section{Methods}{
  \describe{
    \item{getcdfupdates}{\code{signature(chart = "SPCCUSUMNonpar")}: ... }
    \item{Pofdata}{\code{signature(chart = "SPCCUSUMNonpar")}: ... }
    \item{resample}{\code{signature(chart = "SPCCUSUMNonpar")}: ... }
	 }
}
\author{
Axel Gandy
}

\seealso{
 \code{\linkS4class{SPCCUSUMNonparCenterScale}} for an example of an
 instantiation of this class.
}
\examples{
showClass("SPCCUSUMNonpar")
}
\keyword{classes}
