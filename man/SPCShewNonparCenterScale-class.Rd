\name{SPCShewNonparCenterScale-class}
\Rdversion{1.1}
\docType{class}
\alias{SPCShewNonparCenterScale-class}

\title{Class \code{"SPCShewNonparCenterScale"}}
\description{
Implements a Shewart chart, where the observations are centered and
scaled and past data is resampled nonparametrically.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SPCShewNonparCenterScale", ...)}.
}
\section{Extends}{
Class \code{"\linkS4class{SPCShew}"}, directly.
Class \code{"\linkS4class{SPCchart}"}, by class "SPCShew", distance 2.
}
\section{Methods}{
  \describe{
    \item{getcdfupdates}{\code{signature(chart = "SPCShewNonparCenterScale")}: ... }
    \item{Pofdata}{\code{signature(chart = "SPCShewNonparCenterScale")}: ... }
    \item{resample}{\code{signature(chart = "SPCShewNonparCenterScale")}: ... }
    \item{updates}{\code{signature(chart = "SPCShewNonparCenterScale")}: ... }
    \item{xiofP}{\code{signature(chart = "SPCShewNonparCenterScale")}: ... }
	 }
}
\author{
Axel Gandy
}


\seealso{
 \code{\linkS4class{SPCShewNormalCenterScale}}
}
\examples{
showClass("SPCShewNonparCenterScale")
}
\keyword{classes}
