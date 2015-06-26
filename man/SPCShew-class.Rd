\name{SPCShew-class}
\Rdversion{1.1}
\docType{class}
\alias{SPCShew-class}

\title{Class \code{"SPCShew"}}
\description{
Virtual Class for Shewhart charts.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{twosided}:}{TRUE if a two-sided chart should be used. Default:FALSE}
  }
}
\section{Extends}{
Class \code{"\linkS4class{SPCchart}"}, directly.
}
\section{Methods}{
  \describe{
    \item{runchart}{\code{signature(chart = "SPCShew")}: ... }
	 }
}
\author{
  Axel Gandy
}


\seealso{
 \code{\linkS4class{SPCShewNormalCenterScale}},
 \code{\linkS4class{SPCShewNonparCenterScale}}.
}
\examples{
showClass("SPCShew")
}
\keyword{classes}
