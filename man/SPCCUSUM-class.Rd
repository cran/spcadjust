\name{SPCCUSUM-class}
\Rdversion{1.1}
\docType{class}
\alias{SPCCUSUM-class}

\title{Class \code{"SPCCUSUM"}}
\description{
Virtual class from which CUSUM charts should descend.
}
\section{Objects from the Class}{
Virtual class - should not be instantiated.
}
\section{Extends}{
Class \code{"\linkS4class{SPCchart}"}, directly.
}
\section{Methods}{
  \describe{
    \item{runchart}{\code{signature(chart = "SPCCUSUM")}: evaluates the
      chart, relying on the method updates, which descendants must provide. }
	 }
}
\author{
Axel Gandy
}
\seealso{\code{\linkS4class{SPCCUSUMNormal}},\code{\linkS4class{SPCCUSUMNonpar}}
}
\keyword{classes}
