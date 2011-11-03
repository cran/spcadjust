\name{SPCchart-class}
\Rdversion{1.1}
\docType{class}
\alias{SPCchart-class}

\title{Class \code{"SPCchart"}}
\description{
  Virtual class from which all control charts are being derived.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Methods}{
  \describe{
    \item{xiofdata}{\code{signature(chart = "SPCchart")}: method that
      concatenates PofData and xiofdata, which descendants must provide.}
  }
}
\author{Axel Gandy}

\seealso{\code{\linkS4class{SPCCUSUM}},\code{\linkS4class{SPCShew}}
}
\keyword{classes}
