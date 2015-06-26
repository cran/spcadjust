\name{SPCCUSUMNonparCenterScale-class}
\Rdversion{1.1}
\docType{class}
\alias{SPCCUSUMNonparCenterScale-class}

\title{Class \code{"SPCCUSUMNonparCenterScale"}}
\description{
CUSUM chart with centering and scaling of the observations. Resampling
of the data is done none-parametrically.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SPCCUSUMNonparCenterScale", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{Delta}:}{Amount subtracted at each step in the control chart. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{SPCCUSUMNonpar}"}, directly.
Class \code{"\linkS4class{SPCCUSUM}"}, by class "SPCCUSUMNonpar", distance 2.
Class \code{"\linkS4class{SPCchart}"}, by class "SPCCUSUMNonpar", distance 3.
}
\section{Methods}{
  \describe{
    \item{updates}{\code{signature(chart = "SPCCUSUMNonparCenterScale")}: ... }
    \item{xiofP}{\code{signature(chart = "SPCCUSUMNonparCenterScale")}: ... }
	 }
}
\author{Axel Gandy}
\seealso{
 \code{\linkS4class{SPCCUSUMNonparCenterScale}}.
}
\examples{
chartnp <- new("SPCCUSUMNonparCenterScale",Delta=1)
X<- rnorm(100)
SPCproperty(data=X,
            nrep=100,property=new("ARLCUSUM",chart=chartnp,threshold=4))

}
\keyword{classes}
