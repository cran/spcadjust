\name{SPCCUSUMlogreg-class}
\Rdversion{1.1}
\docType{class}
\alias{SPCCUSUMlogreg-class}
\alias{updates,SPCCUSUMlogreg-method}
\alias{xiofP,SPCCUSUMlogreg-method}

\title{Class \code{"SPCCUSUMlogreg"}}
\description{
Risk-adjusted CUSUM chart for binary responses based on a logistic
regression model.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SPCCUSUMlogreg", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{Delta}:}{Object of class \code{"numeric"}  }
    \item{\code{formula}:}{Object of class \code{"character"}  }
  }
}
\section{Extends}{
Class \code{"\linkS4class{SPCCUSUMNonpar}"}, directly.
Class \code{"\linkS4class{SPCCUSUM}"}, by class "SPCCUSUMNonpar", distance 2.
Class \code{"\linkS4class{SPCchart}"}, by class "SPCCUSUMNonpar", distance 3.
}
\section{Methods}{
  \describe{
    \item{updates}{\code{signature(chart = "SPCCUSUMlogreg")}: ... }
    \item{xiofP}{\code{signature(chart = "SPCCUSUMlogreg")}: ... }
	 }
}
\author{
Axel Gandy
}
\examples{
n <- 1000
Xlogreg <- data.frame(x1=rbinom(n,1,0.4), x2=runif(n,0,1), x3=rnorm(n))
xbeta <- -1+Xlogreg$x1*100+Xlogreg$x2+Xlogreg$x3
Xlogreg$y <- rbinom(n,1,exp(xbeta)/(1+exp(xbeta)))
chartlogreg <- new("SPCCUSUMlogreg",Delta= 1, formula="y~x1+x2+x3")
SPCproperty(data=Xlogreg,
            nrep=10,
            property=new("calARLCUSUM",chart=chartlogreg,target=100))
#increase nrep for real applications.
}
\keyword{classes}
