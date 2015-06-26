\name{SPCCUSUMlm-class}
\Rdversion{1.1}
\docType{class}
\alias{SPCCUSUMlm-class}
\alias{updates,SPCCUSUMlm-method}
\alias{xiofP,SPCCUSUMlm-method}

\title{Class \code{"SPCCUSUMlm"}}
\description{
Risk-adjusted CUSUM chart based on a linear model.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SPCCUSUMlm", ...)}.
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
    \item{updates}{\code{signature(chart = "SPCCUSUMlm")}: ... }
    \item{xiofP}{\code{signature(chart = "SPCCUSUMlm")}: ... }
	 }
}
\author{
Axel Gandy
}
\examples{
n <- 1000
Xlinreg <- data.frame(x1= rbinom(n,1,0.4), x2= runif(n,0,1), x3= rnorm(n))
Xlinreg$y <- 2 + Xlinreg$x1 + Xlinreg$x2 + Xlinreg$x3 + rnorm(n)

chartlinreg <- new("SPCCUSUMlm",Delta=1,formula="y~x1+x2+x3")
#increase nrep in real applications
SPCproperty(data=Xlinreg,
            nrep=10,
            property=new("calARLCUSUM",chart=chartlinreg,target=100))
chartlinreg2 <- new("SPCCUSUMlm",Delta=1,formula="y~x1")
#increase nrep in real applications
SPCproperty(data=Xlinreg,
            nrep=10,
            property=new("calARLCUSUM",chart=chartlinreg2,target=100))
}
\keyword{classes}
