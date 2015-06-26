\name{qtrafo-methods}
\docType{methods}
\alias{qtrafo}
\alias{qtrafo-methods}
\alias{qtrafo,SPCARL-method}
\alias{qtrafo,SPCcalARL-method}
\alias{qtrafo,SPCcalhitprob-method}
\alias{qtrafo,SPChitprob-method}
\title{Methods for Function \code{qtrafo} in Package \pkg{spcadjust}}
\description{
Transformations used in the properties of charts.
}
\section{Methods}{
\describe{

\item{\code{signature(property = "SPCARL")}}{
Returns exp(x), where \code{x} is the second argument of the method.
}

\item{\code{signature(property = "SPCcalARL")}}{
Returns exp(x), where \code{x} is the second argument of the method.
}

\item{\code{signature(property = "SPCcalhitprob")}}{
Returns exp(x), where \code{x} is the second argument of the method.
}

\item{\code{signature(property = "SPChitprob")}}{
Returns exp(x)/(1+exp(x))
}
}}
\keyword{methods}
