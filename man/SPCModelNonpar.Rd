% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{SPCModelNonpar}
\alias{SPCModelNonpar}
\title{Generic Model for Nonparametric Resampling}
\usage{
SPCModelNonpar(updates, xiofP)
}
\arguments{
\item{updates}{function that computes updates.}

\item{xiofP}{function that computes xi given P.}
}
\value{
An object of class SPCDataModel.
}
\description{
Generic model that allows nonparametric resampling (with
replacement) of the data. The transformation of data into updates
needs to be defined via the arguments.
}
\details{
The parameters to the functions being returned have the following meaning.
\itemize{
\item data: a numeric vector or a matrix where the rows contain the observations.
\item xi: depends on the parameter xiofP.
\item P: The \code{data} with no modifications (thus either a numeric vector or a matrix).
}

The main operations are defined as follows:
\itemize{
\item resample(P): generates a new data set of the same size by either resampling the values (if the data is a vector) or by resampling the rows of the data matrix (both use resampling with replacement).
}
}

