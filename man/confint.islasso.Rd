\name{confint.islasso}
\alias{confint.islasso}
\alias{print.confint.islasso}
\alias{plot.confint.islasso}
\title{confint method for \code{islasso} objects}
\usage{
\method{confint}{islasso}(object, parm, level = 0.95, ...)
}
\arguments{
\item{object}{a fitted model object of class "islasso".}

\item{parm}{a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names.
If missing, all parameters are considered.}

\item{level}{the confidence level required.}

\item{...}{additional argument(s) for methods.}
}
\description{
confint method for \code{islasso} objects
}
