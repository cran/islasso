\name{summary.islasso}
\alias{summary.islasso}
\alias{print.summary.islasso}
\title{summary method for islasso fitted objects}

\description{
summary method for islasso fitted objects
}

\usage{
\method{summary}{islasso}(object, pval = 1, use.t=FALSE, ...)
}

\arguments{
\item{object}{fitted \code{"islasso"} object}

\item{pval}{a threshold p-value value indicating which coefficients should be printed. if \code{pval = 0.10}, say, only the variables/coefficients
with \eqn{p-value\le 0.10}{p-value\le 0.10} are printed.}

\item{use.t}{if \code{TRUE}, the p-values are computed using the t-distribution with residual model degrees of freedom}

\item{...}{not used}
}
\examples{
\dontrun{
#continues example from ?islasso
summary(o, pval=.1) #print just the "borderline" significant coefficients
}
}