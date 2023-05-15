\name{summary.islasso}
\alias{summary.islasso}
\alias{print.summary.islasso}
\title{summary method for islasso fitted objects}

\description{
summary method for islasso fitted objects
}

\usage{
\method{summary}{islasso}(object, pval = 1, which, use.t = FALSE,
  type.pval = "wald", ...)
}

\arguments{
\item{object}{fitted \code{"islasso"} object}

\item{pval}{a threshold p-value value indicating which coefficients should be printed. If \code{pval = 0.10}, say, only the variables/coefficients with \eqn{p-value\le 0.10}{p-value \le 0.10} are printed. Possible unpenalized coefficients (including the intercept if in the model) are always printed, regardless of their p-value.}

\item{which}{a specification of which parameters are to be given p-values. If missing, all parameters are considered.}

\item{use.t}{if \code{TRUE}, the p-values are computed using the t-distribution with residual model degrees of freedom}

\item{type.pval}{Only Wald-type confidence intervals are implemented yet! type.pval = "wald" (default) estimates and standard errors are used to build confidence interval}

\item{...}{not used}
}
\author{
Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>
}

\seealso{
\code{\link{islasso.fit}}, \code{\link{summary.islasso}}, \code{\link{residuals.islasso}}, \code{\link{logLik.islasso}}, \code{\link{predict.islasso}} and \code{\link{deviance.islasso}} methods.
}
\examples{
\dontrun{
#continues example from ?islasso
summary(o, pval = .1) #print just the "borderline" significant coefficients
}
}
