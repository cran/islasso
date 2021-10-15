\name{summary.islasso.path}
\alias{summary.islasso.path}
\alias{print.summary.islasso.path}
\title{summary method for islasso.path fitted objects}

\description{
summary method for islasso.path fitted objects
}

\usage{
\method{summary}{islasso.path}(object, pval = 1, use.t = FALSE, lambda, ...)
}

\arguments{
\item{object}{fitted \code{"islasso.path"} object}

\item{pval}{a threshold p-value value indicating which coefficients should be printed. If \code{pval = 0.10}, say, only the variables/coefficients with \eqn{p-value\le 0.10}{p-value \le 0.10} are printed. Possible unpenalized coefficients (including the intercept if in the model) are always printed, regardless of their p-value.}

\item{use.t}{if \code{TRUE}, the p-values are computed using the t-distribution with residual model degrees of freedom}

\item{lambda}{Value of the penalty parameter lambda at which summary are required.}

\item{...}{not used}
}
\author{
Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>
}
\seealso{
\code{\link{islasso.path}}, \code{\link{islasso.path.fit}}, \code{\link{coef.islasso.path}}, \code{\link{residuals.islasso.path}}, \code{\link{GoF.islasso.path}}, \code{\link{logLik.islasso.path}}, \code{\link{fitted.islasso.path}}, \code{\link{predict.islasso.path}} and \code{\link{deviance.islasso.path}} methods.
}
\examples{
\dontrun{
#continues example from ?islasso.path
summary(o, pval = .1, lambda = 5) #print just the "borderline" significant coefficients
}
}
