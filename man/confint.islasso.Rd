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
\author{
Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>
}

\seealso{
\code{\link{islasso.fit}}, \code{\link{coef.islasso}}, \code{\link{summary.islasso}}, \code{\link{residuals.islasso}}, \code{\link{AIC.islasso}}, \code{\link{logLik.islasso}}, \code{\link{fitted.islasso}}, \code{\link{predict.islasso}} and \code{\link{deviance.islasso}} methods.
}
\examples{
\dontrun{
  set.seed(1)
  n <- 100
  p <- 100
  p1 <- 20  #number of nonzero coefficients
  coef.veri <- sort(round(c(seq(.5, 3, l=p1/2), seq(-1, -2, l=p1/2)), 2))
  sigma <- 1

  coef <- c(coef.veri, rep(0, p-p1))

  X <- matrix(rnorm(n*p), n, p)
  mu <- drop(X\%*\%coef)
  y <- mu + rnorm(n, 0,sigma)

  lambda <- 2
  o <- islasso(y ~ ., data = data.frame(y = y, X), 
               family = gaussian(), lambda = lambda)
  confint(o)
  plot(confint(o))
}
}