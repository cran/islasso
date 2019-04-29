\name{islasso-package}
\alias{islasso-package}
\docType{package}
\title{
The Induced Smoothed lasso: A practical framework for hypothesis testing in high dimensional regression
}
\description{
This package implements an induced smoothed approach for hypothesis testing in lasso regression.
}
\details{
\tabular{ll}{
Package: \tab islasso\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2019-04-20\cr
License: \tab GPL-2\cr
}
\code{\link{islasso}} is used to fit a l1-penalized models in glm family. It can be used to make inference in lasso framework. This function calls another function \code{\link{islasso.fit}}islasso.fit to estimate coefficients and standard errors. Many auxiliary functions are implemented \code{\link{summary.islasso}},  \code{\link{predict.islasso}},  \code{\link{coef.islasso}}, \code{\link{fitted.islasso}},  \code{\link{logLik.islasso}}, \code{\link{AIC.islasso}}, \code{\link{deviance.islasso}},  \code{\link{residuals.islasso}}.
}
\author{
Gianluca Sottile

Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>
}
\references{
Cilluffo, G, Sottile, G, S, La Grutta, S and Muggeo, VMR (2019). \emph{The Induced Smoothed lasso: A practical framework for hypothesis testing in high dimensional regression}.
Statistical Methods in Medical Research, DOI: 10.1177/0962280219842890.

Cilluffo, G, Fasola, S, La Grutta, S and Muggeo, VMR (2016). \emph{The Induced Smoothed Lasso}.
Proceedings 31st International Workshop of Statistical Modeling, Rennes, vol.1 69-74.
}
\keyword{ package }
\examples{

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
o <- islasso(y~-1+X, family=gaussian, lambda=lambda)
o
summary(o)

}
