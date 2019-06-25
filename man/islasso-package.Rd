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
Version: \tab 1.1.0\cr
Date: \tab 2019-05-27\cr
License: \tab GPL-2\cr
}
\code{\link{islasso}} is used to fit generalized linear models with a L1-penalty on (some) regression coefficients. Along with point estimates, the main advantage is to return the full covariance matrix of estimate. The resulting standard errors can be used to make inference in the lasso framework. The main function is \code{\link{islasso}} and the correspoinding fitter function \code{\link{islasso.fit}}, and many auxiliary functions are implemented to summarize and visualize results: \code{\link{summary.islasso}},  \code{\link{predict.islasso}},  \code{\link{coef.islasso}}, \code{\link{fitted.islasso}},  \code{\link{logLik.islasso}}, \code{\link{AIC.islasso}}, \code{\link{deviance.islasso}},  \code{\link{residuals.islasso}}.
}
\author{
Gianluca Sottile based on some preliminary functions by Vito Muggeo.

Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>
}
\references{
Cilluffo, G, Sottile, G, S, La Grutta, S and Muggeo, VMR (2019). \emph{The Induced Smoothed lasso: A practical framework for hypothesis testing in high dimensional regression}.
Statistical Methods in Medical Research, DOI: 10.1177/0962280219842890.

Sottile, G, Cilluffo, G, Muggeo, VMR (2019). \emph{The R package islasso: estimation and
  hypothesis testing in lasso regression}. Technical Report on ResearchGate.
  doi:10.13140/RG.2.2.16360.11521.
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
