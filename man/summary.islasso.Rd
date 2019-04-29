\name{summary.islasso}
\alias{summary.islasso}
\alias{print.summary.islasso}
\title{summary method for islasso fitted objects}
\usage{
\method{summary}{islasso}(object, type = c("wald", "gradient"),
  trace = FALSE, pval = 1, ...)
}
\arguments{
\item{object}{fitted \code{"islasso"} object}

\item{type}{the statistics used to compute p-values. \code{type = "wald"} used the Wald statistic and \code{type = "gradient"} used the Gradient statistic}

\item{trace}{if TRUE and \code{type = "score"} iteration bar is printed}

\item{pval}{one is the default value and it means all variables are printed; if \code{pval = 0.05} only the variables
with a p-value lower or equal than 0.05 are printed}

\item{...}{not used}
}
\description{
summary method for islasso fitted objects
}
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