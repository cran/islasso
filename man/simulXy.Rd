\name{simulXy}
\alias{simulXy}
\alias{modelX}
\title{Simulate model matrix and response}
\usage{
simulXy(n, p, interc = 0, beta, family = gaussian(), prop =
  0.1, lim.b = c(-3, 3), sigma = 1, size = 1, rho = 0,
  scale = TRUE, seed, X)
}
\arguments{
\item{n}{number of observations.}

\item{p}{total number of covariates in the model matrix.}

\item{interc}{the model intercept.}

\item{beta}{the vector of p coefficients in the linear predictor.}

\item{family}{a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. Only gaussian, binomial or poisson are allowed.}

\item{prop}{if \code{beta} is missing, \code{prop} represent the quote of non-null coefficients out of \code{p}. The default is \code{0.10 p}.}

\item{lim.b}{if \code{beta} is missing, the coefficients come from uniform variates in \code{lim.b}. The default is (-3,3).}
\item{sigma}{if family is 'gaussian', the standard deviation of the response. The default is 1.}
\item{size}{if family is 'binomial', the number of trials to build the response vector. The default is 1.}
\item{rho}{correlation value to define the variance covariance matrix to build the model matrix, i.e., rho^|i-j| i,j = 1,...,p and i different from j. The default is 0.}
\item{scale}{Should the columns of the mdoel matrix be scaled? The default is TRUE}
\item{seed}{optional, the seed to generate the data.}
\item{X}{optional, the model matrix.}
}
\description{
Simulate model matrix and response from a specified distribution.
}
\examples{

n <- 100
p <- 100
beta <- c(runif(10, -3, 3), rep(0, p-10))
dat <- simulXy(n, p,  beta = beta, seed=1234)
}