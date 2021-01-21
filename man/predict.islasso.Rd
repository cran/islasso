\name{predict.islasso}
\alias{predict.islasso}
\title{Prediction method for islasso fitted objects}
\usage{
\method{predict}{islasso}(object, newdata, type = c("link", "response",
  "coefficients", "class"), se.fit = FALSE, ci = NULL, level = .95, ...)
}
\arguments{
\item{object}{a fitted object of class "islasso".}

\item{newdata}{optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used.}

\item{type}{the type of prediction required. The default is on the scale of the linear predictors; the alternative "response" is on the scale of the response variable.
Thus for a default binomial model the default predictions are of log-odds (probabilities on logit scale) and type = "response" gives the predicted probabilities. The \code{coefficients} option returns coefficients. Type "class" applies only to "binomial" models, and produces the class label.}

\item{se.fit}{logical switch indicating if confidence intervals are required.}

\item{ci}{optionally, a two columns matrix of estimated confidence intervals for the estimated coefficients.}

\item{level}{the confidence level required.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
An object depending on the type argument
}
\description{
Prediction method for islasso fitted objects
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
 o <- islasso(y~X, family=gaussian, lambda=lambda)
 predict(o)
}
}
