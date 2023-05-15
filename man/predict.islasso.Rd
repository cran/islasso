\name{predict.islasso}
\alias{predict.islasso}
\title{Prediction method for islasso fitted objects}
\usage{
\method{predict}{islasso}(object, newdata = NULL, 
  type = c("link", "response", "coefficients", "class", "terms"), 
  se.fit = FALSE, ci = NULL, type.ci = "wald", 
  level = .95, terms = NULL, na.action = na.pass, ...)
}
\arguments{
\item{object}{a fitted object of class "islasso".}

\item{newdata}{optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used.}

\item{type}{the type of prediction required. The default is on the scale of the linear predictors; the alternative "response" is on the scale of the response variable.
Thus for a default binomial model the default predictions are of log-odds (probabilities on logit scale) and type = "response" gives the predicted probabilities. The \code{coefficients} option returns coefficients. Type "class" applies only to "binomial" models, and produces the class label. The "terms" option returns a matrix giving the fitted values of each term in the model formula on the linear predictor scale.}

\item{se.fit}{logical switch indicating if confidence intervals are required.}

\item{ci}{optionally, a two columns matrix of estimated confidence intervals for the estimated coefficients.}

\item{type.ci}{Only Wald-type confidence intervals are implemented yet! type.ci = "wald" estimates and standard errors are used to build confidence interval}

\item{level}{the confidence level required.}

\item{terms}{with type = "terms" by default all terms are returned. A character vector specifies which terms are to be returned.}

\item{na.action}{function determining what should be done with missing values in newdata. The default is to predict NA.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
An object depending on the type argument
}
\description{
Prediction method for islasso fitted objects
}
\author{
Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>
}

\seealso{
\code{\link{islasso.fit}}, \code{\link{summary.islasso}}, \code{\link{residuals.islasso}}, \code{\link{logLik.islasso}}, \code{\link{predict.islasso}} and \code{\link{deviance.islasso}} methods.
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
 o <- islasso(y ~ ., data = data.frame(y = y, X), lambda = lambda)
 predict(o, type = "response")
}
