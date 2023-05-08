\name{predict.islasso.path}
\alias{predict.islasso.path}
\title{Prediction method for islasso.path fitted objects}
\usage{
\method{predict}{islasso.path}(object, newdata, type = c("link", "response",
  "coefficients", "class"), lambda, ...)
}
\arguments{
\item{object}{a fitted object of class \code{"islasso.path"}.}

\item{newdata}{optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used.}

\item{type}{the type of prediction required. The default is on the scale of the linear predictors; the alternative "response" is on the scale of the response variable.
Thus for a default binomial model the default predictions are of log-odds (probabilities on logit scale) and type = "response" gives the predicted probabilities. The \code{coefficients} option returns coefficients. Type "class" applies only to "binomial" models, and produces the class label.}

\item{lambda}{Value(s) of the penalty parameter lambda at which predictions are required. Default is the entire sequence used to create the model.}

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
\code{\link{islasso.path}}, \code{\link{islasso.path.fit}}, \code{\link{coef.islasso.path}}, \code{\link{residuals.islasso.path}}, \code{\link{GoF.islasso.path}}, \code{\link{logLik.islasso.path}}, \code{\link{fitted.islasso.path}}, \code{\link{summary.islasso.path}} and \code{\link{deviance.islasso.path}} methods.
}
\examples{
\dontrun{
 set.seed(1)
 n <- 100
 p <- 30
 p1 <- 10  #number of nonzero coefficients
 coef.veri <- sort(round(c(seq(.5, 3, l=p1/2), seq(-1, -2, l=p1/2)), 2))
 sigma <- 1

 coef <- c(coef.veri, rep(0, p-p1))

 X <- matrix(rnorm(n*p), n, p)
 mu <- drop(X\%*\%coef)
 y <- mu + rnorm(n, 0,sigma)

 o <- islasso.path(y ~ ., data = data.frame(y = y, X), 
                   family = gaussian())
 temp <- GoF.islasso.path(o)
 predict(o, type = "response", lambda = temp$lambda.min)
}
}
