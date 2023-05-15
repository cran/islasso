\name{GoF.islasso.path}
\alias{GoF.islasso.path}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimization for the selection of the tuning parameter
}
\description{
This function extracts the value of the tuning parameter which minimizes the AIC/BIC/AICc/eBIC/GCV/GIC criterion in \dQuote{\code{islasso.path}}.
}
\usage{
GoF.islasso.path(object, plot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{a fitted model object of class \code{"islasso.path".}}
\item{plot}{a logical flag indicating if each criterion have to be plotted}
\item{...}{further arguments passed to or from other methods.}
}
\details{
Minimization of the Akaike Information Criterion (AIC), or Bayesian Information Criterion (BIC) or several other criteria are sometimes employed to select the tuning parameter as an alternative to the cross validation. The model degrees of freedom (not necessarly integers as in the plain lasso) used in all methods are computed as trace of the hat matrix at convergence. 
}
\value{
A list of

\item{gof}{the goodness of fit measures}
\item{minimum}{the position of the optimal lambda values}
\item{lambda.min}{the optimal lambda values}
}

\author{
Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>
}

\seealso{
\code{\link{islasso.path}}, \code{\link{islasso.path.fit}}, \code{\link{coef.islasso.path}}, \code{\link{residuals.islasso.path}}, \code{\link{summary.islasso.path}}, \code{\link{logLik.islasso.path}}, \code{\link{fitted.islasso.path}}, \code{\link{predict.islasso.path}} and \code{\link{deviance.islasso.path}} methods.
}

\examples{
set.seed(1)
n <- 100
p <- 30
p1 <- 10  #number of nonzero coefficients
coef.veri <- sort(round(c(seq(.5, 3, l=p1/2), seq(-1, -2, l=p1/2)), 2))
sigma <- 1

coef <- c(coef.veri, rep(0, p-p1))

X <- matrix(rnorm(n*p), n, p)
mu <- drop(X\%*\%coef)
y <- mu + rnorm(n, 0, sigma)

o <- islasso.path(y ~ ., data = data.frame(y = y, X))
GoF.islasso.path(o)
}
\keyword{models}
\keyword{regression}

