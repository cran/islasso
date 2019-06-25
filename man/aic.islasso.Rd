\name{aic.islasso}
\alias{aic.islasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimization for the selection of the tuning parameter
}
\description{
This function performs a minimization of the AIC/BIC criterion for selecting the tuning parameter in \dQuote{\code{islasso}}.
}
\usage{
aic.islasso(object, method = c("aic", "bic"), interval, y, X, 
  intercept = FALSE, family = gaussian(), alpha = 1, offset, 
  weights, unpenalized, control = is.control())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{a fitted model object of class "islasso".}
\item{method}{the criterion to optimize, AIC or BIC.}
\item{interval}{the  lower and upper limits of \eqn{\lambda}{\lambda} wherein the AIC/BIC criterion should be optimized. Can be missing, if \code{object} has been obtained via cross-validation (and therefore includes the range of lambdas)  }
\item{y}{if \code{object} is missing, the response vector of length n.}
\item{X}{if \code{object} is missing, the design matrix of dimension n * p.}
\item{intercept}{if \code{object} is missing, if TRUE the intercept is added to the model matrix.}
\item{family}{if \code{object} is missing, a description of the error distribution, \code{family=gaussian}, \code{family=binomial} and \code{family=poisson} are implemented with canonical link.}
\item{alpha}{The elasticnet mixing parameter, with \eqn{0\le\alpha\le 1}. The penalty is defined
    as \deqn{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1.} \code{alpha=1}
    is the lasso penalty, and \code{alpha=0} the ridge penalty.}
\item{offset}{this can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be NULL or a numeric vector of length equal to the number of cases.}
\item{weights}{observation weights. Default is 1 for each observation.}
\item{unpenalized}{a vector used to specify the unpenalized estimators; unpenalized has to be a vector of logicals.}
\item{control}{a list of parameters for controlling the fitting process (see \code{islasso.control} for more details).}
}
\details{
Minimization of the Akaike Information Criterion (AIC) or Bayesian Information Criterion (BIC) is sometimes employed to select the tuning parameter as an alternative to the cross validation. The model degrees of freedom (not necessarly integers as in the plain lasso) used in AIC/BIC are computed as trace of the hat matrix at convergence. 
}
\value{
the optimal lambda value is returned
}

\author{
Gianluca Sottile

Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>

}

\seealso{
\code{\link{islasso.fit}}, \code{\link{coef.islasso}}, \code{\link{summary.islasso}}, \code{\link{residuals.islasso}}, \code{\link{AIC.islasso}}, \code{\link{logLik.islasso}}, \code{\link{fitted.islasso}}, \code{\link{predict.islasso}} and \code{\link{deviance.islasso}} methods.
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
y <- mu + rnorm(n, 0, sigma)

o <- islasso(y ~ X)

\dontrun{
#use the evaluation interval of the fit
lambda_aic <- aic.islasso(o, method="aic") 

#overwrites the evaluation interval for lambda
lambda_bic <- aic.islasso(o, interval=c(.1, 30), method="bic") 
}
}
\keyword{models}
\keyword{regression}

