\name{optim.islasso}
\alias{optim.islasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimization for the selection of the tuning parameter
}
\description{
This function performs a minimization of the AIC/BIC criterion for the model \dQuote{\code{islasso}} for the optimal choice of the tuning parameter.
}
\usage{
optim.islasso(object, y, X, intercept = FALSE, family = gaussian(), offset, 
  weights, unpenalized, control = islasso.control(), interval, method = c("aic", "bic"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{a fitted model object of class "islasso".}
\item{y}{vector of observations of length n.}
\item{X}{design matrix of dimension n * p.}
\item{intercept}{if TRUE the intercept is added to the model matrix.}
\item{family}{a description of the error distribution, \code{family=gaussian}, \code{family=binomial} and \code{family=poisson} are implemented with canonical link.}
\item{offset}{this can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be NULL or a numeric vector of length equal to the number of cases.}
\item{weights}{observation weights. Default is 1 for each observation.}
\item{unpenalized}{a vector used to specify the unpenalized estimators; unpenalized has to be a vector of logicals.}
\item{control}{a list of parameters for controlling the fitting process (see \code{islasso.control} for more details).}
\item{interval}{a 2-vector of lower and upper bound in which lambda lie.}
\item{method}{criterion to optimize, if \code{method=}\dQuote{aic} the Akaike Information Criterion is optimized, if \code{method=}\dQuote{bic} the Bayesian Information Criterion is optimized.}
}
\details{
The optimal lambda choice is based on the optimization of a concave function of AIC or BIC.
%%  ~~ If necessary, more details than the description above ~~
}
\value{

only the optimal lambda value is returned
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

o <- islasso(y ~ X, family=gaussian())

\donttest{
lambda_aic <- optim.islasso(o, interval=c(.1, 30), method="aic"); lambda_aic
lambda_bic <- optim.islasso(o, interval=c(.1, 30), method="bic"); lambda_bic
}
}
\keyword{models}
\keyword{regression}
