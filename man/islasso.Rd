\name{islasso}
\alias{islasso}
\alias{print.islasso}
\alias{islasso.fit}
\alias{coef.islasso}
\alias{vcov.islasso}
\alias{deviance.islasso}
\alias{AIC.islasso}
\alias{fitted.islasso}
\alias{residuals.islasso}
\alias{logLik.islasso}
\alias{print.logLik.islasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Induced Smoothing Lasso Model
}
\description{
\code{islasso} is used to fit induced smoothing lasso model. Fit linear, logistic, poisson and gamma regression models.

}
\usage{
islasso(formula, family = gaussian, lambda, data, weights, subset,
        offset, unpenalized, contrasts = NULL, control = islasso.control())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
an object of class \dQuote{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
}
  \item{family}{
a description of the error distribution, \code{family=gaussian}, \code{family=binomial} and \code{family=quasibinomial}, \code{family=poisson} and \code{quasipoisson}, \code{family=Gamma} are implemented. \code{family=gaussian} is implemented with \code{identity} link, \code{family=binomial} is implemented with \code{logit} and \code{probit} links, \code{family=poisson} is implemented with \code{log} link and \code{family=Gamma} is implemented with \code{inverse}, \code{log} and \code{identity} links.
}
  \item{lambda}{
Value of lambda used to compute beta. If lambda is missing the optimal lambda is computed using \code{\link{cv.glmnet}}.
}
  \item{data}{
an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which \code{islasso} is called.
}
  \item{weights}{
observation weights. Default is 1 for each observation.
}
  \item{subset}{
an optional vector specifying a subset of observations to be used in the fitting process.
}
\item{offset}{
this can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be NULL or a numeric vector of length equal to the number of cases.
}
\item{unpenalized}{
a vector used to specify the unpenalized estimators; unpenalized can be a vector of integers or characters specifying the names of the predictors with unpenalized estimators.
}
\item{contrasts}{an optional list. See the contrasts.arg of \code{model.matrix.default}.}
  \item{control}{
a list of parameters for controlling the fitting process (see \code{islasso.control} for more details).
}
}
\details{
Models for \code{islasso} are specified symbolically. A typical model has the form \code{response ~ terms} where response is the response vector and terms is a series of terms which specifies a linear predictor for response. A terms specification of the form \code{first + second} indicates all the terms in first together with all the terms in second with duplicates removed. A specification of the form \code{first*second} indicates the set of terms obtained by taking the interactions of all terms in first with all terms in second.
}
\value{
%%  ~Describe the value returned
a list of

\item{coefficients}{a named vector of coefficients}
\item{se}{a named vector of standard errors}
\item{res}{the working residuals}
\item{fitted.values}{the fitted values}
\item{linear.predictors}{the linear predictors}
\item{rank}{the estimated degrees of freedom}
\item{family}{the family object used}
\item{deviance}{the family deviance}
\item{null.deviance}{the family null deviance}
\item{aic}{the Akaike Information Criterion}
\item{df.null}{the degrees of freedom of a null model}
\item{phi}{the estimated dispersion parameter}
\item{beta.unbias}{unbiased coefficients}
\item{se.unbias}{unbiased standard errors}
\item{internal}{internal elements}
\item{control}{the value of the control argument used}
\item{model}{if requested (the default), the model frame used.}
\item{terms}{the terms object used.}
\item{contrasts}{(only where relevant) the contrasts used.}
\item{call}{the matched call}
\item{formula}{the formula supplied}
}
\references{
Cilluffo, G, Sottile, G, S, La Grutta, S and Muggeo, VMR (2019). \emph{The Induced Smoothed lasso: A practical framework for hypothesis testing in high dimensional regression}.
Statistical Methods in Medical Research, DOI: 10.1177/0962280219842890.

Cilluffo, G, Fasola, S, La Grutta, S and Muggeo, VMR (2016). \emph{The Induced Smoothed Lasso}.
Proceedings 31st International Workshop of Statistical Modeling, Rennes, vol.1 69-74.
}
\author{
The main function of the same name was inspired by the R function previously implemented by Vito MR Muggeo.
Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>
}

\seealso{
\code{\link{islasso.fit}}, \code{\link{coef.islasso}}, \code{\link{summary.islasso}}, \code{\link{residuals.islasso}}, \code{\link{AIC.islasso}}, \code{\link{logLik.islasso}}, \code{\link{fitted.islasso}}, \code{\link{predict.islasso}} and \code{\link{deviance.islasso}} methods.
}
\examples{

set.seed(1)
n <- 100
p <- 100
p1 <- 10  #number of nonzero coefficients
coef.veri <- sort(round(c(seq(.5, 3, l=p1/2), seq(-1, -2, l=p1/2)), 2))
sigma <- 1

coef <- c(coef.veri, rep(0, p-p1))

X <- matrix(rnorm(n*p), n, p)
eta <- drop(X\%*\%coef)

##### gaussian ######
mu <- eta
y <- mu + rnorm(n, 0, sigma)

o <- islasso(y~-1+X, family=gaussian)
o
summary(o)
coef(o)
fitted(o)
predict(o, type="response")
plot(o)
residuals(o)
deviance(o)
AIC(o)
logLik(o)

\donttest{
# for the interaction
o <- islasso(y~-1+X[,1]*X[,2], family=gaussian)

##### binomial ######
coef <- c(c(1,1,1), rep(0, p-3))
X <- matrix(rnorm(n*p), n, p)
eta <- drop(cbind(1, X)\%*\%c(1, coef))
mu <- binomial()$linkinv(eta)
y <- rbinom(n, 100, mu)
y <- cbind(y, 100-y)

o <- islasso(y~X, family=binomial)

##### poisson ######
coef <- c(c(1,1,1), rep(0, p-3))
X <- matrix(rnorm(n*p), n, p)
eta <- drop(cbind(1, X)\%*\%c(1, coef))
mu <- poisson()$linkinv(eta)
y <- rpois(n, mu)

o <- islasso(y~X, family=poisson)

##### Gamma ######
coef <- c(c(1,1,1), rep(0, p-3))
X <- matrix(rnorm(n*p), n, p)
eta <- drop(cbind(1, X)\%*\%c(1, coef))
mu <- Gamma(link="log")$linkinv(eta)
shape <- 10
phi <- 1 / shape
y <- rgamma(n, scale = mu / shape, shape = shape)

o <- islasso(y~X, family=Gamma(link="log"))

}
}
\keyword{models}
\keyword{regression}
