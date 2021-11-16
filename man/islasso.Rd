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
\alias{model.matrix.islasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Induced Smoothed lasso
}
\description{
\code{islasso} is used to fit lasso regression models wherein the nonsmooth \eqn{L_1}{L1} norm penalty is replaced by a smooth approximation justified under the induced smoothing paradigm. Simple lasso-type or elastic-net penalties are permitted and Linear, Logistic, Poisson and Gamma responses are allowed.

}
\usage{
islasso(formula, family = gaussian, lambda, alpha = 1, data, weights, subset,
        offset, unpenalized, contrasts = NULL, control = is.control())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
an object of class \dQuote{formula} (or one that can be coerced to that class): the `usual'  symbolic description of the model to be fitted.
}
  \item{family}{
the assumed response distribution. Gaussian, (quasi) Binomial, (quasi) Poisson, and Gamma are allowed. \code{family=gaussian} is implemented with \code{identity} link, \code{family=binomial} is implemented with \code{logit} or \code{probit} links, \code{family=poisson} is implemented with \code{log} link, and \code{family=Gamma} is implemented with \code{inverse}, \code{log} and \code{identity} links.
}
  \item{lambda}{
Value of the tuning parameter in the objective. If missing, the optimal lambda is computed using \code{\link{cv.glmnet}}.
}
  \item{alpha}{
The elastic-net mixing parameter, with \eqn{0\le\alpha\le 1}. The penalty is defined
    as \deqn{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1.} \code{alpha=1}
    is the lasso penalty, and \code{alpha=0} the ridge penalty.
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
optional. A vector of integers or characters indicating any covariate (in the formula) with coefficients not to be penalized. The intercept, if included in the model, is \emph{always} unpenalized.
}
\item{contrasts}{an optional list. See the contrasts.arg of \code{model.matrix.default}.}
  \item{control}{
a list of parameters for controlling the fitting process (see \code{islasso.control} for more details).
}
}
\details{
\code{islasso} estimates regression models by imposing a lasso-type penalty on some or all regression coefficients. However the nonsmooth \eqn{L_1}{L1} norm penalty is replaced by a smooth approximation justified under the induced smoothing paradigm. The advantage is that reliable standard errors are returned as model output and hypothesis testing on linear combinantions of the regression parameters can be carried out straightforwardly via the Wald statistic. Simulation studies provide evidence that the proposed approach   controls type-I errors and exhibits good power in different scenarios.}
\value{
%%  ~Describe the value returned
A list of

\item{coefficients}{a named vector of coefficients}
\item{se}{a named vector of standard errors}
\item{residuals}{the working residuals}
\item{fitted.values}{the fitted values}
\item{rank}{the estimated degrees of freedom}
\item{family}{the family object used}
\item{linear.predictors}{the linear predictors}
\item{deviance}{the family deviance}
\item{aic}{the Akaike Information Criterion}
\item{null.deviance}{the family null deviance}
\item{iter}{the number of iterations of IWLS used}
\item{weights}{the working weights, that is the weights in the final iteration of the IWLS fit}
\item{df.residual}{the residual degrees of freedom}
\item{df.null}{the degrees of freedom of a null model}
\item{converged}{logical. Was the IWLS algorithm judged to have converged?}
\item{model}{if requested (the default), the model frame used.}
\item{call}{the matched call}
\item{formula}{the formula supplied}
\item{terms}{the terms object used}
\item{data}{he data argument.}
\item{offset}{the offset vector used.}
\item{control}{the value of the control argument used}
\item{xlevels}{(where relevant) a record of the levels of the factors used in fitting.}
\item{lambda}{the lambda value used in the islasso algorithm}
\item{alpha}{the elasticnet mixing parameter}
\item{dispersion}{the estimated dispersion parameter}
\item{internal}{internal elements}
\item{contrasts}{(only where relevant) the contrasts used.}
}
\references{
Cilluffo, G, Sottile, G, S, La Grutta, S and Muggeo, VMR (2019). \emph{The Induced Smoothed lasso: A practical framework for hypothesis testing in high dimensional regression}.
Statistical Methods in Medical Research, DOI: 10.1177/0962280219842890.

Sottile, G, Cilluffo, G, Muggeo, VMR (2019). \emph{The R package islasso: estimation and
  hypothesis testing in lasso regression}. Technical Report on ResearchGate.
  doi:10.13140/RG.2.2.16360.11521.
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

o <- islasso(y ~ ., data = data.frame(y = y, X), 
             family = gaussian())
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

\dontrun{
confint(o, parm = 1:10)

# for the interaction
o <- islasso(y ~ X1 * X2, data = data.frame(y = y, X), 
             family = gaussian())

##### binomial ######
coef <- c(c(1,1,1), rep(0, p-3))
X <- matrix(rnorm(n*p), n, p)
eta <- drop(cbind(1, X)\%*\%c(-1, coef))
mu <- binomial()$linkinv(eta)
y <- rbinom(n, 100, mu)
y <- cbind(y, 100-y)

o <- islasso(cbind(y1, y2) ~ ., 
             data = data.frame(y1 = y[,1], y2 = y[,2], X), 
             family = binomial())
summary(o, pval = .05)

##### poisson ######
coef <- c(c(1,1,1), rep(0, p-3))
X <- matrix(rnorm(n*p), n, p)
eta <- drop(cbind(1, X)\%*\%c(1, coef))
mu <- poisson()$linkinv(eta)
y <- rpois(n, mu)

o <- islasso(y ~ ., data = data.frame(y = y, X), 
             family = poisson())
summary(o, pval = .05)

##### Gamma ######
coef <- c(c(1,1,1), rep(0, p-3))
X <- matrix(rnorm(n*p), n, p)
eta <- drop(cbind(1, X)\%*\%c(-1, coef))
mu <- Gamma(link="log")$linkinv(eta)
shape <- 10
phi <- 1 / shape
y <- rgamma(n, scale = mu / shape, shape = shape)

o <- islasso(y ~ ., data = data.frame(y = y, X), 
             family = Gamma(link = "log"))
summary(o, pval = .05)
}
}
\keyword{models}
\keyword{regression}
