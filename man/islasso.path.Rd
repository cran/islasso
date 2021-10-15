\name{islasso.path}
\alias{islasso.path}
\alias{print.islasso.path}
\alias{islasso.path.fit}
\alias{coef.islasso.path}
\alias{deviance.islasso.path}
\alias{fitted.islasso.path}
\alias{residuals.islasso.path}
\alias{logLik.islasso.path}
\alias{print.logLik.islasso.path}
\alias{model.matrix.islasso.path}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Induced Smoothed lasso path
}
\description{
\code{islasso.path} is used to fit a generalized linear model via induced smoothed lasso method. The regularization path is computed for the lasso or elasticnet penalty at a grid of values for the regularization parameter lambda. Fits linear, logistic, poisson and gamma regression models.

}
\usage{
islasso.path(formula, family = gaussian(), lambda = NULL, nlambda = 100, 
        lambda.min.ratio = ifelse(nobs < nvars, 1E-3, 1E-05), alpha = 1, data, 
        weights, subset, offset, contrasts = NULL, unpenalized, control = is.control())
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
A user supplied lambda sequence. Typical usage is to have the program compute its own lambda sequence based on nlambda and lambda.min.ratio. Supplying a value of lambda overrides this.
}
\item{nlambda}{
The number of lambda values - default is 100.
}
\item{lambda.min.ratio}{
Smallest value for lambda, as a fraction of lambda.max, the (data derived) entry value (i.e. the smallest value for which all coefficients are zero). The default depends on the sample size \code{nobs} relative to the number of variables \code{nvars}. If \code{nobs > nvars}, the default is 0.00001, close to zero. If \code{nobs < nvars}, the default is 0.001. A very small value of lambda.min.ratio will lead to a saturated fit in the \code{nobs < nvars} case. 
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
\item{contrasts}{an optional list. See the contrasts.arg of \code{model.matrix.default}.}
  \item{control}{
a list of parameters for controlling the fitting process (see \code{islasso.control} for more details).
}
\item{unpenalized}{
optional. A vector of integers or characters indicating any covariate (in the formula) with coefficients not to be penalized. The intercept, if included in the model, is \emph{always} unpenalized.
}

}
\details{
The sequence of models implied by lambda is fit the islasso method. \code{islasso} estimates regression models by imposing a lasso-type penalty on some or all regression coefficients. However the nonsmooth \eqn{L_1}{L1} norm penalty is replaced by a smooth approximation justified under the induced smoothing paradigm. The advantage is that reliable standard errors are returned as model output and hypothesis testing on linear combinantions of the regression parameters can be carried out straightforwardly via the Wald statistic. Simulation studies provide evidence that the proposed approach   controls type-I errors and exhibits good power in different scenarios.}
\value{
%%  ~Describe the value returned
A list of

\item{call}{the matched call.}
\item{Info}{a named matrix containing information about lambda values, estimated degrees of freedom, estimated dispersion parameters, deviance, loglikelhood, number of iterations and convergence criteria.}
\item{GoF}{a named matrix containing information criteria, i.e., AIC, BIC, AICc, eBIC, GCV, GIC.}
\item{Coef}{a \code{length(lambda) x nvars} matrix of coefficients.}
\item{SE}{a \code{length(lambda) x nvars} matrix of standard errors.}
\item{Weights}{a \code{length(lambda) x nvars} matrix of the weight of the mixture in the induced smoothed lasso.}
\item{Linear.predictors}{a \code{length(lambda) x nvars} matrix of linear predictors}
\item{Fitted.values}{a \code{length(lambda) x nvars} matrix of fitted values}
\item{Residuals}{a \code{length(lambda) x nvars} matrix of working residuals}
\item{Input}{a named list containing several input arguments, i.e., the numbers of observations and predictors, if an intercept ha to be estimated, the model matrix and the response vector, the observation weights, the offset, the family object used, The elasticnet mixing parameter and the vector used to specify the unpenalized estimators.}
\item{control}{the value of the control argument used.}
\item{formula}{the formula supplied.}
\item{model}{if requested (the default), the model frame used.}
\item{terms}{the terms object used.}
\item{data}{he data argument.}
\item{xlevels}{(where relevant) a record of the levels of the factors used in fitting.}
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
Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>
}

\seealso{
\code{\link{islasso.path.fit}}, \code{\link{coef.islasso.path}}, \code{\link{summary.islasso.path}}, \code{\link{residuals.islasso.path}}, \code{\link{GoF.islasso.path}}, \code{\link{logLik.islasso.path}}, \code{\link{fitted.islasso.path}}, \code{\link{predict.islasso.path}} and \code{\link{deviance.islasso.path}} methods.
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
eta <- drop(X\%*\%coef)

##### gaussian ######
mu <- eta
y <- mu + rnorm(n, 0, sigma)

o <- islasso.path(y ~ ., data = data.frame(y = y, X), 
             family = gaussian(), nlambda = 30L)
o
summary(o, lambda = 10)
coef(o, lambda = 10)
fitted(o, lambda = 10)
predict(o, type="response", lambda = 10)
plot(o, xvar = "coef")
residuals(o, lambda = 10)
deviance(o, lambda = 10)
logLik(o, lambda = 10)
GoF.islasso.path(o)

\dontrun{
##### binomial ######
coef <- c(c(1,1,1), rep(0, p-3))
X <- matrix(rnorm(n*p), n, p)
eta <- drop(cbind(1, X)\%*\%c(-1, coef))
mu <- binomial()$linkinv(eta)
y <- rbinom(n, 100, mu)
y <- cbind(y, 100-y)

o <- islasso.path(cbind(y1, y2) ~ ., 
             data = data.frame(y1 = y[,1], y2 = y[,2], X), 
             family = binomial(), nlambda = 30L)
temp <- GoF.islasso.path(o)
summary(o, pval = .05, lambda = temp$lambda.min["BIC"])

##### poisson ######
coef <- c(c(1,1,1), rep(0, p-3))
X <- matrix(rnorm(n*p), n, p)
eta <- drop(cbind(1, X)\%*\%c(1, coef))
mu <- poisson()$linkinv(eta)
y <- rpois(n, mu)

o <- islasso.path(y ~ ., data = data.frame(y = y, X), 
             family = poisson(), nlambda = 30L)
temp <- GoF.islasso.path(o)
summary(o, pval = .05, lambda = temp$lambda.min["BIC"])

##### Gamma ######
coef <- c(c(1,1,1), rep(0, p-3))
X <- matrix(rnorm(n*p), n, p)
eta <- drop(cbind(1, X)\%*\%c(-1, coef))
mu <- Gamma(link="log")$linkinv(eta)
shape <- 10
phi <- 1 / shape
y <- rgamma(n, scale = mu / shape, shape = shape)

o <- islasso.path(y ~ ., data = data.frame(y = y, X), 
             family = Gamma(link = "log"), nlambda = 30L)
temp <- GoF.islasso.path(o)
summary(o, pval = .05, lambda = temp$lambda.min["BIC"])
}
}
\keyword{models}
\keyword{regression}
