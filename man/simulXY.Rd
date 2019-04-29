\name{simulXY}
\alias{simulXY}
\alias{modelX}
\title{Simulate model matrix and response}
\usage{
simulXY(n, p, beta0, beta, family=gaussian(), prop=.10, a=-3, b=3, sigma=1, size=1, 
        rho=0, scale=TRUE, seed)
}
\arguments{
\item{n}{number of observations.}

\item{p}{number of covariates in the model matrix.}

\item{beta0}{if not missing is used as coefficient to the intercept.}

\item{beta}{if not missing is used as true vector of p coefficients in the linear predictor.}

\item{family}{a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function.}

\item{prop}{if beta is missing, this proportion is used to define the length of the non-null coefficients. The default is 10\% of p.}

\item{a}{if beta is missing, it sets the lower limit of the uniform distribution. The default is -3.}
\item{b}{if beta is missing, it sets the upper limit of the uniform distribution. The default is 3.}
\item{sigma}{if family is 'gaussian', it sets the standard deviation of the response vector. The default is 1.}
\item{size}{if family is 'binomial', it sets the number of trials to buildthe response vector. The default is 1.}
\item{rho}{correlation value to define the variance covariance matrix to build the model matrix, i.e., rho^|i-j| i,j = 1,...,p and i different from j. The default is 0.}
\item{scale}{to centers and/or scales the columns of the mdoel matrix. The default is TRUE}
\item{seed}{if it is not missing, it sets the random number generetor.}
}
\description{
Simulate model matrix and response from a specified distribution.
}
\examples{

n <- 100
p <- 100
beta0 <- 1
beta <- c(runif(10, -3, 3), rep(0, p-10))
seed <- 1234
dat <- simulXY(n = n, p = p, beta0 = beta0, beta = beta, family = gaussian(), seed=seed)

}