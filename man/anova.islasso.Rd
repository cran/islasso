\name{anova.islasso}
\alias{anova.islasso}
\alias{print.anova.islasso}

\title{General Linear Hypotheses method for \code{islasso} objects}

\description{
General linear hypotheses for linear combinantions of the regression coefficients in  \code{islasso} fits
}

\usage{
\method{anova}{islasso}(object, A, b, ...)
}
\arguments{
\item{object}{a fitted model object of class "islasso".}

\item{A}{a specification of the linear hypotheses to be tested. Linear functions can be specified by either a single vector (of length p) or by a matrix (of dimension k x p) of one or more linear hypotheses. }

\item{b}{an optional numeric vector specifying the right hand side of the hypothesis. Can be a scalar.}

\item{...}{not used.}
}

\details{
For the islasso regression model with coefficients \eqn{\beta}{beta}, the null hypothesis is \deqn{H_0:A\beta=b}{H_0: A \beta=b} where A and b are known matrix and vector. A can be a vector and b can be a scalar
}

\examples{
\dontrun{
set.seed(1)
n <- 100
p <- 100
p1 <- 10  #number of nonzero coefficients
coef.true <- sort(round(c(seq(.5, 3, l=p1/2), seq(-1, -2, l=p1/2)), 2))
sigma <- 1

coef <- c(coef.true, rep(0, p-p1))

X <- matrix(rnorm(n*p), n, p)
eta <- drop(X \%*\% coef)
mu <- eta
y <- mu + rnorm(n, 0, sigma)

o <- islasso(y~-1+X, family=gaussian)
A <- rbind(rep(c(1,0), c(10, p-10)),
           rep(c(0,1), c(10, p-10)))
anova(o, A)

A <- cbind(diag(10), matrix(0, 10, p-10))
b <- coef.true
anova(o, A, b)
}
}