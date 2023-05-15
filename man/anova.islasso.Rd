\name{anova.islasso}
\alias{anova.islasso}
\alias{print.anova.islasso}

\title{General Linear Combination method for \code{islasso} objects}

\description{
General linear hypotheses and confidence intervals estimation for linear combinantions of the regression coefficients in  \code{islasso} fits
}

\usage{
\method{anova}{islasso}(object, A, b = NULL, ci, ...)
}
\arguments{
\item{object}{a fitted model object of class "islasso".}

\item{A}{matrix (or vector) giving linear combinations of coefficients by rows, or a character vector giving the hypothesis in symbolic form (see Details). }

\item{b}{right-hand-side vector for hypothesis, with as many entries as rows in the hypothesis matrix A; can be omitted, in which case it defaults to a vector of zeroes.}

\item{ci}{optionally, a two columns matrix of estimated confidence intervals for the estimated coefficients.}

\item{...}{not used.}
}

\details{
For the islasso regression model with coefficients \eqn{\beta}{beta}, the null hypothesis is \deqn{H_0:A\beta=b}{H_0: A \beta=b} where A and b are known matrix and vector. 

The hypothesis matrix A can be supplied as a numeric matrix (or vector), the rows of which specify linear combinations of the model coefficients, which are tested equal to the corresponding entries in the right-hand-side vector b, which defaults to a vector of zeroes.

Alternatively, the hypothesis can be specified symbolically as a character vector with one or more elements, each of which gives either a linear combination of coefficients, or a linear equation in the coefficients (i.e., with both a left and right side separated by an equals sign). Components of a linear expression or linear equation can consist of numeric constants, or numeric constants multiplying coefficient names (in which case the number precedes the coefficient, and may be separated from it by spaces or an asterisk); constants of 1 or -1 may be omitted. Spaces are always optional. Components are separated by plus or minus signs. Newlines or tabs in hypotheses will be treated as spaces. See the examples below.


}
\value{
An object of class "anova.islasso" which contains the estimates, the standard errors, the Wald statistics and corresponding p value of each linear combination and of the restriced model.
}
\author{
The main function of the same name was inspired by the R function previously implemented by Vito MR Muggeo.

Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>
}

\examples{
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

o <- islasso(y ~ . - 1, data = data.frame(y = y, X), 
             family = gaussian())
anova(o, A = diag(p), b = coef)
anova(o, A = c("X1 + X2 + X3 + X4 + X5 = -7.5"))
anova(o, A = c("X1 + X2 + X3 + X4 + X5 = 0"))
anova(o, A = c("X6 + X7 + X8 + X9 + X10"), b = 8.75)
anova(o, A = c("X6 + X7 + X8 + X9 + X10"), b = 0)
anova(o, A = c("X1 + X2 + X3 + X4 + X5 = -7.5",
               "X6 + X7 + X8 + X9 + X10 = 8.75"))
anova(o, A = c("X1 + X2 + X3 + X4 + X5",
               "X6 + X7 + X8 + X9 + X10"), b = c(-7.5, 8.75))
anova(o, A = c("X1 + X2 + X3 + X4 + X5",
               "X6 + X7 + X8 + X9 + X10"))

}
