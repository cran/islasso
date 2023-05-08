\name{plot.islasso}
\alias{plot.islasso}
\title{Diagnostics plots for Induced Smoothing Lasso Model}
\usage{
\method{plot}{islasso}(x, ...)
}
\arguments{
\item{x}{an object of class \code{islasso}, usually, a result of a call to \code{\link{islasso}}.}

\item{...}{other graphical parameters for the plot}
}
\description{
Diagnostics plots for Induced Smoothing Lasso Model
}
\author{
Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>
}

\seealso{
\code{\link{islasso.fit}}, \code{\link{summary.islasso}}, \code{\link{residuals.islasso}}, \code{\link{logLik.islasso}}, \code{\link{predict.islasso}} and \code{\link{deviance.islasso}} methods.
}
\details{
The plot on the top left is a plot of the standard deviance residuals against the fitted values.
The plot on the top right is a normal QQ plot of the standardized deviance residuals.
The red line is the expected line if the standardized residuals are normally distributed, i.e. it is the line with intercept 0 and slope 1.
The bottom two panels are plots of link and variance functions. On the left is squared standardized Pearson residuals against the fitted values.
On the right working vector against the linear predictor.
}
\examples{
\dontrun{
  set.seed(1)
  n <- 100
  p <- 100
  p1 <- 20  #number of nonzero coefficients
  coef.veri <- sort(round(c(seq(.5, 3, l=p1/2), seq(-1, -2, l=p1/2)), 2))
  sigma <- 1

  coef <- c(coef.veri, rep(0, p-p1))

  X <- matrix(rnorm(n*p), n, p)
  mu <- drop(X\%*\%coef)
  y <- mu + rnorm(n, 0,sigma)

  lambda <- 2
  o <- islasso(y ~ ., data = data.frame(y = y, X), 
               family = gaussian(), lambda = lambda)
  plot(o)
}
}
