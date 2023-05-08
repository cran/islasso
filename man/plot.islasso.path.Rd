\name{plot.islasso.path}
\alias{plot.islasso.path}
\title{plot coefficient profile from a fitted "islasso.path" object.}
\description{
Produces a coefficient profile plot of the coefficient paths for a fitted "islasso.path" object.
}
\usage{
\method{plot}{islasso.path}(x, 
  yvar = c("coefficients", "se", "gradient", "weight", "gof"), 
  gof = c("none", "AIC", "BIC", "AICc", "eBIC", "GCV", "GIC"), 
  label = FALSE, legend = FALSE, ...)
}
\arguments{
\item{x}{an object of class \code{islasso}, usually, a result of a call to \code{\link{islasso.path}}.}
\item{yvar}{What is on the Y-axis. "coef" plot the log-lambda sequence against the coefficients; "se" plot the log-lambda sequence against the standard deviations; "gradient" plot the log-lambda sequence against the gradient; "weight" plot the log-lambda sequence against the mixture weight of the islasso method; "gof" plot the log-lambda sequence against the chosen criterion.}
\item{gof}{the chosen criterion to highlight the active variables. "none" doesn't highlight active variables.}
\item{label}{a logical flag indicating if some labels have to be added.}
\item{legend}{a logical flag indicating if legend has to be shown.}
\item{...}{other graphical parameters for the plot, i.e., main, xlab, ylab, xlim, ylim, lty, col, lwd, cex.axis, cex.lab, cex.main, gof_lty, gof_col and gof_lwd. The last three parameters are used to modify aspects of the legend, and of the goodness of fit measure used.}
}

\details{
A coefficient profile plot is produced for Induced Smoothing Lasso Model path.
}
\examples{
\dontrun{
  set.seed(1)
  n <- 100
  p <- 30
  p1 <- 10  #number of nonzero coefficients
  coef.veri <- sort(round(c(seq(.5, 3, l=p1/2), seq(-1, -2, l=p1/2)), 2))
  sigma <- 1

  coef <- c(coef.veri, rep(0, p-p1))

  X <- matrix(rnorm(n*p), n, p)
  mu <- drop(X\%*\%coef)
  y <- mu + rnorm(n, 0,sigma)

  o <- islasso.path(y ~ ., data = data.frame(y = y, X), 
                    family = gaussian())
  par(mfrow = c(2, 2))
  plot(o, yvar = "coefficients", gof = "AICc", label = TRUE)
  plot(o, yvar = "se", gof = "AICc")
  plot(o, yvar = "gradient", gof = "AICc")
  plot(o, yvar = "gof", gof = "AICc")
}
}
