\name{plot.islasso.path}
\alias{plot.islasso.path}
\title{plot coefficient profile from a fitted "islasso.path" object.}
\description{
Produces a coefficient profile plot of the coefficient paths for a fitted "islasso.path" object.
}
\usage{
\method{plot}{islasso.path}(x, 
  yvar = c("coefficients", "se", "gradient", "weight", 
    "AIC", "BIC", "AICc", "eBIC", "GCV", "GIC"), 
  label = FALSE, cex.lab = 1, ...)
}
\arguments{
\item{x}{an object of class \code{islasso}, usually, a result of a call to \code{\link{islasso.path}}.}
\item{yvar}{What is on the Y-axis. "coef" plot the log-lambda sequence against the coefficients; "se" plot the log-lambda sequence against the standard deviations; "gradient" plot the log-lambda sequence against the gradient; "weight" plot the log-lambda sequence against the mixture weight of the islasso method; "AIC", "BIC", "AICc", "eBIC", "GCV" and "GIC" plot the log-lambda sequence against the chosen criterion.}
\item{label}{a logical flag indicating if some labels have to be added on each plot.}
\item{cex.lab}{a numerical scalar giving the amount by which plotting labels should be scaled.}
\item{...}{other graphical parameters for the plot}
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
  plot(o, yvar = "coefficients")
  plot(o, yvar = "se")
  plot(o, yvar = "gradient")
  plot(o, yvar = "weight")
  
  par(mfrow = c(2, 3))
  plot(o, yvar = "AIC")
  plot(o, yvar = "BIC")
  plot(o, yvar = "AICc")
  plot(o, yvar = "eBIC")
  plot(o, yvar = "GCV")
  plot(o, yvar = "GIC")
}
}
