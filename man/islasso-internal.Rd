\name{islasso-internal}
\title{Internal Functions}

\alias{format.perc}
\alias{qqNorm}
\alias{ginv2}
\alias{printCoefmat2}
\description{Internal islasso functions.}
\usage{
format.perc(probs, digits)
qqNorm(x, probs=seq(.005, .995, l=200), centre=FALSE, scale=FALSE, 
  leg=TRUE, mean=0, sd=1, add=FALSE, p.col=1, dF=FALSE, ylab, ...)
ginv2(X, tol = sqrt(.Machine$double.eps))
printCoefmat2(x, digits = max(3L, getOption("digits") - 2L), 
  signif.stars = getOption("show.signif.stars"), 
  signif.legend = signif.stars, dig.tst = max(1L, min(5L, digits - 1L)), 
  cs.ind = 1:k, tst.ind = k + 1, zap.ind = integer(), P.values = NULL, 
  has.Pvalue = nc >= 4L && length(cn <- colnames(x)) && 
    substr(cn[nc], 1L, 3L) \%in\% c("Pr(", "p-v"), 
  eps.Pvalue = .Machine$double.eps, na.print = "NA", quote = FALSE, 
  right = TRUE, ...)
}
\author{Gianluca Sottile (\email{gianluca.sottile@unipa.it})}
\details{These functions are not intended for users.}
\keyword{internal}
