\name{islasso-internal}
\title{Internal Functions}

\alias{fmt.perc}
\alias{qqNorm}
\alias{ginv2}
\alias{modelX}
\alias{printCoefmat2}

\alias{.checkinput}
\alias{.startpoint}
\alias{.islasso}

\alias{checkinput.islasso.path}
\alias{startpoint.islasso.path}
\alias{islasso.path.fit.glm}
\alias{islasso.path.fitter2}
\alias{interpolate}

\alias{makeHyp}
\alias{printHyp}

\alias{cislasso}
\alias{ci.fitted.islasso}

\alias{grad.islasso}
\alias{hess.islasso}
\alias{is.influence}
\alias{islasso.diag}
\alias{islasso.diag.plots}
\alias{predislasso}

\description{Internal islasso functions.}
\usage{
fmt.perc(probs, digits)
qqNorm(x, probs=seq(.005, .995, l=200), centre=FALSE, scale=FALSE, 
  leg=TRUE, mean=0, sd=1, add=FALSE, p.col=1, dF=FALSE, ylab, ...)
ginv2(X, tol = sqrt(.Machine$double.eps))
modelX(n, p, rho=.5, scale=TRUE)
printCoefmat2(x, digits = max(3L, getOption("digits") - 2L), 
  signif.stars = getOption("show.signif.stars"), 
  signif.legend = signif.stars, dig.tst = max(1L, min(5L, digits - 1L)), 
  cs.ind = 1:k, tst.ind = k + 1, zap.ind = integer(), P.values = NULL, 
  has.Pvalue = nc >= 4L && length(cn <- colnames(x)) && 
    substr(cn[nc], 1L, 3L) \%in\% c("Pr(", "p-v"), 
  eps.Pvalue = .Machine$double.eps, na.print = "NA", quote = FALSE, 
  right = TRUE, ...)
.checkinput(X, y, family, alpha, intercept, weights, offset, 
  unpenalized, control)
.startpoint(X, y, lambda, alpha, weights, offset, mustart, 
  family, intercept, setting)
.islasso(prep, start, Lambda, fam, link)

checkinput.islasso.path(X, y, family, lambda, nlambda, lambda.min.ratio, 
  alpha, intercept, weights, offset, unpenalized, control)
startpoint.islasso.path(X, y, lambda, alpha, weights, offset, mustart, 
  family, intercept, setting)
islasso.path.fit.glm(prep, start, lambda, fam, link)
islasso.path.fitter2(b, se, c, fix.c, l, alpha, x, y, weights, sigma2, 
  cov.unscaled, n, intercept, tol, unpenalized, offset, family, maxIter)
interpolate(y1, y2, x1, x2, x.new)

makeHyp(cnames, hypothesis, rhs = NULL) 
printHyp(L, b, nms)

cislasso(object, a, ci)
ci.fitted.islasso(object, newx, ci = NULL, type.ci = "wald", 
  conf.level=.95, only.ci = FALSE)
  
grad.islasso(object, ...)
hess.islasso(object, ...)
is.influence(model, do.coef = TRUE)
islasso.diag(glmfit)
islasso.diag.plots(glmfit, glmdiag = islasso.diag(glmfit), 
  subset = NULL, iden = FALSE, labels = NULL, ret = FALSE)
predislasso(object, newdata, type = c("response", "terms"), 
  terms = NULL, na.action = na.pass, ...)
}
\author{Gianluca Sottile (\email{gianluca.sottile@unipa.it})}
\details{These functions are not intended for users.}
\keyword{internal}
