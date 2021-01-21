print.islasso <- function(x, digits=max(3L, getOption("digits") - 3L), ...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if(length(coef(x))){
    cat("Coefficients:\n")
    print.default(format(round(x$coef, 6), digits = digits), print.gap = 2, quote = FALSE)
  }else{cat("No coefficients\n\n")}
  cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ", format(signif(x$df.null-(x$rank-1*x$internal$intercept*x$internal$hi[1]), digits)), "Residual\n")
  cat("Null Deviance:", format(signif(x$null.deviance, digits)), "\nResidual Deviance:", format(signif(x$deviance, digits)), "\nAIC:", format(signif(x$aic, digits)), "\nLambda:", format(signif(x$lambda, digits)))
  
  cat("\n")
  invisible(x)
}

summary.islasso <- function(object, pval=1, use.t=FALSE, ...){
  temp <- list(...)
  
  #if(is.null(temp$unbias)) temp$unbias <- FALSE
  #if(temp$unbias & is.null(object$beta.unbias)) stop("No unbias estimates in the model")
  
  coef <- object$coef #if(!temp$unbias) object$coef else object$beta.unbias
  se <- object$se #if(!temp$unbias) object$se else object$se.unbias
  aic <- object$aic
  n <- object$internal$n
  p <- object$internal$p
  
  res <- residuals(object, type="deviance")
  dispersion <- object$dispersion
  rdf <- object$df.residual
  df <- c(object$internal$p, rdf)
  h <- object$internal$hi
  if(object$family$family != "gaussian") use.t <- FALSE
  chival <- (coef / se)
  coefficients <- round(cbind(coef, se, h, chival), 6)
  type <- if(use.t) "t value" else "z value"
  pvalue <- if(type == "t value") 2 * pt(abs(chival), rdf, lower.tail = FALSE) else 2 * pnorm(-abs(chival))
  ptype <- if(type == "t value") "Pr(>|t|)" else "Pr(>|z|)"
  coefficients <- cbind(coefficients, pvalue)
  
  dimnames(coefficients) <- list(object$internal$nms, c("Estimate", "Std. Error", "Df", type, ptype))
  
  out <- list(coefficients=coefficients, dispersion=dispersion, df=df, res=res, aic=aic, lambda=object$lambda, nulldev=object$null.deviance,
              dev=object$deviance, df.null=object$df.null, df.res=object$df.null-(object$rank-1*object$internal$intercept*object$internal$hi[1]),
              family=object$family, iter=object$iter, pval=pval, unpenalized=object$internal$unpenalized, 
              temp=temp, call=object$call)
  
  class(out) <- "summary.islasso"
  return(out)
}

print.summary.islasso <- function(x, digits=max(3L, getOption("digits") - 3L), ...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if(!is.null(x$temp$digits)) digits <- x$temp$digits
  resid <- x$res
  df <- x$df
  rdf <- df[2L]
  cat("Residuals:\n", sep="")
  if(rdf > 5L){
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    zz <- zapsmall(quantile(resid), digits + 1L)
    rq <- structure(zz, names = nam)
    print(rq, digits = digits, ...)
  }else{
    if(rdf > 0L){
      print(resid, digits = digits, ...)
    }else{
      cat("ALL", df[1L], "residuals are 0: no residual degrees of freedom!")
    }
  }
  cat("\n")
  coefs <- x$coefficients
  if(sum(coefs[, 5] <= x$pval) == 0) warning("No coefficients lower than the selected p-value. The lowest p-value is printed.")
  temp <- (coefs[, 5] <= x$pval) | x$unpenalized
  if(sum(temp) == 0) temp <- coefs[, 5] == min(coefs[, 5])
  coefs <- coefs[temp, , drop = FALSE]
  # if(sum(temp) == 1){
  #   coefs <- t(coefs)
  #   rownames(coefs) <- rownames(x$coefficients)[temp]
  # }
  printCoefmat(coefs, digits=digits, signif.stars=TRUE, has.Pvalue=TRUE, na.print="NA", cs.ind = 1L:2L, tst.ind = 3L:4L,, ...)
  cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ", format(x$dispersion), ")\n\n",
      apply(cbind(paste(format(c("Null", "Residual"), justify = "right"), "deviance:"),
                  format(c(x$nulldev, x$dev), digits = max(5L, digits + 1L)), " on",
                  format(c(x$df.null, x$df.res), digits = digits), " degrees of freedom\n"), 1L, paste, collapse = " "), sep = "")
  cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)), "\nLambda: ", format(x$lambda, digits = max(4L, digits + 1L)), "\n\n", 
      "Number of Newton-Raphson iterations: ", x$iter, "\n", sep = "")
  
  cat("\n")
  invisible(x)
}

plot.islasso <- function(x, ...){
  L <- list(...)
  if(is.null(L$lty)) L$lty <- 3
  if(is.null(L$cex)) L$cex <- 1
  if(is.null(L$pch)) L$pch <- 19
  if(is.null(L$col)) L$col <- 2
  if(is.null(L$lwd)) L$lwd <- 2
  if(is.null(L$cex.lab)) L$cex.lab <- 1
  if(is.null(L$cex.axis)) L$cex.axis <- 1
  
  family <- x$family
  nX <- model.matrix(x)
  y <- x$internal$y.internal
  eta <- x$linear.predictors
  mu <- x$fitted.values
  dev <- residuals(x, type="deviance")
  sdev <- dev/sqrt(x$dispersion)
  
  pea <- residuals(x, type="pearson")
  spea <- pea/sqrt(x$dispersion)
  
  w2 <- x$internal$wt
  invH <- x$internal$invH
  
  w <- sqrt(w2)
  WX <- nX*w
  
  H <- WX %*% invH %*% t(WX)
  h <- diag(H)
  p <- x$rank
  rp <- spea/sqrt(1 - h)
  rd <- sdev/sqrt(1 - h)
  res <- sign(sdev) * sqrt(sdev^2 + h * rp^2)
  
  opar <- par()$mfrow
  on.exit(par(mfrow=opar))
  par(mfrow=c(2,2))
  plot(eta, res, cex=L$cex, pch=L$pch, cex.lab=L$cex.lab, cex.axis=L$cex.axis, xlab="Linear predictor", ylab="Residuals")
  
  qqNorm(na.omit(rd), cex.lab=L$cex.lab, cex.axis=L$cex.axis, ylab="Ordered deviance residuals")
  
  plot(spea^2 ~ mu, cex=L$cex, pch=L$pch, cex.lab=L$cex.lab, cex.axis=L$cex.axis, xlab="Fitted values", ylab="Squared stand Pearson residuals", main="Variance function")
  
  zmod <- eta + (y-mu)/(w2)
  plot(zmod ~ eta, cex=L$cex, pch=L$pch, cex.lab=L$cex.lab, cex.axis=L$cex.axis, xlab="Linear predictor", ylab="Working vector", main="Link function")
  abline(lm(zmod ~ eta), col=L$col, lwd=L$lwd)
}

predict.islasso <- function(object, newdata, type=c("link", "response", "coefficients", "class"), 
                            se.fit = FALSE, ci = NULL, level = .95, ...){
  # beta <- as.matrix(coef(object))
  # intercept <- object$internal$intercept
  type <- match.arg(type)
  family <- object$family
  
  tt <- terms(object)
  if (missing(newdata) || is.null(newdata)) {
    X <- model.matrix(object)
    offset <- object$offset
  }
  else {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.pass, xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    offset <- rep(0, nrow(X))
    if (!is.null(off.num <- attr(tt, "offset"))) 
      for (i in off.num) offset <- offset + eval(attr(tt, "variables")[[i + 1]], newdata)
    if (!is.null(object$call$offset)) 
      offset <- offset + eval(object$call$offset, newdata)
  }
  beta <- object$coefficients
  predictor <- drop(X %*% beta)
  if (!is.null(offset)) predictor <- predictor + offset
  
  # if(missing(newdata)){
  #   X <- model.matrix(object)
  # }else{
  #   X <- as.matrix(newdata)
  #   if(!all(colnames(X) %in% rownames(beta))) stop("'newdata' must contain all X-variables")
  #   if(intercept) X <- cbind(1, X)
  # }
  
  # if(nrow(beta) != ncol(X)) stop("the length of the vector 'beta' is not equal to the number of columns of the matrix 'X'")
  if(type == "class" & family$family != "binomial") stop(gettextf("Type 'class' is available only for the %s family", sQuote(family$family)), domain=NA)
  
  out <- switch(type,
                "link"={drop(predictor)},
                "response"={family$linkinv(predictor)},
                "coefficients"={coef(object)},
                "class"={
                  mu <- family$linkinv(predictor)
                  1*I(mu >= .5)})
  
  if(se.fit){
    ci.fit <- NULL
    if(type %in% c("link", "response")) ci.fit <- ci.fitted.islasso(object, X, ci, conf.level = level)
    if(type == "coefficients") ci.fit <- ci.fitted.islasso(object, X, ci, conf.level = level, only.ci = TRUE)
    
    out.ci <- switch(type,
                     "link"={ci.fit},
                     "response"={family$linkinv(ci.fit)},
                     "coefficients"={ci.fit})
    out <- cbind("Fit" = out, out.ci)
  }
  
  return(out)
}

model.matrix.islasso <- function (object, ...) {
  model.matrix(object$terms, object$model, object$contrasts, ...)
}

coef.islasso <- function(object, ...){
  temp <- list(...)
  if(is.null(temp$unbias)) temp$unbias <- FALSE
  if(!temp$unbias) return(object$coef) else return(object$beta.unbias)
}

vcov.islasso <- function(object, ...){
  object$internal$vcov
}

fitted.islasso <- function(object, ...){
  return(object$fitted.values)
}

residuals.islasso <- function(object, type=c("working", "response", "deviance", "pearson"), ...){
  type <- match.arg(type)
  y <- object$internal$y.internal
  mu <- object$fitted.values
  weights <- object$internal$weights
  family <- object$family
  
  switch(type,
         "working"=return(object$residuals),
         "response"=return(y-mu),
         "deviance"=return(sign(y-mu)*sqrt(family$dev.resids(y, mu, weights))),
         "pearson"=return((y-mu)*sqrt(weights)/sqrt(family$variance(mu))))
}

deviance.islasso <- function(object, ...){
  return(object$deviance)
}

logLik.islasso <- function(object, ...){
  p <- object$rank
  val <- p - object$aic/2
  
  out <- list(loglik=val, df=p, object=object, phi=object$phi)
  class(out) <- "logLik.islasso"
  out
}

print.logLik.islasso <- function(x, digits=max(3L, getOption("digits") - 3L), ...){
  cat("\n'log Lik.' ", paste(format(c(x$loglik), digits=digits), collapse = ", "),
      " (df=", format(x$df), ")\n\n", sep = "")
  invisible(x)
}

AIC.islasso <- function(object, k=2, ...){
  df <- object$rank
  ll <- logLik(object)$loglik
  aic <- -2*ll + k*df
  return(aic)
}
