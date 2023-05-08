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

summary.islasso <- function(object, pval = 1, which, use.t = FALSE, type.pval = "wald", ...){
  temp <- list(...)
  type.pval <- match.arg(type.pval)
  if(type.pval != "wald") stop("Only Wald-type confidence intervals are implemented yet!")
  
  n <- object$internal$n
  p <- object$internal$p
  if(missing(which)) which <- seq_len(p)
  coef <- object$coef[which]
  se <- object$se[which]
  aic <- object$aic
  
  res <- residuals(object, type="deviance")
  dispersion <- object$dispersion
  rdf <- object$df.residual
  df <- c(object$internal$p, rdf)
  h <- object$internal$hi[which]
  if(object$family$family != "gaussian") use.t <- FALSE
  
  if(type.pval == "wald"){
    chival <- (coef / se)
  }
  
  coefficients <- round(cbind(coef, se, h, chival), 6)
  type <- if(use.t) "t value" else "z value"
  pvalue <- if(type == "t value") 2 * pt(abs(chival), rdf, lower.tail = FALSE) else 2 * pnorm(-abs(chival))
  ptype <- if(type == "t value") "Pr(>|t|)" else "Pr(>|z|)"
  coefficients <- cbind(coefficients, pvalue)
  
  dimnames(coefficients) <- list(object$internal$nms[which], c("Estimate", "Std. Error", "Df", type, ptype))
  
  out <- list(coefficients = coefficients, dispersion = dispersion, df = df, 
              res = res, aic = aic, lambda = object$lambda, nulldev = object$null.deviance,
              dev = object$deviance, df.null = object$df.null, 
              df.res = object$df.null - (object$rank - 1*object$internal$intercept*object$internal$hi[1]),
              family = object$family, iter = object$iter, pval = pval, unpenalized = object$internal$unpenalized[which], 
              temp = temp, call = object$call)
  
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
  # if(rdf > 5L){
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    zz <- zapsmall(quantile(resid), digits + 1L)
    rq <- structure(zz, names = nam)
    print(rq, digits = digits, ...)
  # }else{
  #   if(rdf > 0L){
  #     print(resid, digits = digits, ...)
  #   }else{
  #     cat("ALL", df[1L], "residuals are 0: no residual degrees of freedom!")
  #   }
  # }
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
  y <- x$y
  eta <- x$linear.predictors
  mu <- x$fitted.values
  dev <- residuals(x, type = "deviance")
  sdev <- dev / sqrt(x$dispersion)
  
  pea <- residuals(x, type = "pearson")
  spea <- pea / sqrt(x$dispersion)
  
  w2 <- x$weights
  invH <- x$internal$invH
  
  w <- sqrt(w2)
  WX <- nX * w
  
  H <- WX %*% invH %*% t(WX)
  h <- diag(H)
  p <- x$rank
  rp <- spea / sqrt(1 - h)
  rd <- sdev / sqrt(1 - h)
  res <- sign(sdev) * sqrt(sdev^2 + h * rp^2)
  
  opar <- par()$mfrow
  on.exit(par(mfrow=opar))
  par(mfrow = c(2, 2))
  plot(eta, res, cex=L$cex, pch=L$pch, cex.lab=L$cex.lab, cex.axis=L$cex.axis, 
       xlab="Linear predictor", ylab="Residuals")
  
  qqNorm(na.omit(rd), cex.lab=L$cex.lab, cex.axis=L$cex.axis, 
         ylab="Ordered deviance residuals")
  
  plot(spea^2 ~ mu, cex=L$cex, pch=L$pch, cex.lab=L$cex.lab, cex.axis=L$cex.axis, 
       xlab="Fitted values", ylab="Squared stand Pearson residuals", main="Variance function")
  
  zmod <- eta + (y - mu) / w2
  plot(zmod ~ eta, cex=L$cex, pch=L$pch, cex.lab=L$cex.lab, cex.axis=L$cex.axis, 
       xlab="Linear predictor", ylab="Working vector", main="Link function")
  abline(lm(zmod ~ eta), col=L$col, lwd=L$lwd)
  
  invisible()
}

predict.islasso <- function (object, newdata = NULL, 
                             type = c("link", "response", "coefficients", "class", "terms"), 
                             se.fit = FALSE, ci = NULL, type.ci = "wald", level = 0.95,
                             terms = NULL, na.action = na.pass, ...) {
  type <- match.arg(type)
  type.ci <- match.arg(type.ci)
  if(type.ci != "wald") stop("Only Wald-type confidence intervals are implemented yet!")
  ci.fit <- NULL
  na.act <- object$na.action
  object$na.action <- NULL
  
  if(type == "class" & family(object)$family != "binomial") 
    stop(gettextf("Type 'class' is available only for the %s family", sQuote(binomial()$family)), domain=NA)
  
  if(type == "coefficients") {
    fit <- coef(object)
    if(se.fit) ci.fit <- ci.fitted.islasso(object, X, ci, type.ci = type.ci, conf.level = level, only.ci = TRUE)
  } 
  else {
    if (missing(newdata)) {
      fit <- switch(type, 
                    "response" = , "class" = , 
                    "link" = object$linear.predictors, 
                    "terms" = predislasso(object, type = "terms", terms = terms))
      X <- model.matrix(object)
    }
    else {
      fit <- predislasso(object, newdata, 
                          type = if (type %in% c("link", "class")) 
                            "response"
                          else type, terms = terms, na.action = na.action)
      X <- attr(fit, "X")
    }
    attr(fit, "X") <- NULL
    if (missing(newdata) && !is.null(na.act)) fit <- napredict(na.act, fit)
    if(se.fit & !(type %in% c("class", "terms"))) {
      ci.fit <- ci.fitted.islasso(object, X, ci, type.ci = type.ci, conf.level = level)
      if (missing(newdata) && !is.null(na.act)) {
        ci.fit[, 1L] <- napredict(na.act, ci.fit[, 1L])
        ci.fit[, 2L] <- napredict(na.act, ci.fit[, 2L])
      }
    }
  }
  out <- if(type != "terms") cbind("Fit" = fit, ci.fit) else fit
  
  out <- switch(type, 
         "response" = family(object)$linkinv(out), 
         "class" = {
           mu <- family$linkinv(out)
           1 * I(mu >= .5)
         },
         "terms" = out,
         "coefficients" =, "link" = drop(out))
  
  drop(out)
}

# model.matrix.islasso <- function (object, ...) {
#   model.matrix(object$terms, object$model, object$contrasts, ...)
# }

model.matrix.islasso <- function (object, ...) {
  data <- model.frame(object, xlev = object$xlevels, ...)
  if (exists(".GenericCallEnv", inherits = FALSE)) 
    NextMethod("model.matrix", data = data, contrasts.arg = object$contrasts)
  else {
    dots <- list(...)
    dots$data <- dots$contrasts.arg <- NULL
    do.call("model.matrix.default", c(list(object = object, 
                                           data = data, contrasts.arg = object$contrasts), 
                                      dots))
  }
}

# coef.islasso <- function(object, ...){
#   temp <- list(...)
#   if(is.null(temp$unbias)) temp$unbias <- FALSE
#   if(!temp$unbias) return(object$coef) else return(object$beta.unbias)
# }

vcov.islasso <- function(object, ...) object$internal$vcov

# fitted.islasso <- function(object, ...) object$fitted.values

residuals.islasso <- function (object, type = c("deviance", "pearson", 
                                                "working", "response", "partial"), ...) {
  type <- match.arg(type)
  y <- object$y
  r <- object$residuals
  mu <- object$fitted.values
  wts <- object$prior.weights
  switch(type, deviance = , pearson = , 
         response = if (is.null(y)) {
           mu.eta <- object$family$mu.eta
           eta <- object$linear.predictors
           y <- mu + r * mu.eta(eta)
  })
  res <- switch(type, 
                deviance = if (object$df.residual > 0) {
                  d.res <- sqrt(pmax((object$family$dev.resids)(y, mu, wts), 0))
                  ifelse(y > mu, d.res, -d.res)
                  } else rep.int(0, length(mu)), 
                pearson = (y - mu) * sqrt(wts) / sqrt(object$family$variance(mu)), 
                working = r, response = y - mu, partial = r)
  if (!is.null(object$na.action)) 
    res <- naresid(object$na.action, res)
  if (type == "partial") 
    res <- res + predict(object, type = "terms")
  res
}

deviance.islasso <- function(object, ...) object$deviance

# logLik.islasso <- function(object, ...){
#   p <- object$rank
#   val <- p - object$aic/2
#   
#   out <- list(loglik=val, df=p, object=object, phi=object$phi)
#   class(out) <- "logLik.islasso"
#   out
# }

logLik.islasso <- function (object, ...) {
  if (!missing(...)) warning("extra arguments discarded")
  fam <- family(object)$family
  p <- object$rank
  if (fam %in% c("gaussian", "Gamma", "inverse.gaussian")) 
    p <- p + 1
  val <- p - object$aic/2
  attr(val, "nobs") <- sum(!is.na(object$residuals))
  attr(val, "df") <- p
  class(val) <- "logLik"
  val
}

# print.logLik.islasso <- function(x, digits=max(3L, getOption("digits") - 3L), ...){
#   cat("\n'log Lik.' ", paste(format(c(x$loglik), digits=digits), collapse = ", "),
#       " (df=", format(x$df), ")\n\n", sep = "")
#   invisible(x)
# }

# AIC.islasso <- function(object, k = 2, ...){
#   df <- object$rank
#   ll <- logLik(object)$loglik
#   aic <- -2*ll + k*df
#   return(aic)
# }

extractAIC.islasso <- function (fit, scale = 0, k = 2, ...) {
  n <- length(fit$residuals)
  edf <- n - fit$df.residual
  aic <- fit$aic
  c(edf, aic + (k - 2) * edf)
}

family.islasso <- function (object, ...) object$family

formula.islasso <- function (x, ...) {
  form <- x$formula
  if (!is.null(form)) {
    form <- formula(x$terms)
    environment(form) <- environment(x$formula)
    form
  }
  else formula(x$terms)
}

model.frame.islasso <- function (formula, ...) {
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]
  if (length(nargs) || is.null(formula$model)) {
    fcall <- formula$call
    fcall$method <- "model.frame"
    fcall[[1L]] <- quote(stats::glm)
    fcall[names(nargs)] <- nargs
    env <- environment(formula$terms) %||% parent.frame()
    eval(fcall, env)
  }
  else formula$model
}

nobs.islasso <- function (object, ...) if (!is.null(w <- object$prior.weights)) sum(w != 0) else length(object$residuals)

weights.islasso <- function (object, type = c("prior", "working"), ...) {
  type <- match.arg(type)
  res <- if (type == "prior") 
    object$prior.weights
  else object$weights
  if (is.null(object$na.action)) 
    res
  else naresid(object$na.action, res)
}

variable.names.islasso <- function (object, ...) object$internal$nms

influence.islasso <- function (model, do.coef = TRUE, ...) {
  res <- is.influence(model, do.coef = do.coef, ...)
  pRes <- na.omit(residuals(model, type = "pearson"))[model$prior.weights != 0]
  pRes <- naresid(model$na.action, pRes)
  names(res)[names(res) == "wt.res"] <- "dev.res"
  c(res, list(pear.res = pRes))
}

cooks.distance.islasso <- function (model, infl = influence(model, do.coef = FALSE), res = infl$pear.res, 
                                    dispersion = summary(model)$dispersion, hat = infl$hat, ...) {
  p <- model$rank
  res <- (res/(1 - hat))^2 * hat/(dispersion * p)
  res[is.infinite(res)] <- NaN
  res
}

rstandard.islasso <- function (model, infl = influence(model, do.coef = FALSE), 
                               type = c("deviance", "pearson"), ...) {
  type <- match.arg(type)
  res <- switch(type, pearson = infl$pear.res, infl$dev.res)
  res <- res/sqrt(summary(model)$dispersion * (1 - infl$hat))
  res[is.infinite(res)] <- NaN
  res
}

rstudent.islasso <- function (model, infl = influence(model, do.coef = FALSE), ...) {
  r <- infl$dev.res
  r <- sign(r) * sqrt(r^2 + (infl$hat * infl$pear.res^2)/(1 - infl$hat))
  r[is.infinite(r)] <- NaN
  if (any(family(model)$family == c("binomial", "poisson"))) 
    r
  else r/infl$sigma
}
