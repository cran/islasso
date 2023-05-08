### S3 methods for islasso.path class

interpolate <- function(y1, y2, x1, x2, x.new) {
  m <- (y2 - y1) / (log(x2) - log(x1))
  return(y1 + m * (log(x.new) - log(x1)))
}

print.islasso.path <-  function(x, digits = max(3L, getOption("digits") - 3L), ...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if(length(x$Coef)){
    cat("Coefficients:\n")
    Info <- data.frame(x$Info[, 1:5, drop = FALSE])
    Info$lambda <- formatC(Info$lambda, format = "f", digits = 4, width = 4 + nchar(round(max(Info$lambda))))
    Info$df <- formatC(Info$df, format = "f", digits = 4, width = 4 + nchar(round(max(Info$df))))
    Info$phi <- formatC(Info$phi, format = "f", digits = 4, width = 4 + nchar(round(max(Info$phi))))
    Info$deviance <- formatC(Info$deviance, format = "f", digits = 4, width = 4 + nchar(round(max(Info$deviance))))
    Info$logLik <- formatC(Info$logLik, format = "f", digits = 4, width = 4 + nchar(round(max(Info$logLik))))
    # Info$iterB <- formatC(Info$iterB, format = "d", digits = 0, width = 0 + nchar(round(max(Info$iterB))))
    # Info$iterSE <- formatC(Info$iterSE, format = "f", digits = 0, width = 0 + nchar(round(max(Info$iterSE))))
    # Info$errB <- formatC(Info$errB, format = "f", digits = 7, width = 7 + nchar(round(max(Info$errB))))
    # Info$errSE <- formatC(Info$errSE, format = "f", digits = 7, width = 7 + nchar(round(max(Info$errSE))))
    print(Info)
  }else{cat("No coefficients\n\n")}
  
  cat("\n")
  invisible(x)
}

summary.islasso.path <- function(object, pval = 1, use.t = FALSE, lambda, ...){
  temp <- list(...)
  
  lambda.seq <- object$Info[,"lambda"]
  nlambda <- length(lambda.seq)
  if(nlambda == 1){
    lambda <- lambda.seq
    coef <- object$Coef[1,]
    se <- object$SE[1,]
    aic <- object$GoF[1, "AIC"]
    dispersion <- object$Info[1, "phi"]
    df0 <- object$Info[1, "df"]
    logLik <- object$Info[1, "logLik"]
    iter <- object$Info[1, "iter"]
  }
  else{
    if(missing(lambda)) return(print(object))
    if(lambda < min(lambda.seq) | lambda > max(lambda.seq)) stop("value of lambda out of bound")
    id1 <- rev(which(lambda >= lambda.seq))[1]
    id2 <- which(lambda < lambda.seq)[1]
    if(any(is.na(c(id1, id2)))) {
      id1 <- c(id1, id2)[!is.na(c(id1, id2))]
      coef <- object$Coef[id1,]
      se <- object$SE[id1,]
      aic <- object$GoF[id1, "AIC"]
      dispersion <- object$Info[id1, "phi"]
      df0 <- object$Info[id1, "df"]
      logLik <- object$Info[id1, "logLik"]
      iter <- object$Info[id1, "iter"]
    }
    else {
      coef <- interpolate(object$Coef[id1,], object$Coef[id2,], lambda.seq[id1], lambda.seq[id2], lambda)
      se <- interpolate(object$SE[id1,], object$SE[id2,], lambda.seq[id1], lambda.seq[id2], lambda)
      aic <- interpolate(object$GoF[id1, "AIC"], object$GoF[id2, "AIC"], lambda.seq[id1], lambda.seq[id2], lambda)
      dispersion <- interpolate(object$Info[id1, "phi"], object$Info[id2, "phi"], lambda.seq[id1], lambda.seq[id2], lambda)
      df0 <- interpolate(object$Info[id1, "df"], object$Info[id2, "df"], lambda.seq[id1], lambda.seq[id2], lambda)
      logLik <- interpolate(object$Info[id1, "logLik"], object$Info[id2, "logLik"], lambda.seq[id1], lambda.seq[id2], lambda)
      iter <- object$Info[id1, "iter"]
    }
  }
  
  
  # coef <- object$Coef[lambda,]
  # se <- object$SE[lambda,]
  # aic <- object$Info[lambda, "AIC"]
  n <- object$Input$n
  p <- object$Input$p
  
  x <- model.matrix(object)
  y <- object$y
  offset <- object$offset
  family <- object$family
  
  eta <- drop(x %*% coef) + offset
  mu <- family$linkinv(eta)
  res <- drop(y - mu)
  # dispersion <- object$Info[lambda, "phi"]
  rdf <- n - df0
  df <- c(p, rdf)
  
  chival <- (coef / se)
  coefficients <- round(cbind(coef, se, chival), 6)
  type <- if(use.t) "t value" else "z value"
  pvalue <- if(type == "t value") 2 * pt(abs(chival), rdf, lower.tail = FALSE) else 2 * pnorm(-abs(chival))
  ptype <- if(type == "t value") "Pr(>|t|)" else "Pr(>|z|)"
  coefficients <- cbind(coefficients, pvalue)
  
  dimnames(coefficients) <- list(colnames(object$Coef), c("Estimate", "Std. Error", type, ptype))
  
  # null model
  intercept <- object$Input$intercept
  weights <- object$prior.weights
  wtdmu <- if(intercept) sum(weights * y)/sum(weights) else family$linkinv(offset)
  nulldev <- sum(family$dev.resids(y, wtdmu, weights))
  n.ok <- n - sum(weights == 0)
  nulldf <- n.ok - as.integer(intercept)
  
  out <- list(coefficients = coefficients, dispersion = dispersion, df = df, res = res, aic = aic, lambda = lambda, 
              nulldev = nulldev, dev = logLik, df.null = nulldf, df.res = nulldf - (df0 - 1*object$Input$intercept),
              family = object$Input$family$family, iter = iter, pval = pval, temp = temp, call = object$call)
  
  class(out) <- "summary.islasso.path"
  return(out)
}

print.summary.islasso.path <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
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
  if(sum(coefs[, 4] <= x$pval) == 0) warning("No coefficients lower than the selected p-value. The lowest p-value is printed.")
  temp <- (coefs[, 4] <= x$pval)
  if(sum(temp) == 0) temp <- coefs[, 4] == min(coefs[, 4])
  coefs <- coefs[temp, , drop = FALSE]
  
  printCoefmat(coefs, digits=digits, signif.stars=TRUE, has.Pvalue=TRUE, na.print="NA", cs.ind = 1L:2L, tst.ind = 3L,, ...)
  cat("\n(Dispersion parameter for ", x$family, " family taken to be ", format(x$dispersion), ")\n\n",
      apply(cbind(paste(format(c("Null", "Residual"), justify = "right"), "deviance:"),
                  format(c(x$nulldev, x$dev), digits = max(5L, digits + 1L)), " on",
                  format(c(x$df.null, x$df.res), digits = digits), " degrees of freedom\n"), 1L, paste, collapse = " "), sep = "")
  cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)), "\nLambda: ", format(x$lambda, digits = max(4L, digits + 1L)), "\n\n", 
      "Number of Newton-Raphson iterations: ", x$iter, "\n", sep = "")
  
  cat("\n")
  invisible(x)
}

plot.islasso.path <- function (x, yvar = c("coefficients", "se", "gradient", "weight", "gof"), 
                               gof = c("none", "AIC", "BIC", "AICc", "eBIC", "GCV", "GIC"), 
                               label = FALSE, legend = FALSE, ...) {
  object <- x
  if(!inherits(object, "islasso.path")) stop("x is not an object of islasso.path class")
  label <- abs(label)
  yvar <- match.arg(yvar)
  gof <- match.arg(gof)
  if(yvar == "gof" & gof == "none") stop("You have to select one criterion between AIC, BIC, AICc, eBIC, GCV, GIC")
  id.best <- apply(object$GoF, 2, which.min)
  lambda <- object$Info[, "lambda"]
  loglambda <- log(lambda)
  nlambda <- length(lambda)
  if(nlambda == 1) stop("no path in islasso.path fit!")
  
  intercept <- object$Input$intercept
  coef1 <- object$Coef
  if(intercept) coef1 <- coef1[, -1]
  unactive <- abs(coef1[id.best[gof], ]) <= 1E-6
  active <- !unactive
  if(gof == "none") active <- rep(TRUE, length(active))
  
  dots <- list(...)
  if(is.null(dots$main)) dots$main <- ""
  if(is.null(dots$lty)) dots$lty <- if(gof != "none") ifelse(unactive, 2L, 1L) else 1L
  if(is.null(dots$col)) dots$col <- if(gof != "none") ifelse(unactive, "gray80", "gray20") else "gray50"
  if(is.null(dots$lwd)) dots$lwd <- 1L
  if(is.null(dots$gof_lty)) dots$gof_lty <- 3L
  if(is.null(dots$gof_col)) dots$gof_col <- "red"
  if(is.null(dots$gof_lwd)) dots$gof_lwd <- 1L
  if(is.null(dots$cex.axis)) dots$cex.axis <- 1L
  if(is.null(dots$cex.lab)) dots$cex.lab <- 1L
  if(is.null(dots$cex.main)) dots$cex.main <- 1L
  
  if (yvar == "coefficients") {
    if(is.null(dots$xlab)) dots$xlab <- "log Lambda"
    if(is.null(dots$ylab)) dots$ylab <- "Coefficients"
    if(is.null(dots$xlim)) dots$xlim <- range(loglambda) + c(-label, legend)
    if(is.null(dots$ylim)) dots$ylim <- range(coef1[1, ])
    
    matplot(loglambda, coef1, type = "l", xlab = dots$xlab, ylab = dots$ylab, main = dots$main, 
            xlim = dots$xlim, ylim = dots$ylim, lty = dots$lty, col = dots$col, lwd = dots$lwd,
            cex.axis = dots$cex.axis, cex.lab = dots$cex.lab, cex.main = dots$cex.main)
    abline(h = 0, lty = 1, col = "gray50")
    if (label) text(x = min(loglambda) - label, y = coef1[1, active], 
                    labels = colnames(coef1)[active], cex = dots$cex.lab, xpd = TRUE, pos = 4)
    if(gof != "none") abline(v = loglambda[id.best[gof]], col = dots$gof_col, lty = dots$gof_lty, lwd = dots$gof_lwd)
    if(legend & gof != "none") legend('topright', legend = paste0("min.", gof), lty = dots$gof_lty, col = dots$gof_col, 
                                      lwd = dots$gof_lwd, bty = "n", cex = dots$cex.lab)
  }
  if (yvar == "se") {
    se1 <- object$SE
    if(intercept) se1 <- se1[, -1]
    
    if(is.null(dots$xlab)) dots$xlab <- "log Lambda"
    if(is.null(dots$ylab)) dots$ylab <- "Std.Errs"
    if(is.null(dots$xlim)) dots$xlim <- range(loglambda) + c(-label, legend)
    if(is.null(dots$ylim)) dots$ylim <- range(se1)
    
    matplot(loglambda, se1, type = "l", xlab = dots$xlab, ylab = dots$ylab, main = dots$main, 
            xlim = dots$xlim, ylim = dots$ylim, lty = dots$lty, col = dots$col, lwd = dots$lwd,
            cex.axis = dots$cex.axis, cex.lab = dots$cex.lab, cex.main = dots$cex.main)
    abline(h = 0, lty = 1, col = "gray50")
    if (label) text(x = min(loglambda) - label, y = se1[1, active], 
                    labels = colnames(coef1)[active], cex = dots$cex.lab, xpd = TRUE, pos = 4)
    if(gof != "none") abline(v = loglambda[id.best[gof]], col = dots$gof_col, lty = dots$gof_lty, lwd = dots$gof_lwd)
    if(legend & gof != "none") legend('topright', legend = paste0("min.", gof), lty = dots$gof_lty, col = dots$gof_col, 
                                      lwd = dots$gof_lwd, bty = "n", cex = dots$cex.lab)
  }
  if (yvar == "weight") {
    weight1 <- object$Weight
    if(intercept) weight1 <- weight1[, -1]
    
    if(is.null(dots$xlab)) dots$xlab <- "log Lambda"
    if(is.null(dots$ylab)) dots$ylab <- "Mixture weight"
    if(is.null(dots$xlim)) dots$xlim <- range(loglambda) + c(-label, legend)
    if(is.null(dots$ylim)) dots$ylim <- range(weight1)
    
    matplot(loglambda, weight1, type = "l", xlab = dots$xlab, ylab = dots$ylab, main = dots$main, 
            xlim = dots$xlim, ylim = dots$ylim, lty = dots$lty, col = dots$col, lwd = dots$lwd,
            cex.axis = dots$cex.axis, cex.lab = dots$cex.lab, cex.main = dots$cex.main)
    abline(h = 0, lty = 1, col = "gray50")
    if (label) text(x = min(loglambda) - label, y = weight1[1, active], 
                    labels = colnames(coef1)[active], cex = dots$cex.lab, xpd = TRUE, pos = 4)
    if(gof != "none") abline(v = loglambda[id.best[gof]], col = dots$gof_col, lty = dots$gof_lty, lwd = dots$gof_lwd)
    if(legend & gof != "none") legend('topright', legend = paste0("min.", gof), lty = dots$gof_lty, col = dots$gof_col, 
                                      lwd = dots$gof_lwd, bty = "n", cex = dots$cex.lab)
  }
  if (yvar == "gradient") {
    fn0 <- function(x, y, b, s = 1, c = .5, lambda, alpha, unpenalized, family, offset, weights) {
      eta <- x %*% b + offset
      mu <- family$linkinv(eta)
      v <- family$variance(mu)
      m.e <- family$mu.eta(eta)
      r <- y - mu
      rv <- r / v * m.e * weights
      grad <- drop(t(rv) %*% x)
      
      bsc <- (b / s)
      # grad <- crossprod(x, y - drop(x %*% b))
      r <- alpha * (c*(2 * pnorm(bsc, 0, 1) - 1) + (1 - c)*(2 * pnorm(bsc, 0, 1E-5) - 1)) + (1 - alpha) * b
      if(any(unpenalized)) r[unpenalized] <- 0
      return(- grad + lambda * r)
    }
    grad <- t(sapply(seq_len(nlambda), function(i) {
      fn0(x = model.matrix(object), y = object$y, b = object$Coef[i, ], 
          s = object$SE[i, ], c = object$Weight[i, ], 
          lambda = lambda[i], alpha = object$alpha, 
          unpenalized = object$Input$unpenalized, family = object$family, 
          offset = object$offset, weights = object$prior.weights)
    }))
    if(intercept) grad <- grad[, -1]
    
    if(is.null(dots$xlab)) dots$xlab <- "Lambda"
    if(is.null(dots$ylab)) dots$ylab <- "Gradient"
    if(is.null(dots$xlim)) dots$xlim <- range(lambda) + c(-legend, label)
    if(is.null(dots$ylim)) dots$ylim <- range(grad)
    
    matplot(lambda, grad, type = "l", xlab = dots$xlab, ylab = dots$ylab, main = dots$main, 
            xlim = dots$xlim, ylim = dots$ylim, lty = dots$lty, col = dots$col, lwd = dots$lwd,
            cex.axis = dots$cex.axis, cex.lab = dots$cex.lab, cex.main = dots$cex.main)
    abline(h = 0, lty = 1, col = "gray50")
    if (label) text(x = max(lambda) + label, y = grad[nlambda, active], 
                    labels = colnames(coef1)[active], cex = dots$cex.lab, xpd = TRUE, pos = 4)
    if(gof != "none") abline(v = lambda[id.best[gof]], col = dots$gof_col, lty = dots$gof_lty, lwd = dots$gof_lwd)
    if(legend & gof != "none") legend('topleft', legend = paste0("min.", gof), lty = dots$gof_lty, col = dots$gof_col, 
                                      lwd = dots$gof_lwd, bty = "n", cex = dots$cex.lab)
  }
  if (yvar == "gof") {
    if(is.null(dots$xlab)) dots$xlab <- "log Lambda"
    if(is.null(dots$ylab)) dots$ylab <- gof
    
    plot(loglambda, object$GoF[, gof, drop = FALSE], type = "l", ylab = dots$ylab, 
         xlab = dots$xlab, col = "gray50", lwd = dots$lwd, lty = 1L,
         cex.axis = dots$cex.axis, cex.lab = dots$cex.lab, cex.main = dots$cex.main)
    abline(v = loglambda[id.best[gof]], col = dots$gof_col, lty = dots$gof_lty, lwd = dots$gof_lwd)
  }
  
  invisible()
}

# plot.islasso.path <- function (x, yvar = c("coefficients", "se", "gradient", "weight", "AIC", "BIC", "AICc", "eBIC", "GCV", "GIC"), 
#                                label = FALSE, cex.lab = 1, ...) {
#   object <- x
#   label <- abs(label)
#   yvar <- match.arg(yvar)
#   id.best <- apply(object$GoF, 2, which.min)
#   lambda <- object$Info[, "lambda"]
#   nlambda <- length(lambda)
#   if(nlambda == 1) stop("no path in islasso.path fit!")
#   if (yvar == "AIC") {
#     # AIC
#     plot(log(lambda), object$GoF[, "AIC", drop = FALSE], type = 'l', ylab = 'AIC', xlab = 'log Lambda', col = 'gray', ...)
#     abline(v = log(lambda[id.best["AIC"]]), col = 2, lty = 2, lwd = 2)
#   }
#   if (yvar == "BIC") {
#     # BIC
#     plot(log(lambda), object$GoF[, "BIC", drop = FALSE], type = 'l', ylab = 'BIC', xlab = 'log Lambda', col = 'gray', ...)
#     abline(v = log(lambda[id.best["BIC"]]), col = 2, lty = 2, lwd = 2)
#   }
#   if (yvar == "AICc") {
#     # AICc
#     plot(log(lambda), object$GoF[, "AICc", drop = FALSE], type = 'l', ylab = 'AICc', xlab = 'log Lambda', col = 'gray', ...)
#     abline(v = log(lambda[id.best["AICc"]]), col = 2, lty = 2, lwd = 2)
#   }
#   if (yvar == "eBIC") {
#     # BIC
#     plot(log(lambda), object$GoF[, "eBIC", drop = FALSE], type = 'l', ylab = 'eBIC', xlab = 'log Lambda', col = 'gray', ...)
#     abline(v = log(lambda[id.best["eBIC"]]), col = 2, lty = 2, lwd = 2)
#   }
#   if (yvar == "GCV") {
#     # GCV
#     plot(log(lambda), object$GoF[, "GCV", drop = FALSE], type = 'l', ylab = 'GCV', xlab = 'log Lambda', col = 'gray', ...)
#     abline(v = log(lambda[id.best["GCV"]]), col = 2, lty = 2, lwd = 2)
#   }
#   if (yvar == "GIC") {
#     # GCV
#     plot(log(lambda), object$GoF[, "GIC", drop = FALSE], type = 'l', ylab = 'GIC', xlab = 'log Lambda', col = 'gray', ...)
#     abline(v = log(lambda[id.best["GIC"]]), col = 2, lty = 2, lwd = 2)
#   }
#   intercept <- object$Input$intercept
#   coef1 <- object$Coef
#   if(intercept) coef1 <- coef1[, -1]
#   nlambda <- length(lambda)
#   if (yvar == "coefficients") {
#     matplot(log(lambda), coef1, type = "l", xlab = "log Lambda", ylab = "Coefficients", xlim = range(log(lambda)) + c(-label, 0), ...)
#     abline(h = 0, lty = 1, col = 'gray')
#     if (label) text(x = min(log(lambda)) - label, y = coef1[1, ], labels = colnames(coef1), cex = cex.lab, xpd = TRUE, pos = 4)
#     abline(v = log(lambda[id.best["AIC"]]), col = 2, lty = 3, lwd = 2)
#     abline(v = log(lambda[id.best["AICc"]]), col = 3, lty = 3, lwd = 2)
#     abline(v = log(lambda[id.best["BIC"]]), col = 4, lty = 3, lwd = 2)
#     abline(v = log(lambda[id.best["eBIC"]]), col = 5, lty = 3, lwd = 2)
#     abline(v = log(lambda[id.best["GCV"]]), col = 6, lty = 3, lwd = 2)
#     abline(v = log(lambda[id.best["GIC"]]), col = 7, lty = 3, lwd = 2)
#     legend('topright', legend = c('min.AIC', 'min.AICc', 'min.BIC', 'min.eBIC', 'min.GCV', 'min.GIC'), lty = 3, col = 2:7, lwd = 1.5, bty = "n")
#   }
#   if (yvar == "se") {
#     se1 <- object$SE
#     if(intercept) se1 <- se1[, -1]
#     matplot(log(lambda), se1, type = "l", xlab = "log Lambda", ylab = "Std.Errs", xlim = range(log(lambda)) + c(-label,0 ), ...)
#     abline(h = 0, lty = 1, col = 'gray')
#     if (label) text(x = min(log(lambda)) - label, y = se1[1, ], labels = colnames(coef1), cex = cex.lab, xpd = TRUE, pos = 4)
#     abline(v = log(lambda[id.best["AIC"]]), col = 2, lty = 3, lwd = 2)
#     abline(v = log(lambda[id.best["AICc"]]), col = 3, lty = 3, lwd = 2)
#     abline(v = log(lambda[id.best["BIC"]]), col = 4, lty = 3, lwd = 2)
#     abline(v = log(lambda[id.best["eBIC"]]), col = 5, lty = 3, lwd = 2)
#     abline(v = log(lambda[id.best["GCV"]]), col = 6, lty = 3, lwd = 2)
#     abline(v = log(lambda[id.best["GIC"]]), col = 7, lty = 3, lwd = 2)
#     legend('topright', legend = c('min.AIC', 'min.AICc', 'min.BIC', 'min.eBIC', 'min.GCV', 'min.GIC'), lty = 3, col = 2:7, lwd = 1.5, bty = "n")
#   }
#   if (yvar == "weight") {
#     weight1 <- object$Weight
#     if(intercept) weight1 <- weight1[, -1]
#     matplot(log(lambda), weight1, type = "l", xlab = "log Lambda", ylab = "Mixture weight", xlim = range(log(lambda)) + c(-label, 0), ...)
#     abline(h = 0, lty = 1, col = 'gray')
#     if (label) text(x = min(log(lambda)) - label, y = weight1[1, ], labels = colnames(coef1), cex = cex.lab, xpd = TRUE, pos = 4)
#     abline(v = log(lambda[id.best["AIC"]]), col = 2, lty = 3, lwd = 2)
#     abline(v = log(lambda[id.best["AICc"]]), col = 3, lty = 3, lwd = 2)
#     abline(v = log(lambda[id.best["BIC"]]), col = 4, lty = 3, lwd = 2)
#     abline(v = log(lambda[id.best["eBIC"]]), col = 5, lty = 3, lwd = 2)
#     abline(v = log(lambda[id.best["GCV"]]), col = 6, lty = 3, lwd = 2)
#     abline(v = log(lambda[id.best["GIC"]]), col = 7, lty = 3, lwd = 2)
#     legend('topright', legend = c('min.AIC', 'min.AICc', 'min.BIC', 'min.eBIC', 'min.GCV', 'min.GIC'), lty = 3, col = 2:7, lwd = 1.5, bty = "n")
#   }
#   if (yvar == "gradient") {
#     fn0 <- function(x, y, b, s = 1, c = .5, lambda, alpha, unpenalized, family, offset, weights) {
#       eta <- x %*% b + offset
#       mu <- family$linkinv(eta)
#       v <- family$variance(mu)
#       m.e <- family$mu.eta(eta)
#       r <- y - mu
#       rv <- r / v * m.e * weights
#       grad <- drop(t(rv) %*% x)
#       
#       bsc <- (b / s)
#       # grad <- crossprod(x, y - drop(x %*% b))
#       r <- alpha * (c*(2 * pnorm(bsc, 0, 1) - 1) + (1 - c)*(2 * pnorm(bsc, 0, 1E-5) - 1)) + (1 - alpha) * b
#       if(any(unpenalized)) r[unpenalized] <- 0
#       return(- grad + lambda * r)
#     }
#     grad <- t(sapply(seq_len(nlambda), function(i) {
#       fn0(x = object$Input$x, y = object$Input$y, b = object$Coef[i, ], 
#           s = object$SE[i, ], c = object$Weight[i, ], 
#           lambda = lambda[i], alpha = object$Input$alpha, 
#           unpenalized = object$Input$unpenalized, family = object$Input$family, 
#           offset = object$Input$offset, weights = object$Input$weights)
#     }))
#     if(intercept) grad <- grad[, -1]
#     matplot(lambda, grad, type = "l", xlab = "Lambda", ylab = "Gradient", xlim = range(lambda) + c(0, label), ...)
#     abline(h = 0, lty = 1, col = 'gray')
#     if (label) text(x = max(lambda), y = grad[nlambda, ], labels = colnames(coef1), cex = cex.lab, xpd = TRUE, pos = 4)
#     # abline(v = lambda[id.best["AIC"]], col = 2, lty = 2, lwd = 2)
#     # abline(v = lambda[id.best["AICc"]], col = 3, lty = 3, lwd = 2)
#     # abline(v = lambda[id.best["BIC"]], col = 4, lty = 4, lwd = 2)
#     # abline(v = lambda[id.best["GCV"]], col = 5, lty = 5, lwd = 2)
#     # legend('topleft', legend = c('min.AIC', 'min.AICc', 'min.BIC', 'min.GCV'), lty = 2:5, col = 2:5, lwd = 1)
#   }
#   
#   invisible()
# }

predict.islasso.path <- function(object, newdata, 
                                 type = c("link", "response", "coefficients", "class"), lambda, ...){
  type <- match.arg(type)
  family <- object$family
  if(type == "class" & family$family != "binomial") 
    stop(gettextf("Type 'class' is available only for the %s family", sQuote(binomial()$family)), domain=NA)
  
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
  lambda.seq <- object$Info[,"lambda"]
  nlambda <- length(lambda.seq)
  if(missing(lambda) | nlambda == 1){
    beta <- object$Coef
    predictor <- tcrossprod(X, beta)
    lambda <- lambda.seq
    nlmb <- nlambda
  }
  else{
    if(any(lambda < min(lambda.seq)) | any(lambda > max(lambda.seq))) stop("value of lambda out of bound")
    lambda <- sort(lambda)
    nlmb <- length(lambda)
    beta <- NULL
    for(i in seq_len(nlmb)){
      id1 <- rev(which(lambda[i] >= lambda.seq))[1]
      id2 <- which(lambda[i] < lambda.seq)[1]
      beta <- cbind(beta, interpolate(object$Coef[id1,], object$Coef[id2,], lambda.seq[id1], lambda.seq[id2], lambda[i]))
    }
    predictor <- drop(X %*% beta)
  }
  if (!is.null(offset)) predictor <- predictor + offset
  
  out <- switch(type,
                "link" = { predictor },
                "response" = { family$linkinv(predictor) },
                "coefficients" = { t(beta) },
                "class" = {
                  mu <- family$linkinv(predictor)
                  1*I(mu >= .5)} )
  if(is.matrix(out)) colnames(out) <- paste0("lambda = ", formatC(lambda, format = "f", digits = 3))
  return(out)
}

coef.islasso.path <- function(object, lambda, ...){
  lambda.seq <- object$Info[,"lambda"]
  nlambda <- length(lambda.seq)
  if(nlambda == 1){
    coef <- object$Coef[1,]
  }
  else{
    if(missing(lambda)) return(object$Coef)
    if(lambda < min(lambda.seq) | lambda > max(lambda.seq)) stop("value of lambda out of bound")
    id1 <- rev(which(lambda >= lambda.seq))[1]
    id2 <- which(lambda < lambda.seq)[1]
    coef <- interpolate(object$Coef[id1,], object$Coef[id2,], lambda.seq[id1], lambda.seq[id2], lambda)  
  }
  
  return(coef)
  
}

fitted.islasso.path <- function(object, lambda, ...){
  lambda.seq <- object$Info[,"lambda"]
  nlambda <- length(lambda.seq)
  if(nlambda == 1){
    fit <- object$Fitted.values[1,]
  }
  else{
    if(missing(lambda)) return(object$Fitted.values)
    if(lambda < min(lambda.seq) | lambda > max(lambda.seq)) stop("value of lambda out of bound")
    id1 <- rev(which(lambda >= lambda.seq))[1]
    id2 <- which(lambda < lambda.seq)[1]
    fit <- interpolate(object$Fitted.values[id1,], object$Fitted.values[id2,], lambda.seq[id1], lambda.seq[id2], lambda)
  }
  return(fit)
}

residuals.islasso.path <- function(object, type = c("working", "response", "deviance", "pearson"), lambda, ...){
  lambda.seq <- object$Info[,"lambda"]
  nlambda <- length(lambda.seq)
  if(nlambda == 1){
    res <- object$Residuals[1,]
    mu <- object$Fitted.values[1,]
  }
  else{
    if(missing(lambda)) return(object$Residuals)
    if(lambda < min(lambda.seq) | lambda > max(lambda.seq)) stop("value of lambda out of bound")
    id1 <- rev(which(lambda >= lambda.seq))[1]
    id2 <- which(lambda < lambda.seq)[1]
    res <- interpolate(object$Residuals[id1,], object$Residuals[id2,], lambda.seq[id1], lambda.seq[id2], lambda)
    mu <- interpolate(object$Fitted.values[id1,], object$Fitted.values[id2,], lambda.seq[id1], lambda.seq[id2], lambda)
  }
  type <- match.arg(type)
  y <- object$y
  weights <- object$priorweights
  family <- object$family
  
  switch(type,
         "working" = return(res),
         "response" = return(y - mu),
         "deviance"=return(sign(y - mu) * sqrt(family$dev.resids(y, mu, weights))),
         "pearson"=return((y - mu) * sqrt(weights) / sqrt(family$variance(mu))))
}

deviance.islasso.path <- function(object, lambda, ...){
  lambda.seq <- object$Info[,"lambda"]
  nlambda <- length(lambda.seq)
  if(nlambda == 1){
    dev <- object$Info[1, "deviance"][[1]]
  }
  else{
    if(missing(lambda)) return(object$Info[, "deviance"])
    if(lambda < min(lambda.seq) | lambda > max(lambda.seq)) stop("value of lambda out of bound")
    id1 <- rev(which(lambda >= lambda.seq))[1]
    id2 <- which(lambda < lambda.seq)[1]
    dev <- interpolate(object$Info[id1, "deviance"], object$Info[id2, "deviance"], lambda.seq[id1], lambda.seq[id2], lambda)[[1]]
  }
  return(dev)
}

logLik.islasso.path <- function(object, lambda, ...){
  lambda.seq <- object$Info[,"lambda"]
  nlambda <- length(lambda.seq)
  if(nlambda == 1){
    ll <- object$Info[1, "logLik"]
    df <- object$Info[1, "df"]
    phi <- object$Info[1, "phi"]
  }
  else{
    if(missing(lambda)) return(object$Info[, "logLik"])
    if(lambda < min(lambda.seq) | lambda > max(lambda.seq)) stop("value of lambda out of bound")
    id1 <- rev(which(lambda >= lambda.seq))[1]
    id2 <- which(lambda < lambda.seq)[1]
    ll <- interpolate(object$Info[id1, "logLik"], object$Info[id2, "logLik"], lambda.seq[id1], lambda.seq[id2], lambda)
    df <- interpolate(object$Info[id1, "df"], object$Info[id2, "df"], lambda.seq[id1], lambda.seq[id2], lambda)
    phi <- interpolate(object$Info[id1, "phi"], object$Info[id2, "phi"], lambda.seq[id1], lambda.seq[id2], lambda)
  }
  out <- list(loglik = ll, df = df, object = object, phi = object$phi)
  class(out) <- "logLik.islasso.path"
  out
}

print.logLik.islasso.path <- function(x, digits=max(3L, getOption("digits") - 3L), ...){
  cat("\n'log Lik.' ", paste(format(c(x$loglik), digits = digits), collapse = ", "),
      " (df = ", format(x$df), ")\n\n", sep = "")
  invisible(x)
}
