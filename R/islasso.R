islasso <- function(formula, family=gaussian, lambda, alpha=1, data, weights, subset, offset, unpenalized, contrasts = NULL, control = is.control()){
  this.call <- match.call()

  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "offset"), names(mf), 0L)
  ioff <- if(m[5] != 0) TRUE else FALSE
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  if(ioff) off <- mf$offset
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "any")
  X <- if(!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else stop("Model matrix is empty!")
  temp <- which(attr(mt, "dataClasses")[names(attr(mt, "dataClasses")) %in% attr(mt, "term.labels")] %in% c("factor", "character"))
  temp <- which(attr(X, "assign") %in% temp)

  if(ioff){
    noff <- match(unlist(lapply(off, as.character)), colnames(X))
    if(!all(is.na(noff))) X <- X[, -noff[which(noff!=0)]]
  }
  offset <- as.vector(model.offset(mf))
  if(!is.null(offset)){
    if(length(offset) != NROW(y)) stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                                                length(offset), NROW(y)), domain = NA)
  }
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0))
    stop("negative weights not allowed")
  if(alpha > 1 | alpha < 0)
    stop("alpha must be in [0, 1] (0 for ridge penalty, 1 for lasso penalty)")
  
  if(attr(mt,"intercept") == 0){
    intercept <- FALSE
  }else{
    intercept <- TRUE
    X <- X[, -1, drop = FALSE]
  }
  attributes(X)$dataClasses <- temp
  setting <- control
  if(missing(unpenalized)) unpenalized <- NULL

  # fit <- switch(setting$algorithm,
  #               "smoothing"=islasso.fit(X=X, y=y, family=family, lambda=lambda, intercept=intercept, weights=weights, offset=offset, unpenalized=unpenalized, control=control))#,
  #               "hyperbolic"=hlasso.fit(X=X, y=y, family=family, lambda=lambda, intercept=intercept, weights=weights, offset=offset, unpenalized=unpenalized, control=control))
  
  fit <- islasso.fit(X=X, y=y, family=family, lambda=lambda, alpha=alpha, intercept=intercept, weights=weights, offset=offset, unpenalized=unpenalized, control=control)
  fit$offset <- offset
  fit$formula <- formula
  fit$call <- this.call
  fit$model <- mf
  fit$terms <- mt
  fit$data <- data
  fit$contrasts <- contrasts
  fit$xlevels <- .getXlevels(mt, mf)
  class(fit) <- "islasso"

  return(fit)
}

model.matrix.islasso <- function (object, ...) {
  model.matrix(object$terms, object$model, object$contrasts, ...)
}

is.control <- function(sigma2 = -1, tol = 1e-04, itmax = 500, stand = TRUE,
                       trace = 0, nfolds = 5, seed=NULL, debias = FALSE, adaptive = FALSE, 
                       b0 = NULL, V0 = NULL, c = -1){
  
  list(sigma2=sigma2, tol=tol, itmax=itmax, trace=trace,
       debias=debias, stand=stand, nfolds=nfolds, seed=seed,
       adaptive=adaptive, b0=b0, V0=V0, c=c)
}

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

  if(is.null(temp$unbias)) temp$unbias <- FALSE
  if(temp$unbias & is.null(object$beta.unbias)) stop("No unbias estimates in the model")

  coef <- if(!temp$unbias) object$coef else object$beta.unbias
  se <- if(!temp$unbias) object$se else object$se.unbias
  aic <- object$aic
  n <- object$internal$n
  p <- object$internal$p

  res <- residuals(object, type="deviance")
  dispersion <- object$phi
  rdf <- object$internal$n - object$rank
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
              family=object$family, iter=object$internal$iter, pval=pval, temp=temp, call=object$call)

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
  temp <- coefs[, 5] <= x$pval
  if(sum(temp) != 0) coefs <- coefs[temp, ]
  if(sum(temp) == 1){
    coefs <- t(coefs)
    rownames(coefs) <- rownames(x$coefficients)[temp]
  }
  printCoefmat(coefs, digits=digits, signif.stars=TRUE, has.Pvalue=TRUE, na.print="NA", cs.ind = 1L:2L, tst.ind = 3L:4L,, ...)
  cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ", format(x$dispersion), ")\n\n",
      apply(cbind(paste(format(c("Null", "Residual"), justify = "right"), "deviance:"),
                  format(c(x$nulldev, x$dev), digits = max(5L, digits + 1L)), " on",
                  format(c(x$df.null, x$df.res), digits = digits), " degrees of freedom\n"), 1L, paste, collapse = " "), sep = "")
  cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)), "\nLambda: ", format(x$lambda, digits = max(4L, digits + 1L)), "\n\n", "Number of Newton-Raphson iterations: ", x$iter, "\n", sep = "")

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
  sdev <- dev/sqrt(x$phi)

  pea <- residuals(x, type="pearson")
  spea <- pea/sqrt(x$phi)

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

format.perc <- function (probs, digits){
  paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), "%")
}

predict.islasso <- function(object, type=c("link", "response", "coefficients", "class"), newdata, ...){
  beta <- as.matrix(coef(object))
  intercept <- object$internal$intercept
  type <- match.arg(type)
  family <- object$family

  if(missing(newdata)){
    X <- model.matrix(object)
  }else{
    X <- as.matrix(newdata)
    if(!all(colnames(X) %in% rownames(beta))) stop("'newdata' must contain all X-variables")
    if(intercept) X <- cbind(1, X)
  }

  if(nrow(beta) != ncol(X)) stop("the length of the vector 'beta' is not equal to the number of columns of the matrix 'X'")
  if(type == "class" & family$family != "binomial") stop(gettextf("Type 'class' is available only for the %s family", sQuote(family$family)), domain=NA)

  out <- switch(type,
                "link"={drop(X %*% beta)},
                "response"={family$linkinv(drop(X %*% beta))},
                "coefficients"={drop(beta)},
                "class"={
                  mu <- family$linkinv(drop(X %*% beta))
                  1*I(mu >= .5)})

  return(out)
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

islasso.fit <- function(X, y, family=gaussian, lambda, alpha=1, intercept=FALSE, weights=NULL, 
                        offset=NULL, unpenalized=NULL, control=is.control()){
  this.call <- match.call()

  X <- as.matrix(X)
  nX <- if(intercept) cbind(1, X) else X
  yy <- y
  nobs <- nrow(nX)
  nvars <- ncol(nX)
  nms <- colnames(X)
  if(is.null(nms) & ncol(X) != 0) nms <- paste0("X", 1:ncol(X))
  if(intercept) nms <- c("(Intercept)", nms)
  colnames(nX) <- nms

  if(alpha > 1 | alpha < 0)
    stop("alpha must be in [0, 1] (0 for ridge penalty, 1 for lasso penalty)")
  if(is.null(unpenalized)){
    unpenalized <- rep(FALSE, (nvars-1*intercept))
    if(intercept) unpenalized <- c(TRUE, unpenalized)
  }else{
    if(!is.vector(unpenalized)) stop("''unpenalized is not a vector")
    if(is.list(unpenalized)) stop("'unpenalized' can not be a list")
    if (is.factor(unpenalized)) stop("'unpenalized' can not be a factor")
    if(is.character(unpenalized)){
      temp_nms <- if(intercept) nms[-1] else nms
      unpenalized_id <- pmatch(unpenalized, temp_nms)
      if(any(is.na(unpenalized_id))) stop(gettextf("the following names are not in colnames(X): %s", paste(unpenalized[is.na(unpenalized_id)], collapse = ", ")))
      unpenalized <- sort(unpenalized_id)
    }else unpenalized <- sort(unpenalized)
    if(any(abs(unpenalized - round(unpenalized)) > .Machine$double.eps^0.5)) stop("some element of 'unpenalized' is not an integers")
    if(any(unpenalized <= 0)) stop("some element of 'unpenalized' is smaller than zero")
    if(any(unpenalized > (nvars-1*intercept))) stop("some element of 'unpenalized' is greater than the number of columns of the matrix 'X'")
    temp <- rep(FALSE, (nvars-1*intercept))
    temp[unpenalized] <- TRUE
    if(intercept) temp <- c(TRUE, temp)
    unpenalized <- temp
  }

  n <- rep(1, nobs)
  if(is.null(offset)) offset <- rep(0, nobs)
  if(is.null(weights)){
    weights <- as.double(rep(1, nobs))
  }else{
    if(!is.vector(weights)) stop("argument 'weights' is not a vector")
    if(is.list(weights)) stop("argument 'weights' can not be a list")
    if(is.factor(weights)) stop("argument 'weights' can not be a factor")
    if(is.character(weights)) stop("vector 'weights' can not be a character")
    if(length(weights) != nobs) stop("the length of the vector 'weights' is not equal to ", sQuote(n))
    if(any(is.nan(weights))) weights[is.nan(weights)] <- 0
    if(all(weights == 0)) stop("all the entries of the vector 'weights' are equal to zero")
  }
  if(is.character(family)) family <- get(family, mode="function", envir=parent.frame())
  if(is.function(family)) family <- do.call(family, args=list(), envir=parent.frame())
  if(is.null(family$family)){
    print(family)
    stop("'family' not recognized!")
  }
  tempFamily <- family$family
  tempLink <- family$link

  setting <- control
  if(setting$nfolds < 3) stop("'nfolds' should greater than 3")
  if(setting$tol <= 0) stop("'tol' should be a non-negative value")
  if(setting$itmax < 0) stop("'itmax' should be a non-negative value")
  if(setting$c[1] > 1) stop("'pai' should be fixed in (0,1) or estimated using a negative value")
  if(setting$c[1] < 0){
    estpai <- TRUE
    c <- 1
  }else{
    estpai <- FALSE
  }
  setting$c <- rep(setting$c, nvars)
  if(tempFamily == "binomial") setting$sigma2 <- 1
  if(tempFamily == "quasibinomial") setting$sigma2 <- -1
  if(tempFamily == "poisson") setting$sigma2 <- 1
  if(tempFamily == "quasipoisson") setting$sigma2 <- -1

  okLinks <- c("identity", "logit", "log", "probit", "inverse")
  indLink <- which(okLinks %in% tempLink)
  if(!(tempLink %in% okLinks)) stop(gettextf("%s link not recognized", sQuote(tempLink)), domain = NA)
  if(tempFamily == "gaussian" & !is.element(tempLink, c("identity")))
    stop(gettextf("The %s family does not accept the link %s.\n  The accepted link is: 'identity'", sQuote(tempFamily), sQuote(tempLink)), domain=NA)
  if(any(tempFamily == c("poisson", "quasipoisson")) & !is.element(tempLink, c("log")))
    stop(gettextf("The %s family does not accept the link %s.\n  The accepted link is: 'log'", sQuote(tempFamily), sQuote(tempLink)), domain = NA)
  if(any(tempFamily == c("binomial", "quasibinomial")) & !is.element(tempLink, c("logit", "probit")))
    stop(gettextf("The %s family does not accept the link %s.\n  The accepted link are: 'logit' and 'probit'", sQuote(tempFamily), sQuote(tempLink)), domain=NA)
  if(tempFamily == "Gamma" & !is.element(tempLink, c("log", "identity", "inverse")))
    stop(gettextf("The %s family does not accept the link %s.\n  The accepted link are: 'log', 'inverse', 'identity'", sQuote(tempFamily), sQuote(tempLink)), domain=NA)
  if(is.character(y) & all(tempFamily != c("binomial", "quasibinomial")))
    stop(gettextf("The %s family does not accept a character as responce variable", sQuote(tempFamily)), domain=NA)
  if(is.factor(y) & all(tempFamily != c("binomial", "quasibinomial")))
    stop(gettextf("The %s family does not accept a factor as responce variable", sQuote(tempFamily)), domain=NA)
  if(is.matrix(y) & all(tempFamily != c("binomial", "quasibinomial")))
    stop(gettextf("The %s family does not accept a matrix as responce variable", sQuote(tempFamily)), domain=NA)
  if(any(tempFamily == c("binomial", "quasibinomial"))){
    if(NROW(y) != nobs) stop("the length of the vector 'y' is not equal to the number of rows of the matrix 'X'")
    if(NCOL(y) == 1){
      if(is.character(y)) y <- factor(y)
      if(is.factor(y)){
        if(nlevels(y) != 2) stop("'y' can not be a factor with more than two levels")
        y <- as.numeric(y) - 1
      }
      if(any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")
      mustart <- (weights * y + 0.5)/(weights + 1)
    }else{
      if(NCOL(y) == 2){
        if(is.data.frame(y)) y <- as.matrix(y)
        if(is.character(y)) stop("'y' can not be a character matrix")
        if(any(abs(y - round(y)) > 0.001)) warning("non-integer counts in a binomial glm!")
        n <- y[, 1] + y[, 2]
        y <- ifelse(n == 0, 0, y[, 1]/n)
        weights <- weights * n
        mustart <- (n * y + 0.5)/(n + 1)
      }else{
        stop("for the 'binomial' family, y must be a vector of 0 and 1's\nor a 2 column matrix where col 1 is no. successes and col 2 is no. failures")
      }
    }
  }
  if(tempFamily %in% c("poisson", "quasipoisson") & any(y < 0))
    stop(gettextf("each element of 'y' must be positive when family %s is used", sQuote(tempFamily)), domain=NA)
  if(tempFamily %in% c("Gamma") & any(y <= 0))
    stop(gettextf("each element of 'y' must be positive when family %s is used", sQuote(tempFamily)), domain=NA)
  if(tempFamily == "gaussian") mustart <- y
  if(any(tempFamily == c("poisson", "quasipoisson"))) mustart <- y + 0.1
  if(tempFamily == "Gamma") mustart <- y

  variance <- family$variance
  linkinv <- family$linkinv
  dev.resids <- family$dev.resids
  aic <- family$aic
  mu.eta <- family$mu.eta
  linkfun <- family$linkfun

  if(!missing(lambda)){
    if(is.character(lambda)) stop("\n'lambda' can not be a character\n")
    if(is.factor(lambda)) stop("\n'lambda' can not be a factor\n")
    if(lambda < 0) stop("\n'lambda' is negative\n")
  }

  if(is.null(setting$b0)){
    tempFamily2 <- switch(tempFamily,
                          "gaussian"="gaussian",
                          "binomial"="binomial",
                          "quasibinomial"="binomial",
                          "poisson"="poisson",
                          "quasipoisson"="poisson",
                          "Gamma"="Gamma")
    
    if((nvars - intercept) < 2 | tempFamily2 == "Gamma"){
      if(missing(lambda)) stop("Insert a positive value for lambda")
      # if(missing(lambda)){lambda <- sqrt(log(nvars)/nobs)}
      est <- rep(0.1, nvars)
    }
    else{
      if(missing(lambda)){
        type.measure <- switch(tempFamily,
                               "gaussian"="mse",
                               "binomial"="class",
                               "quasibinomial"="class",
                               "poisson"="deviance",
                               "quasipoisson"="deviance")
        obj <- cv.glmnet(x=as.matrix(X), y=if(is.matrix(yy)) yy[, c(2,1)] else yy, family=tempFamily2, 
                         nfolds=setting$nfolds, type.measure=type.measure, standardize=setting$stand, 
                         intercept=intercept, offset=offset, alpha=alpha)
        lambda <- obj$lambda.min*nobs
        est <- as.vector(coef(obj, s="lambda.min"))
      }else{
        obj <- glmnet(x=as.matrix(X), y=if(is.matrix(yy)) yy[, c(2,1)] else yy, family=tempFamily2, alpha=alpha, 
                      standardize=setting$stand, intercept=intercept, offset=offset)
        est <- as.vector(coef(obj, s=lambda/nobs))
      }
      if(!intercept) est <- est[-1]
    }
  }
  else{
    if(missing(lambda)) stop("Insert a positive value for lambda")
    est <- setting$b0
  }
  interval <- if(exists("obj")) range(rev(obj$lambda*nobs)) else NULL

  Lambda <- rep(lambda, nvars)
  Lambda[unpenalized] <- 0

  ntheta <- est
  Gamcov1 <- if(is.null(setting$V0)) diag(1, nvars) else setting$V0
  se1 <- sqrt(diag(Gamcov1))
  eta <- linkfun(mustart) + offset
  mu <- linkinv(eta)
  res <- y - mu

  fam <-  c("binomial","poisson","Gamma")
  fam2 <-  c("quasibinomial","quasipoisson","Gamma")
  if((tempFamily %in% fam) | tempFamily %in% fam2){
    fam <- if(tempFamily %in% fam) pmatch(tempFamily, fam) else pmatch(tempFamily, fam2)
    link <- switch(fam,
                   "1"=pmatch(tempLink, c("logit","probit")),
                   "2"=pmatch(tempLink, c("log")),
                   "3"=pmatch(tempLink, c("inverse","log","identity")))
  }
  else{
    fam <- 0
    link <- 0
  }

  storage.mode(nX) <- "double"
  storage.mode(y) <- "double"
  storage.mode(nobs) <- "integer"
  storage.mode(nvars) <- "integer"
  storage.mode(ntheta) <- "double"
  storage.mode(se1) <- "double"
  storage.mode(Gamcov1) <- "double"
  storage.mode(Lambda) <- "double"
  storage.mode(alpha) <- "double"
  storage.mode(setting$c) <- "double"
  h <- as.double(1)
  storage.mode(setting$itmax) <- "integer"
  storage.mode(setting$tol) <- "double"
  storage.mode(setting$sigma2) <- "double"
  storage.mode(offset) <- "double"
  storage.mode(eta) <- "double"
  storage.mode(mu) <- "double"
  storage.mode(res) <- "double"
  storage.mode(weights) <- "double"

  fit <- if(tempFamily == "gaussian"){
    .Fortran(C_islasso2, X=nX, y=y, n=nobs, p=nvars, ntheta=ntheta, se=se1, cov=Gamcov1, lambda=Lambda, alpha=alpha, pi=setting$c, estpi=as.integer(estpai),
             h=h, itmax=setting$itmax, tol=setting$tol, phi=setting$sigma2, trace=as.integer(setting$trace), adaptive=as.integer(setting$adaptive),
             offset=offset, conv=integer(1), stand=as.integer(setting$stand), intercept=as.integer(intercept), eta=eta, mu=mu, res=res, dev=double(1),
             weights=weights, hi=double(nvars), edf=double(1), grad=double(nvars))
  }else{
    .Fortran(C_islasso_glm, X=nX, y=y, n=nobs, p=nvars, ntheta=ntheta, se=se1, cov=Gamcov1, lambda=Lambda, alpha=alpha, pi=setting$c, estpi=as.integer(estpai),
             h=h, itmax=setting$itmax, tol=setting$tol, phi=setting$sigma2, trace=as.integer(setting$trace), adaptive=as.integer(setting$adaptive),
             offset=offset, conv=integer(1), stand=as.integer(setting$stand), intercept=as.integer(intercept), eta=eta, mu=mu, dev=double(1),
             weights=weights, hi=double(nvars), edf=double(1), fam=as.integer(fam), link=as.integer(link), grad=double(nvars))
  }
  if(fit$conv == 1) warning("Maximum number of iterations attained!!")

  setting$c <- fit$pi
  ntheta <- fit$ntheta
  se1 <- fit$se
  Gamcov1 <- fit$cov
  dev <- fit$dev
  s2 <- fit$phi
  eta <- drop(nX %*% ntheta + offset)
  mu <- linkinv(eta)
  varmu <- variance(mu)
  mu.eta.val <- mu.eta(eta)
  w <- ((weights*mu.eta.val^2)/varmu)
  W <- diag(w)
  res <- (y - mu)/mu.eta.val
  XtW <- crossprod(nX, W)
  XtX <- XtW %*% nX
  A <- .Fortran(C_hessian, ntheta, se1, Lambda, XtX, setting$c, as.integer(nvars), hess=XtX, alpha)$hess
  invH <- .Fortran(C_inv3, as.integer(nvars), A, invA=A, integer(1))$invA
  gradient <- fit$grad #.Fortran(C_gradient, ntheta, se1, Lambda, XtW, res, setting$c, as.integer(nobs), as.integer(nvars), grad=ntheta, alpha)$grad

  hi <- fit$hi
  edf <- fit$edf

  wtdmu <- if(intercept) sum(weights * y)/sum(weights) else linkinv(offset)
  nulldev <- sum(dev.resids(y, wtdmu, weights))
  n.ok <- nobs - sum(weights == 0)
  nulldf <- n.ok - as.integer(intercept)
  rank <- edf
  resdf <- n.ok - rank
  aic.model <- aic(y, n, mu, weights, dev) + 2 * rank

  iter <- fit$itmax
  conv <- fit$conv

  se.unbias <- beta.unbias <- NULL
  if(setting$debias){
   invH2 <- ginv2(XtX) #.Fortran(C_inv2, as.integer(nvars), XtX, invA=XtX, integer(0))$invA
   # invH2 <- .Fortran(C_inv3, as.integer(nvars), XtX, invA=XtX, integer(1))$invA
   gradient2 <- XtW %*% res
   delta <- invH2 %*% gradient2
   beta.unbias <- drop(ntheta + delta)
   se.unbias <-  sqrt(diag(s2*(invH2)))
   names(beta.unbias) <- names(se.unbias) <- nms
  }

  names(gradient) <- names(se1) <- names(ntheta) <- colnames(XtX) <- rownames(XtX) <- colnames(Gamcov1) <- rownames(Gamcov1) <- nms

  internal <- list(model=NULL, terms=NULL, y.internal=y, offset=offset, n=nobs, p=nvars, weights=weights, wt=w, lambda.seq=fit$lambda,
                   XtW=XtW, res=res, XtX=XtX, invH=invH, vcov=Gamcov1, gradient=gradient, hessian=A, hi=hi, intercept=intercept, unpenalized=unpenalized,
                   fam=fam, link=link, nms=nms, estc=estpai, lmbd.interval=interval, iter=iter, conv=conv)
  out <- list(coefficients=ntheta, se=se1, residuals=res, fitted.values=mu, linear.predictors=eta, rank=edf, family=family, lambda=lambda, alpha=alpha,
              deviance=dev, null.deviance=nulldev, aic=aic.model, df.null=nulldf, phi=s2, beta.unbias=beta.unbias, se.unbias=se.unbias,
              internal=internal, control=setting, call=this.call, formula=NULL)

  return(out)
}

ginv2 <- function (X, tol = sqrt(.Machine$double.eps)){
  if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
    stop("'X' must be a numeric or complex matrix")
  if (!is.matrix(X))
    X <- as.matrix(X)
  Xsvd <- svd(X)
  if (is.complex(X))
    Xsvd$u <- Conj(Xsvd$u)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  if (all(Positive))
    Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive))
    array(0, dim(X)[2L:1L])
  else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
}

qqNorm <- function(x, probs=seq(.005, .995, l=200), centre=FALSE, scale=FALSE, leg=TRUE,
                   mean=0, sd=1, add=FALSE, p.col=1, dF=FALSE, ylab, ...){
  #a simple qqNorm function to compare an emprical distribution with respect to
  # a Normal distribution with given mean and sd
  #x: the data vector
  #probs: percentiles wrt compute the quantiles to be plotted
  #centre: if TRUE the observed data are centered (to have a zero mean)
  #scale: if TRUE the observed data are scaled (to have a unit variance)
  #leg: if TRUE the legend with some information is added
  #mean, sd: mean and st.dev of the theoric Normal distribution
  #add: if TRUE the QQ plot is added (provided graphical device is open..)
  #p.col: the color of dots
  #dF: se TRUE the distribution function (rather than the QQplot) is plotted
  #author: vito.muggeo@unipa.it
  if(centre) x <- x - mean(x)
  if(scale){
    x <- x/sd(x) + mean(x)*(sd(x)-1)/sd(x)
  }
  emp <- quantile(x, probs, names=FALSE)
  emp <- c(min(x), emp, max(x))

  teor <- qnorm(probs, mean=mean, sd=sd)
  teor <- qnorm(c(.0005, probs, .9995), mean=mean, sd=sd)
  if(dF){
    val <- seq(min(x), max(x), l=200)
    plot(sort(x), 1:length(x)/length(x), type="s", xlab=deparse(substitute(x)), ylab="Distribution Function",  col=p.col,...)
    lines(val, pnorm(val, mean=mean, sd=sd), col=2)
  }else{
    if(add){
      points(teor, emp, pch=19, col=p.col)
      return(invisible(NULL))
    }
    if(missing(ylab)) ylab <- "Empirical quantiles"
    plot(teor, emp, xlab="Theoretical quantiles",
         ylab=ylab, pch=19, col=p.col,
         xlim=range(teor), ylim=range(emp), ...)
    abline(0, 1, col=2, lwd=2)
    #	abline(h=mean, v=mean, col=2, lty=3, lwd=1.5)
    segments(par()$usr[1], mean(x), max(teor[emp<=mean(x)]), mean(x), col=1, lty=2)
    segments(par()$usr[1], mean, max(teor[emp<=mean(x)]), mean, col=2, lty=3, lwd=1.5)
    segments(mean, par()$usr[3], mean, mean(x), col=2, lty=3, lwd=1.5)
  }
  if(leg){
    et <- c(paste("emp. mean=", round(mean(x), 3), sep=""), paste("emp. sd=", round(sd(x), 3), sep=""))
    legend("topleft", et, bty="n", cex=.7 )
    etT <- c(paste("theor. mean=", round(mean, 3), sep=""), paste("theor. sd=", round(sd, 3), sep=""))
    legend("bottomright", etT, bty="n", cex=.7, text.col =2 )
  }
}

simulXy <- function(n, p, interc=0, beta, family = gaussian(), prop = 0.1, 
                   lim.b=c(-3,3), sigma = 1, size = 1, rho = 0, scale = TRUE, 
                   seed, X){
  if (!missing(seed)) set.seed(seed)
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- do.call(family, args = list(), envir = parent.frame())
  if (is.null(family$family) | !(family$family %in% c("gaussian", "binomial", "poisson"))) {
    #print(family)
    stop("'family' not recognized!")
  }
  
  if (missing(beta)) {
    if(prop <0 || prop>1) stop("invalid 'prop' ")
    p.true <- trunc(prop * p)
    beta <- c(runif(p.true, lim.b[1], lim.b[2]), rep(0, p - p.true))
  }
  if(missing(X)) {
    X <- modelX(n, p, rho, scale)
    colnames(X) <- paste0("X", 1:p)
  }
  eta <- drop(X %*% beta) + interc
  mu <- family$linkinv(eta)
  y <- switch(family$family, gaussian = rnorm(n, mu, sigma), 
              binomial = rbinom(n, size, mu), poisson = rpois(n, mu))
  if (family$family == "binomial" & size > 1) {
    y <- cbind(size - y, y)
    colnames(y) <- c("failure", "success")
  }
  data.frame(y = y, X)
}

modelX <- function(n, p, rho=.5, scale=TRUE){
  Sigma <- forceSymmetric(t(rho^outer(1:p,1:p,"-")))
  cSigma <- chol(Sigma)
  X <- as.matrix(replicate(p, rnorm(n)) %*% cSigma)
  if(scale) X <- scale(X)
  return(X)
}

anova.islasso <- function(object, A, b, ...){
  beta <- coef(object)
  nms <- names(beta)
  p <- length(beta)
  if(missing(A) & missing(b)) A <- diag(p)
  if(missing(A) & !missing(b)) stop("Constraint matrix A is missing")
  if(is.vector(A)) A <- t(A)
  k <- nrow(A)
  mpc <- if(k > 1) TRUE else FALSE
  if(!missing(A) & missing(b)) b <- rep(0, k)
  if(length(b) == 1 & k > 1) b <- rep(b, k)
  V <- vcov(object)
  Identity <- diag(k)
  AI <- cbind(A, Identity)
  betab <- c(beta, -b)
  beta_new <- AI %*% betab
  V_new <- tcrossprod((A %*% V), A)
  nomi <- character(length = k)
  if(mpc){
  #   for(i in seq_len(k)){
  #     id <- which(A[i, ] != 0)
  #     nomi[i] <- paste(A[i, id], nms[id], sep="*", collapse = " + ")
  #     if(nchar(nomi[i]) > 60) nomi[i] <- paste(strtrim(nomi[i], 60), "...")
  #   }
    statistic2 <- drop(crossprod(beta_new, solve(V_new, beta_new)))
    pval2 <- 1 - pchisq(statistic2, df = k)
  }else{
    statistic2 <- 0
    pval2 <- 0
  }
    pval <- statistic <- double(length = k)
    for(i in seq_len(k)){
      id <- which(A[i, ] != 0)
      nomi[i] <- paste(A[i, id], nms[id], sep="*", collapse = " + ")
      if(nchar(nomi[i]) > 60) nomi[i] <- paste(strtrim(nomi[i], 60), "...")
      statistic[i] <- (beta_new[i]^2) / V_new[1]
      pval[i] <- 1 - pchisq(statistic[i], df = 1)
    }
    if(mpc) nomi <- c(nomi, "Overall")
  # }
  
  object$anova <- list(A=A, b=b, coefficients=beta_new, vcov=V_new, k=k, nms=nomi, mpc=mpc, 
                       tstat=statistic, pvalues=pval, tstat2=statistic2, pvalues2=pval2)
  class(object) <- c("anova.islasso", class(object))
  return(object)
}

print.anova.islasso <- function(x, digits = max(3, getOption("digits") - 3), ...){
  cat("\n\t", "Simultaneous Tests for General Linear Hypotheses\n\n")
  call <- x$call
  if (!is.null(call)) {
    cat("Fit: ")
    print(call)
    cat("\n")
  }
  pq <- x$anova
  
  mtests <- cbind(pq$coefficients, sqrt(diag(pq$vcov)), pq$tstat, pq$pvalues)
  if(pq$mpc) mtests <- rbind(mtests, c(0, 0, pq$tstat2, pq$pvalues2))
  colnames(mtests) <- c("Estimate", "Std. Error", "chi2 value", "Pr(>|chi2|)")
  rownames(mtests) <- paste(pq$nms, "=", format(pq$b, digits=3, trim=TRUE))
  if(pq$mpc) rownames(mtests)[nrow(mtests)] <- pq$nms[nrow(mtests)]
  sig <- .Machine$double.eps
  cat("Linear Hypotheses:\n")
  if(pq$mpc){
    printCoefmat2(mtests, digits = digits, has.Pvalue = TRUE, P.values = TRUE, eps.Pvalue = sig)
  }else{
    printCoefmat(mtests, digits = digits, has.Pvalue = TRUE, P.values = TRUE, eps.Pvalue = sig)
  }
  # if(pq$mpc){
  #   cat("\nMultiple Comparison of Linear Hypotheses:\n")
  #   mtests <- cbind(pq$coefficients, sqrt(diag(pq$vcov)), pq$tstat2, pq$pvalues2)
  #   colnames(mtests) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  #   rownames(mtests) <- paste(pq$nms, "=", format(pq$b, digits=3, trim=TRUE))
  # # if(pq$mpc) 
  #   printCoefmat2(mtests, digits = digits, has.Pvalue = TRUE, P.values = TRUE, eps.Pvalue = sig)
  # # else
  # #   printCoefmat(mtests, digits = digits, has.Pvalue = TRUE, P.values = TRUE, eps.Pvalue = sig)
  # }
  cat("\n")
  invisible(x)
}

printCoefmat2 <- function (x, digits = max(3L, getOption("digits") - 2L), signif.stars = getOption("show.signif.stars"), 
                           signif.legend = signif.stars, dig.tst = max(1L, min(5L, digits - 1L)), cs.ind = 1:k, tst.ind = k + 1, 
                           zap.ind = integer(), P.values = NULL, 
                           has.Pvalue = nc >= 4L && length(cn <- colnames(x)) && substr(cn[nc], 1L, 3L) %in% c("Pr(", "p-v"), 
                           eps.Pvalue = .Machine$double.eps, na.print = "NA", quote = FALSE, right = TRUE, ...){
  if (is.null(d <- dim(x)) || length(d) != 2L) 
    stop("'x' must be coefficient matrix/data frame")
  nc <- d[2L]
  if (is.null(P.values)) {
    scp <- getOption("show.coef.Pvalues")
    if (!is.logical(scp) || is.na(scp)) {
      warning("option \"show.coef.Pvalues\" is invalid: assuming TRUE")
      scp <- TRUE
    }
    P.values <- has.Pvalue && scp
  }
  else if (P.values && !has.Pvalue) 
    stop("'P.values' is TRUE, but 'has.Pvalue' is not")
  if (has.Pvalue && !P.values) {
    d <- dim(xm <- data.matrix(x[, -nc, drop = FALSE]))
    nc <- nc - 1
    has.Pvalue <- FALSE
  }
  else xm <- data.matrix(x)
  k <- nc - has.Pvalue - (if (missing(tst.ind)) 
    1
    else length(tst.ind))
  if (!missing(cs.ind) && length(cs.ind) > k) 
    stop("wrong k / cs.ind")
  Cf <- array("", dim = d, dimnames = dimnames(xm))
  ok <- !(ina <- is.na(xm))
  for (i in zap.ind) xm[, i] <- zapsmall(xm[, i], digits)
  if (length(cs.ind)) {
    acs <- abs(coef.se <- xm[, cs.ind, drop = FALSE])
    if (any(ia <- is.finite(acs))) {
      digmin <- 1 + if (length(acs <- acs[ia & acs != 
                                          0])) 
        floor(log10(range(acs[acs != 0], finite = TRUE)))
      else 0
      Cf[, cs.ind] <- format(round(coef.se, max(1L, digits - 
                                                  digmin)), digits = digits)
    }
  }
  if (length(tst.ind)) 
    Cf[, tst.ind] <- format(round(xm[, tst.ind], digits = dig.tst), 
                            digits = digits)
  if (any(r.ind <- !((1L:nc) %in% c(cs.ind, tst.ind, if (has.Pvalue) nc)))) 
    for (i in which(r.ind)) Cf[, i] <- format(xm[, i], digits = digits)
  ok[, tst.ind] <- FALSE
  okP <- if (has.Pvalue) 
    ok[, -nc]
  else ok
  x1 <- Cf[okP]
  dec <- getOption("OutDec")
  if (dec != ".") 
    x1 <- chartr(dec, ".", x1)
  x0 <- (xm[okP] == 0) != (as.numeric(x1) == 0)
  if (length(not.both.0 <- which(x0 & !is.na(x0)))) {
    Cf[okP][not.both.0] <- format(xm[okP][not.both.0], digits = max(1L, 
                                                                    digits - 1L))
  }
  if (any(ina)) 
    Cf[ina] <- na.print
  if (P.values) {
    if (!is.logical(signif.stars) || is.na(signif.stars)) {
      warning("option \"show.signif.stars\" is invalid: assuming TRUE")
      signif.stars <- TRUE
    }
    if (any(okP <- ok[, nc])) {
      pv <- as.vector(xm[, nc])
      Cf[okP, nc] <- format.pval(pv[okP], digits = dig.tst, 
                                 eps = eps.Pvalue)
      signif.stars <- signif.stars && any(pv[okP] < 0.1)
      if (signif.stars) {
        Signif <- symnum(pv, corr = FALSE, na = FALSE, 
                         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                         symbols = c("***", "**", "*", ".", " "))
        Cf <- cbind(Cf, format(Signif))
      }
    }
    else signif.stars <- FALSE
  }
  else signif.stars <- FALSE
  # Cf2 <- matrix("", nrow(Cf) + 1, ncol(Cf), dimnames = list(c("", rownames(Cf)), colnames(Cf)))
  # Cf2[1, 3:ncol(Cf)] <- Cf[1, 3:ncol(Cf)]
  # Cf2[2:nrow(Cf2), 1:2] <- Cf[, 1:2]
  # Cf <- Cf2
  Cf[nrow(Cf), 1:2] <- ""
  print.default(Cf, quote = quote, right = right, na.print = na.print, 
                ...)
  if (signif.stars && signif.legend) {
    if ((w <- getOption("width")) < nchar(sleg <- attr(Signif, 
                                                       "legend"))) 
      sleg <- strwrap(sleg, width = w - 2, prefix = "  ")
    cat("---\nSignif. codes:  ", sleg, sep = "", fill = w + 
          4 + max(nchar(sleg, "bytes") - nchar(sleg)))
  }
  invisible(x)
}

aic.islasso <- function(object, method = c("aic", "bic"), interval, y, X, intercept = FALSE, family = gaussian(), alpha = 1, offset, weights, unpenalized, 
                          control = is.control()){
  
  if(missing(object)){
    if(missing(interval)) stop("please specify an interval to search for a minimum")
    if(missing(y)) stop("model response is missing")
    if(missing(X)) stop("model matrix is missing")
    n <- NROW(X)
    p <- NCOL(X)
    if(alpha < 0 | alpha > 1) stop("alpha parameter must be in [0,1]")
    if(missing(offset)) offset <- rep(0, n)
    if(missing(weights)) weights <- rep(1, n)
    if(missing(unpenalized)) unpenalized <- NULL
  }else{
    if(missing(interval)) interval <- object$internal$lmbd.interval
    if(is.null(interval)) stop("please specify an interval to search for a minimum")
    n <- object$internal$n
    p <- object$internal$p
    X <- model.matrix(object)
    y <- model.response(object$model)
    intercept <- object$internal$intercept
    nms <- object$internal$nms
    unpenalized <- object$internal$unpenalized
    alpha <- object$alpha
    if(intercept){
      X <- X[,-1]
      nms <- nms[-1]
      unpenalized <- unpenalized[-1]
    }
    unpenalized <- nms[unpenalized]
    family <- object$family
    offset <- object$internal$offset
    weights <- object$internal$weights
    control <- object$control
    if(object$internal$estc) control$pai <- -1
  }
  control$trace <- 0
  
  method <- match.arg(method)
  k <- switch(method,
              "aic"=2,
              "bic"=log(nrow(X)))
  
  fun <- function(lambda, X, y, alpha, family, intercept, weights, offset, unpenalized, control, k, n, p){
    obj <- islasso.fit(X = X, y = y, family = family, lambda = lambda, alpha = alpha, 
                       intercept = intercept, weights = weights, offset = offset, unpenalized = unpenalized, 
                       control = control)
    if((n > p) & (family$family %in% c("gaussian", "binomial", "poisson"))) 
      -2*(obj$rank - obj$aic/2) + k*obj$rank 
    else  
      obj$deviance + k*obj$rank
  }
  
  lambda.min <- optimize(fun, interval=interval, X=X, y=y, alpha=alpha, family=family, intercept=intercept, 
                         weights=weights, offset=offset, unpenalized=unpenalized, control=control, k=k,
                         n=n, p=p)$minimum
  
  return(lambda.min)
}


