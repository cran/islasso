fmt.perc <- function (probs, digits){
  paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), "%")
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

grad.islasso <- function(object, ...) {
  x <- model.matrix(object)
  y <- object$y
  b <- coef(object)
  s <- object$se
  offset <- object$offset
  c <- object$control$c
  lambda <- object$lambda
  alpha <- object$alpha
  unpenalized <- object$internal$unpenalized
  family <- family(object)
  weights <- weights(object)
  
  eta <- x %*% b + offset
  mu <- family$linkinv(eta)
  v <- family$variance(mu)
  m.e <- family$mu.eta(eta)
  r <- y - mu
  rv <- r / v * m.e * weights
  grad <- drop(t(rv) %*% x)
  
  bsc <- (b / s)
  r <- alpha * (c * (2 * pnorm(bsc, 0, 1) - 1) + (1 - c) * (2 * pnorm(bsc, 0, 1E-5) - 1)) + (1 - alpha) * b
  if(any(unpenalized)) r[unpenalized] <- 0
  
  return(- grad + lambda * r)
}

hess.islasso <- function(object, ...) {
  H <- object$internal$XtX
  b <- coef(object)
  s <- object$se
  c <- object$control$c
  lambda <- object$lambda
  alpha <- object$alpha
  
  bsc <- (b / s)
  diag(H) <- diag(H) + 
    lambda * (2 * alpha * (c * dnorm(bsc, 0, 1) + (1 - c) * dnorm(bsc, 0, 1E-5)) / s + (1 - alpha))
  return(H)
}

is.influence <- function (model, do.coef = TRUE) {
  infl <- function(mqr, do.coef, e, tol) {
    qr <- mqr$qr
    qraux <- mqr$qraux
    n <- nrow(qr)
    k <- mqr$rank
    if(length(e) != n) e <- c(e, rep(0, n - length(e)))
    tmp <- .Fortran(C_lminfl, qr, n, n, k, as.integer(1), qraux, e, hat = double(n), 
                    sigma = double(n), tol)
    # tmp <- .Fortran(C_lminfl, qr, n, n, k, as.integer(do.coef), qraux, e, hat = double(n), 
    #                 coefficients = matrix(0.0, n, k), sigma = double(n), tol)
    for (i in seq_len(n)) {
      if (tmp$hat[i] > 1. - tol) tmp$hat[i]  <-  1.
    }
    tmp[c("hat", "sigma")]
  }
  
  wt.res <- weighted.residuals(model)
  e <- na.omit(wt.res)
  n <- length(wt.res)
  is.mlm <- is.matrix(e)
  if (model$rank == 0) {
    n <- length(wt.res)
    sigma <- sqrt(deviance(model)/df.residual(model))
    res <- list(hat = rep(0, n), coefficients = matrix(0, n, 0), sigma = rep(sigma, n))
  }
  else {
    e[abs(e) < 100 * .Machine$double.eps * median(abs(e))] <- 0
    mqr <- model$qr
    do.coef <- as.logical(do.coef)
    tol <- 10 * .Machine$double.eps
    # res2 <- .Call(stats:::C_influence, mqr, e, tol)
    res <- infl(mqr, FALSE, e, tol)
    
    res$hat <- res$hat[1:n]
    res$sigma <- res$sigma[1:n]
    if (do.coef) {
      ok <- seq_len(mqr$rank)
      Q <- qr.Q(mqr)[, ok, drop = FALSE]
      R <- qr.R(mqr)[ok, ok, drop = FALSE]
      hat <- res$hat
      invRQtt <- t(backsolve(R, t(Q)))
      k <- NCOL(Q)
      q <- NCOL(e)
      if (is.mlm) {
        cf <- array(0, c(n, k, q))
        for (j in seq_len(q)) cf[, , j] <- invRQtt[1:n, , drop = FALSE] * ifelse(hat == 1, 0, e[, j]/(1 - hat))
      }
      else cf <- invRQtt[1:n, , drop = FALSE] * ifelse(hat == 1, 0, e/(1 - hat))
      res$coefficients <- cf
    }
    drop1d <- function(a) {
      d <- dim(a)
      if (length(d) == 3L && d[[3L]] == 1L) dim(a) <- d[-3L]
      a
    }
    if (is.null(model$na.action)) {
      if (!is.mlm) {
        res$sigma <- drop(res$sigma)
        if (do.coef) res$coefficients <- drop1d(res$coefficients)
      }
    }
    else {
      hat <- naresid(model$na.action, res$hat)
      hat[is.na(hat)] <- 0
      res$hat <- hat
      if (do.coef) {
        coefficients <- naresid(model$na.action, res$coefficients)
        coefficients[is.na(coefficients)] <- 0
        res$coefficients <- if (is.mlm) 
          coefficients
        else drop1d(coefficients)
      }
      sigma <- naresid(model$na.action, res$sigma)
      sigma[is.na(sigma)] <- sqrt(deviance(model)/df.residual(model))
      res$sigma <- if (is.mlm) 
        sigma
      else drop(sigma)
    }
  }
  res$wt.res <- naresid(model$na.action, e)
  res$hat[res$hat > 1 - 10 * .Machine$double.eps] <- 1
  names(res$hat) <- names(res$sigma) <- names(res$wt.res)
  if (do.coef) {
    cf <- coef(model)
    if (is.mlm) {
      dnr <- dimnames(res$wt.res)
      dimnames(res$coefficients) <- list(dnr[[1L]], rownames(cf)[!apply(cf, 1L, anyNA)], dnr[[2L]])
    }
    else dimnames(res$coefficients) <- list(names(res$wt.res), names(cf)[!is.na(cf)])
  }
  res[c("hat", "coefficients", "sigma", "wt.res")]
}

islasso.diag <- function (glmfit) {
  w <- if (is.null(glmfit$prior.weights)) 
    rep(1, length(glmfit$residuals))
  else glmfit$prior.weights
  sd <- sqrt(glmfit$dispersion)
  dev <- residuals(glmfit, type = "deviance") / sd
  pear <- residuals(glmfit, type = "pearson") / sd
  h <- rep(0, length(w))
  h[w != 0] <- is.influence(glmfit)$hat
  p <- glmfit$rank
  rp <- pear/sqrt(1 - h)
  rd <- dev/sqrt(1 - h)
  cook <- (h * rp^2)/((1 - h) * p)
  res <- sign(dev) * sqrt(dev^2 + h * rp^2)
  list(res = res, rd = rd, rp = rp, cook = cook, h = h, sd = sd)
}

islasso.diag.plots <- function (glmfit, glmdiag = islasso.diag(glmfit), subset = NULL, 
                                iden = FALSE, labels = NULL, ret = FALSE) {
  if (is.null(glmdiag)) glmdiag <- islasso.diag(glmfit)
  if (is.null(subset)) 
    subset <- seq_along(glmdiag$h)
  else if (is.logical(subset)) 
    subset <- seq_along(subset)[subset]
  else if (is.numeric(subset) && all(subset < 0)) 
    subset <- (1L:(length(subset) + length(glmdiag$h)))[subset]
  else if (is.character(subset)) {
    if (is.null(labels)) 
      labels <- subset
    subset <- seq_along(subset)
  }
  
  par(mfrow = c(2, 2))
  x1 <- predict(glmfit)
  plot(x1, glmdiag$res, xlab = "Linear predictor", ylab = "Residuals")
  pars <- vector(4L, mode = "list")
  pars[[1L]] <- par("usr")
  y2 <- glmdiag$rd
  x2 <- qnorm(ppoints(length(y2)))[rank(y2)]
  plot(x2, y2, ylab = "Quantiles of standard normal", xlab = "Ordered deviance residuals")
  abline(0, 1, lty = 2)
  pars[[2L]] <- par("usr")
  hh <- glmdiag$h/(1 - glmdiag$h)
  plot(hh, glmdiag$cook, xlab = "h/(1-h)", ylab = "Cook statistic")
  rx <- range(hh)
  ry <- range(glmdiag$cook)
  rank.fit <- glmfit$rank
  nobs <- rank.fit + glmfit$df.residual
  cooky <- 8/(nobs - 2 * rank.fit)
  hy <- (2 * rank.fit)/(nobs - 2 * rank.fit)
  if ((cooky >= ry[1L]) && (cooky <= ry[2L])) 
    abline(h = cooky, lty = 2)
  if ((hy >= rx[1L]) && (hy <= rx[2L])) 
    abline(v = hy, lty = 2)
  pars[[3L]] <- par("usr")
  plot(subset, glmdiag$cook, xlab = "Case", ylab = "Cook statistic")
  if ((cooky >= ry[1L]) && (cooky <= ry[2L])) 
    abline(h = cooky, lty = 2)
  xx <- list(x1, x2, hh, subset)
  yy <- list(glmdiag$res, y2, glmdiag$cook, glmdiag$cook)
  pars[[4L]] <- par("usr")
  if (is.null(labels)) 
    labels <- names(x1)
  while (iden) {
    cat("****************************************************\n")
    cat("Please Input a screen number (1,2,3 or 4)\n")
    cat("0 will terminate the function \n")
    num <- as.numeric(readline())
    if ((length(num) > 0L) && ((num == 1) || (num == 2) || 
                               (num == 3) || (num == 4))) {
      cat(paste("Interactive Identification for screen", 
                num, "\n"))
      cat("left button = Identify, center button = Exit\n")
      nm <- num + 1
      par(mfg = c(trunc(nm/2), 1 + nm%%2, 2, 2))
      par(usr = pars[[num]])
      identify(xx[[num]], yy[[num]], labels)
    }
    else iden <- FALSE
  }
  par(mfrow = c(1, 1))
  if (ret) 
    glmdiag
  else invisible()
}

predislasso <- function(object, newdata, type = c("response", "terms"), 
                        terms = NULL, na.action = na.pass, ...){
  type <- match.arg(type)
  
  tt <- terms(object)
  if (missing(newdata) || is.null(newdata)) {
    mm <- X <- model.matrix(object)
    mmDone <- TRUE
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
    mmDone <- FALSE
  }
  
  n <- object$internal$n
  p <- object$internal$p
  p1 <- seq_len(p)
  piv <- if (p) object$qr$pivot[p1]
  beta <- object$coefficients
  predictor <- drop(X[, piv, drop = FALSE] %*% beta[piv])
  if (!is.null(offset)) predictor <- predictor + offset
  
  if (type == "terms") {
    if (!mmDone) {
      mm <- model.matrix(object)
      mmDone <- TRUE
    }
    aa <- attr(mm, "assign")
    ll <- attr(tt, "term.labels")
    hasintercept <- attr(tt, "intercept") > 0L
    if (hasintercept) 
      ll <- c("(Intercept)", ll)
    aaa <- factor(aa, labels = ll)
    asgn <- split(order(aa), aaa)
    if (hasintercept) {
      asgn$"(Intercept)" <- NULL
      avx <- colMeans(mm)
      termsconst <- sum(avx[piv] * beta[piv])
    }
    nterms <- length(asgn)
    if (nterms > 0) {
      predictor <- matrix(ncol = nterms, nrow = NROW(X))
      dimnames(predictor) <- list(rownames(X), names(asgn))
      if (hasintercept) 
        X <- sweep(X, 2L, avx, check.margin = FALSE)
      unpiv <- rep.int(0L, NCOL(X))
      unpiv[piv] <- p1
      for (i in seq.int(1L, nterms, length.out = nterms)) {
        iipiv <- asgn[[i]]
        ii <- unpiv[iipiv]
        iipiv[ii == 0L] <- 0L
        predictor[, i] <- if (any(iipiv > 0L)) 
          X[, iipiv, drop = FALSE] %*% beta[iipiv]
        else 0
      }
      if (!is.null(terms)) predictor <- predictor[, terms, drop = FALSE]
    }
    else predictor <- ip <- matrix(0, n, 0L)
    attr(predictor, "constant") <- if (hasintercept) 
      termsconst
    else 0
  }
  if (missing(newdata) && !is.null(na.act <- object$na.action)) predictor <- napredict(na.act, predictor)
  attr(predictor, "X") <- X
  return(predictor)
}

'%||%' <- function (L, R) if (is.null(L)) R else L

