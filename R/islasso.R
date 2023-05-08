islasso <- function(formula, family=gaussian, lambda, alpha=1, data, weights, subset, offset,
                    unpenalized, contrasts = NULL, control = is.control()){
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
  }
  else{
    intercept <- TRUE
    X <- X[, -1, drop = FALSE]
  }
  attributes(X)$dataClasses <- temp
  setting <- control
  if(missing(unpenalized)) unpenalized <- NULL
  if(!missing(lambda)){
    if(is.character(lambda)) stop("\n'lambda' can not be a character\n")
    if(is.factor(lambda)) stop("\n'lambda' can not be a factor\n")
    if(lambda < 0) stop("\n'lambda' is negative\n")
  }
  
  fit <- islasso.fit(X=X, y=y, family=family, lambda=lambda, alpha=alpha, intercept=intercept, 
                     weights=weights, offset=offset, unpenalized=unpenalized, control=control)
  fit$model <- mf
  fit$call <- this.call
  fit$formula <- formula
  fit$terms <- mt
  fit$data <- data
  fit$contrasts <- contrasts
  fit$xlevels <- .getXlevels(mt, mf)
  class(fit) <- "islasso"
  
  return(fit)
}

islasso.fit <- function(X, y, family=gaussian, lambda, alpha=1, intercept=FALSE, weights=NULL, 
                        offset=NULL, unpenalized=NULL, control=is.control()){
  this.call <- match.call()
  
  # call the general input checking function
  prep <- .checkinput(X, y, family, alpha, intercept, weights, offset, unpenalized, control)
  
  # call the starting point function
  start <- .startpoint(prep$X, prep$y, lambda, alpha, prep$weights, prep$offset, prep$mustart, prep$family, intercept, prep$setting)
  
  Lambda <- rep(start$lambda, prep$nvars)
  Lambda[prep$unpenalized] <- 0
  
  fam <-  c("binomial","poisson","Gamma"); fam2 <-  c("quasibinomial","quasipoisson","Gamma")
  if((a <- prep$tempFamily %in% fam) | (prep$tempFamily %in% fam2)){
    fam <- ifelse(a, pmatch(prep$tempFamily, fam), pmatch(prep$tempFamily, fam2))
    link <- switch(fam,
                   "1"=pmatch(prep$tempLink, c("logit","probit")),
                   "2"=pmatch(prep$tempLink, c("log")),
                   "3"=pmatch(prep$tempLink, c("inverse","log","identity")))
  }
  else{
    fam <- 0
    link <- 0
  }
  
  out <- .islasso(prep, start, Lambda, fam, link)
  
  return(out)
}

.checkinput <- function(X, y, family, alpha, intercept, weights, offset, unpenalized, control) {
  X <- as.matrix(X); X <- if(intercept) cbind(1, X) else X
  nobs <- nrow(X); nvars <- ncol(X)
  nms <- colnames(X)
  if(is.null(nms) & ncol(X) != 0) nms <- paste0("X", 1:nvars)
  if(intercept) nms[1] <- "(Intercept)"
  colnames(X) <- nms
  
  # check argument
  if(alpha > 1 | alpha < 0) stop("alpha must be in [0, 1] (0 for ridge penalty, 1 for lasso penalty)")
  if(is.null(unpenalized)){
    unpenalized <- rep(FALSE, (nvars-1*intercept))
    if(intercept) unpenalized <- c(TRUE, unpenalized)
  }
  else{
    if(!is.vector(unpenalized)) stop("''unpenalized is not a vector")
    if(is.list(unpenalized)) stop("'unpenalized' can not be a list")
    if(is.factor(unpenalized)) stop("'unpenalized' can not be a factor")
    if(is.character(unpenalized)){
      temp_nms <- if(intercept) nms[-1] else nms
      unpenalized_id <- pmatch(unpenalized, temp_nms)
      if(any(is.na(unpenalized_id))) stop(gettextf("the following names are not in colnames(X): %s", paste(unpenalized[is.na(unpenalized_id)], collapse = ", ")))
      unpenalized <- sort(unpenalized_id)
    }else unpenalized <- sort(unpenalized)
    if(any(abs(unpenalized - round(unpenalized)) > .Machine$double.eps^0.5)) stop("some element of 'unpenalized' is not an integers")
    if(any(unpenalized < 0)) stop("some element of 'unpenalized' is smaller than zero")
    if(any(unpenalized > (nvars-1*intercept))) stop("some element of 'unpenalized' is greater than the number of columns of the matrix 'X'")
    temp <- rep(FALSE, (nvars-1*intercept))
    temp[unpenalized] <- TRUE
    if(intercept) temp <- c(TRUE, temp)
    unpenalized <- temp
  }
  # pen <- X[, !unpenalized, drop=FALSE]
  # unp <- X[, unpenalized, drop=FALSE]
  # X <- cbind(unp, pen)
  # unpenalized <- sort(unpenalized, decreasing = TRUE)
  
  if(is.null(offset)) offset <- rep(0, nobs)
  if(is.null(weights)){
    weights <- as.double(rep(1, nobs))
  }
  else{
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
  setting$estpai <- FALSE
  if(setting$c[1] > 1) stop("'mixing parameter' should be fixed in (0,1) or estimated using a negative value")
  if(setting$c[1] < 0){
    setting$estpai <- TRUE
    setting$c <- .5
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
  
  ## initialize from family function. Makes y a vector in case of binomial, and possibly changes weights
  ## Expects nobs to be defined, and creates n and mustart (neither used here)
  ## Some cases expect to see things, so we set it up just to make it work
  etastart <- 0; mustart <- NULL; start <- NULL
  if(tempFamily == "gaussian"){
    n <- rep.int(1, nobs)
    mustart <- y
  }
  if(tempFamily %in% c("binomial", "quasibinomial")){
    if (NCOL(y) == 1) {
      if (is.factor(y)) y <- y != levels(y)[1L]
      n <- rep.int(1, nobs)
      y[weights == 0] <- 0
      if (any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")
      mustart <- (weights * y + 0.5)/(weights + 1)
      m <- weights * y
      if ("binomial" == "binomial" && any(abs(m - round(m)) > 0.001)) 
        warning(gettextf("non-integer #successes in a %s glm!", "binomial"), domain = NA)
    }
    else if (NCOL(y) == 2) {
      if ("binomial" == "binomial" && any(abs(y - round(y)) > 0.001)) 
        warning(gettextf("non-integer counts in a %s glm!", "binomial"), domain = NA)
      n <- (y1 <- y[, 1L]) + y[, 2L]
      y <- y1 / n
      if (any(n0 <- n == 0)) y[n0] <- 0
      weights <- weights * n
      mustart <- (n * y + 0.5)/(n + 1)
    }
    else stop(gettextf("for the '%s' family, y must be a vector of 0 and 1's\nor a 2 column matrix where col 1 is no. successes and col 2 is no. failures", "binomial"), domain = NA)
  }
  if(tempFamily %in% c("poisson", "quasipoisson")){
    if (any(y < 0)) stop("negative values not allowed for the 'Poisson' family")
    n <- rep.int(1, nobs)
    mustart <- y + 0.1
  }
  if(tempFamily == "Gamma"){
    if (any(y <= 0)) stop("non-positive values not allowed for the 'Gamma' family")
    n <- rep.int(1, nobs)
    mustart <- y
  }
  
  # eval(family$initialize)
  ##
  ## Just in case this was not done in initialize()
  y <- drop(y)  # we don't like matrix responses
  
  return(list(
    y = y,
    X = X,
    intercept = intercept,
    alpha = alpha,
    mustart = mustart,
    weights = weights,
    offset = offset,
    unpenalized = unpenalized,
    nobs = nobs,
    nvars = nvars,
    n = n,
    nms = nms,
    family = family,
    tempFamily = tempFamily,
    tempLink = tempLink,
    setting = setting
  ))
}

.startpoint <- function(X, y, lambda, alpha, weights, offset, mustart, family, intercept, setting) {
  nobs <- nrow(X); nvars <- ncol(X)
  
  if(family$family %in% c("gaussian", "binomial", "quasibinomial", "poisson", "quasipoisson")){
      tempFamily <-  family$family
      if(tempFamily == "quasibinomial") tempFamily <- "binomial"
      if(tempFamily == "quasipoisson") tempFamily <- "poisson"
  }
  else tempFamily <- family
  
  if(is.null(setting$b0)){
    if((nvars - intercept) < 2){
      if(missing(lambda)) stop("Insert a positive value for lambda")
      # if(missing(lambda)){lambda <- sqrt(log(nvars)/nobs)}
      est <- rep(0.1, nvars)
    }
    else{
      x <- as.matrix(X)
      y2 <- y
      weights2 <- weights
      if(family$family %in% c("binomial", "quasibinomial")) {
          y2 <- weights * cbind(1-y, y)
          weights2 <- rep.int(1, nobs)
      }
      if(intercept) x <- x[,-1]
      if(missing(lambda)){
        obj <- suppressWarnings(cv.glmnet(x=x, y=y2, family=tempFamily, nfolds=setting$nfolds,
                         standardize=setting$stand, intercept=intercept, offset=offset, weights=weights2, alpha=alpha))
        lambda <- obj$lambda.min * nobs
        est <- as.vector(coef(obj, s="lambda.min"))
      }
      else{
        obj <- suppressWarnings(glmnet(x=x, y=y2, family=tempFamily, alpha=alpha, weights=weights2,
                      standardize=setting$stand, intercept=intercept, offset=offset))
        est <- as.vector(coef(obj, s=lambda / nobs))
      }
      if(setting$stand) {
        x_mean <- colMeans(x)
        x_centered <- sweep(x, 2, x_mean, "-")
        x_sd <- apply(x_centered, 2, function(x) sqrt(sum(x^2) / nobs))
        est[1] <- est[1] + sum(x_mean / x_sd * est[1])
        est[-1] <- x_sd * est[-1]
      }
      if(!intercept) est <- est[-1]
    }
  }
  else{
    if(missing(lambda)) stop("Insert a positive value for lambda")
    est <- setting$b0
  }
  interval <- if(exists("obj")) range(rev(obj$lambda))*nobs else NULL
  
  covar <- if(is.null(setting$V0)) diag(.01, nvars) else setting$V0
  se <- sqrt(diag(covar))
  eta <- family$linkfun(mustart) + offset
  mu <- family$linkinv(eta)
  residuals <- (y - mu) / family$mu.eta(eta)
  
  return(list(
    fam = family$family,
    lambda = lambda,
    interval = interval,
    beta = est,
    covar = covar,
    se = se,
    eta = eta,
    mu = mu,
    residuals = residuals
  ))
}

.islasso <- function(prep, start, Lambda, fam, link) {
  setting <- prep$setting
  X <- prep$X
  storage.mode(X) <- "double"
  y <- prep$y
  storage.mode(y) <- "double"
  intercept <- prep$intercept
  storage.mode(intercept) <- "integer"
  nobs <- prep$nobs
  storage.mode(nobs) <- "integer"
  nvars <- prep$nvars
  storage.mode(nvars) <- "integer"
  beta <- start$beta
  if(any(beta == 0)) beta[beta == 0] <- 1E-5
  storage.mode(beta) <- "double"
  se <- start$se
  storage.mode(se) <- "double"
  covar <- start$covar
  storage.mode(covar) <- "double"
  storage.mode(Lambda) <- "double"
  alpha <- prep$alpha
  storage.mode(alpha) <- "double"
  offset <- prep$offset
  storage.mode(offset) <- "double"
  weights <- prep$weights
  storage.mode(weights) <- "double"
  eta <- start$eta
  storage.mode(eta) <- "double"
  mu <- start$mu
  storage.mode(mu) <- "double"
  residuals <- start$residuals
  storage.mode(residuals) <- "double"
  
  storage.mode(setting$c) <- "double"
  storage.mode(setting$estpai) <- "integer"
  storage.mode(setting$trace) <- "integer"
  storage.mode(setting$adaptive) <- "integer"
  storage.mode(setting$stand) <- "integer"
  h <- as.double(1)
  storage.mode(setting$itmax) <- "integer"
  storage.mode(setting$tol) <- "double"
  storage.mode(setting$sigma2) <- "double"
  
  if(prep$tempFamily == "gaussian") {
    # if(setting$sandwich) {
      fit <- .Fortran(C_islasso, X = X, y = y, n = nobs, p = nvars, ntheta = beta, se = se, cov = covar, lambda = Lambda,
                      alpha = alpha, pi = setting$c, estpi = setting$estpai, h = h, itmax = setting$itmax, itmaxse = setting$itmax, tol = setting$tol, 
                      phi = setting$sigma2, trace = setting$trace, adaptive = setting$adaptive, offset = offset, conv = integer(1), 
                      stand = setting$stand, intercept = intercept, eta = eta, mu = mu, res = residuals, dev = double(1), weights = weights, 
                      hi = double(nvars), edf = double(1), grad = double(nvars))
    # }
    # else {
    #   fit <- .Fortran(C_islasso3, X = X, y = y, n = nobs, p = nvars, ntheta = beta, se = se, cov = covar, lambda = Lambda,
    #                   alpha = alpha, pi = setting$c, estpi = setting$estpai, h = h, itmax = setting$itmax, itmaxse = setting$itmax, tol = setting$tol, 
    #                   phi = setting$sigma2, trace = setting$trace, adaptive = setting$adaptive, offset = offset, conv = integer(1), 
    #                   stand = setting$stand, intercept = intercept, eta = eta, mu = mu, res = residuals, dev = double(1), weights = weights, 
    #                   hi = double(nvars), edf = double(1), grad = double(nvars))
    # }
  }
  else {
    # if(setting$sandwich) {
      fit <- .Fortran(C_islasso_glm, X = X, y = y, n = nobs, p = nvars, ntheta = beta, se = se, cov = covar, lambda = Lambda, 
                      alpha = alpha, pi = setting$c, estpi = setting$estpai, h = h, itmax = setting$itmax, itmaxse = setting$itmax, tol = setting$tol, 
                      phi = setting$sigma2, trace = setting$trace, adaptive = setting$adaptive, offset = offset, conv = integer(1), 
                      stand = setting$stand, intercept = intercept, eta = eta, mu = mu, dev = double(1), weights = weights, 
                      hi = double(nvars), edf = double(1), fam = as.integer(fam), link = as.integer(link), grad = double(nvars))
    # }
    # else {
    #   fit <- .Fortran(C_islasso_glm2, X = X, y = y, n = nobs, p = nvars, ntheta = beta, se = se, cov = covar, lambda = Lambda, 
    #                   alpha = alpha, pi = setting$c, estpi = setting$estpai, h = h, itmax = setting$itmax, itmaxse = setting$itmax, tol = setting$tol, 
    #                   phi = setting$sigma2, trace = setting$trace, adaptive = setting$adaptive, offset = offset, conv = integer(1), 
    #                   stand = setting$stand, intercept = intercept, eta = eta, mu = mu, dev = double(1), weights = weights, 
    #                   hi = double(nvars), edf = double(1), fam = as.integer(fam), link = as.integer(link), grad = double(nvars))
    # }
  }
  if(fit$conv == -1) stop("Infinite values attained, try to change lambda value!!")
  if(fit$conv == 1) warning("Maximum number of iterations attained!!")
  # if(fit$conv == 2) stop("Exit from lm algorithm after an inversion problem or a step halving too small, try to change lambda value!!")
  if(fit$conv == 2) stop("Safe exit from ISLASSO algorithm after an inversion problem, try to change lamda value!!")
  
  setting$c <- fit$pi
  
  # extract family functions
  family <- prep$family
  variance <- family$variance
  linkinv <- family$linkinv
  dev.resids <- family$dev.resids
  aic <- family$aic
  mu.eta <- family$mu.eta
  linkfun <- family$linkfun
  
  # extract from fit
  rank <- fit$edf
  iter <- fit$itmax
  conv <- fit$conv
  
  beta <- fit$ntheta
  eta <- drop(X %*% beta + offset)
  mu <- linkinv(eta)
  dev <- sum(dev.resids(y, mu, weights))
  #s2 <- if(prep$tempFamily %in% c("binomial", "poisson")) 1.0 else dev / (nobs - rank)
  
  s2 <- fit$phi
  covar <- s2 * fit$cov
  #se <- sqrt(diag(covar))
  se <- fit$se
  # dev <- fit$dev
  
  # aic.model <- aic(y, nobs, mu, weights, dev) + 2 * rank
  aic.model <- if((nobs > nvars) & (prep$tempFamily %in% c("gaussian", "binomial", "poisson", "Gamma")))
    aic(y, prep$n, mu, weights, dev) + 2 * rank
  else
    dev + 2 * rank
  
  varmu <- variance(mu)
  mu.eta.val <- mu.eta(eta)
  w <- (weights * mu.eta.val^2) / varmu
  residuals <- (y - mu) / mu.eta.val
  z <- (eta - offset) + residuals
  
  bsc <- (beta / se)
  
  P <- alpha * (setting$c * (2 * pnorm(bsc, 0, 1) - 1) + (1 - setting$c) * (2 * pnorm(bsc, 0, 1E-5) - 1)) + 
    (1 - alpha) * beta
  if(any(prep$unpenalized)) P[prep$unpenalized] <- 0
  
  # tmp <- lm.wfit(rbind(X, sqrt(Lambda) * diag(sqrt(P / beta))), 
  #                    c(z, rep(0, nvars)), 
  #                    c(w, rep(0, nvars)))
  # tmp2 <- .Call(stats:::Cdqrls,
  #              rbind(X * sqrt(w), sqrt(Lambda) * diag(sqrt(P / beta))),
  #              c(z * sqrt(w), rep(0, nvars)), setting$tol, check = FALSE)
  fit$qr <- .lm.fit(x = rbind(X * sqrt(w), sqrt(Lambda) * diag(sqrt(P / beta))), 
                 y = c(z * sqrt(w), rep(0, nvars)))
  # fit$qr <- structure(tmp[c("qr", "rank", "qraux", "pivot", "tol")], class = "qr")
  nr <- min(nobs, nvars)
  if (nr < nvars) {
    Rmat <- diag(nvars)
    Rmat[1L:nr, 1L:nvars] <- fit$qr$qr[1L:nr, 1L:nvars]
  }
  else Rmat <- fit$qr$qr[1L:nvars, 1L:nvars]
  Rmat <- as.matrix(Rmat)
  Rmat[row(Rmat) > col(Rmat)] <- 0
  fit$effects <- fit$qr$effects[1:nobs]
  
  # fn0 <- function(grad, b, s, c, lambda, alpha, unpenalized) {
  #   bsc <- (b / s)
  #   # grad <- crossprod(x, y - drop(x %*% b))
  #   r <- alpha * (c*(2 * pnorm(bsc, 0, 1) - 1) + (1 - c)*(2 * pnorm(bsc, 0, .00001) - 1)) + (1 - alpha) * b
  #   if(any(unpenalized)) r[unpenalized] <- 0
  #   return(drop(- grad + lambda * r))
  # }
  # gradient <- fn0(grad = crossprod(X, w * residuals), beta, se, setting$c, start$lambda, alpha, prep$unpenalized)
  gradient <- drop(- crossprod(X, w * residuals) + Lambda * P)
  
  XtW <- t(w * X)
  XtX <- XtW %*% X
  H <- .Fortran(C_hessian, beta, se, Lambda, XtX, setting$c, nvars, hess = XtX, alpha)$hess
  invH <- .Fortran(C_inv, nvars, H, invA = H, integer(1))$invA
  
  # null model
  wtdmu <- if(intercept) sum(weights * y)/sum(weights) else linkinv(offset)
  nulldev <- sum(dev.resids(y, wtdmu, weights))
  
  n.ok <- nobs - sum(weights == 0)
  nulldf <- n.ok - as.integer(intercept)
  resdf <- n.ok - rank
  
  # # compute unbias estimate (beta code)
  # se.unbias <- beta.unbias <- NULL
  # if(setting$debias){
  #   invH2 <- ginv2(XtX) #.Fortran(C_inv2, as.integer(nvars), XtX, invA=XtX, integer(0))$invA
  #   # invH2 <- .Fortran(C_inv3, as.integer(nvars), XtX, invA=XtX, integer(1))$invA
  #   gradient2 <- XtW %*% residuals
  #   delta <- invH2 %*% gradient2
  #   beta.unbias <- drop(beta + delta)
  #   se.unbias <-  sqrt(diag(s2*(invH2)))
  #   names(beta.unbias) <- names(se.unbias) <- prep$nms
  # }
  
  names(gradient) <- names(se) <- names(beta) <- colnames(XtX) <- 
    rownames(XtX) <- colnames(covar) <- rownames(covar) <- prep$nms
  
  # internal argument
  internal <- list(n = nobs, p = nvars, lambda.seq = fit$lambda, 
                   XtW = XtW, XtX = XtX, invH = invH, vcov = covar, 
                   gradient = gradient, hessian = H, hi = fit$hi, 
                   intercept = intercept, unpenalized = prep$unpenalized, 
                   fam = fam, link = link, nms = prep$nms, 
                   estc = setting$estpai, lmbd.interval = start$interval)
  
  out <- list(coefficients = beta, se = se, dispersion = s2, residuals = residuals, fitted.values = mu, effects = fit$effect, 
              R = Rmat, rank = rank, qr = structure(fit$qr[c("qr", "qraux", "pivot", "tol", "rank")], class = "qr"), family = family, 
              linear.predictors = eta, deviance = dev, 
              aic = aic.model, null.deviance = nulldev, iter = fit$itmax, weights = w, prior.weights = weights, 
              df.residual = resdf, df.null = nulldf, y = y, converged = fit$conv, model = NULL, call = NULL, 
              formula = NULL, terms = NULL, data = NULL, offset = offset, contrasts = NULL, 
              control = setting, internal = internal, contrasts = NULL, xlevels = NULL, lambda = start$lambda, alpha = alpha)
  
  out
}

is.control <- function(sigma2 = -1, tol = 1E-05, itmax = 1E+3, stand = TRUE,
                       trace = 0, nfolds = 5, seed = NULL, adaptive = FALSE, g = .5,
                       b0 = NULL, V0 = NULL, c = .5) {
  
  list(sigma2 = sigma2, tol = tol, itmax = itmax, trace = trace, stand = stand, nfolds = nfolds, seed = seed,#, debias=debias
       adaptive = adaptive, g = g, b0 = b0, V0 = V0, c = c)
}

aic.islasso <- function(object, method = c("AIC", "BIC", "AICc", "GCV", "GIC"), interval, g = 0,
                        y, X, intercept = FALSE, family = gaussian(), alpha = 1, offset, weights, 
                        unpenalized, control = is.control(), trace = TRUE){
  
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
  }
  else{
    if(missing(interval)) interval <- object$internal$lmbd.interval
    if(is.null(interval)) stop("please specify an interval to search for a minimum")
    n <- object$internal$n
    p <- object$internal$p
    X <- model.matrix(object)
    y <- object$y#model.response(object$model)
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
    offset <- object$offset
    weights <- object$weights
    control <- object$control
    if(control$estpai) control$c <- -1
    if(length(control$c) != 1L) control$c <- max(control$c)
  }
  control$trace <- 0
  
  method <- match.arg(method)
  k <- switch(method, 
              "AIC" = "ll + 2 * df", 
              "BIC" = "ll + (log(n) + 2 * g * log(p)) * df", 
              "AICc" = "ll + 2 *  n * df * (df + 1) / (n - df - 1) + 2 * df", 
              "GCV" = "ll / ((1 - df / n)^2)",
              "GIC" = "ll + log(log(n)) * log(p) * df")
  if(method == "BIC" & g != 0) method <- "eBIC"
  if(g < 0 | g > 1) stop("gamma parameter have to be set in (0,1)")
  
  if(trace) cat(paste0("\nOptimization through ", method, "\n\n"))
  fun <- function(lambda, X, y, alpha, family, intercept, weights, offset, unpenalized, control, k, n, p, trace){
    obj <- try(islasso.fit(X = X, y = y, family = family, lambda = lambda, alpha = alpha, intercept = intercept,
                       weights = weights, offset = offset, unpenalized = unpenalized, control = control), silent = TRUE)
    #if(trace) cat(".")
    if(inherits(obj, "try-error")) {
        temp <- Inf #.Machine$double.xmax
    }
    else{
      # temp <- if(k == 2){ obj$aic }else{ obj$aic - 2 * obj$rank + k * obj$rank }
      ll <- obj$aic - 2 * obj$rank
      df <- obj$rank
      temp <- eval(parse(text = k))
        # temp <- if((n > p) & (family$family %in% c("gaussian", "binomial", "poisson", "Gamma")))
        #   if(k == 2){ obj$aic }else{ obj$aic - 2 * obj$rank + k * obj$rank }
        # else
        #   obj$deviance + k * obj$rank
    }
    #if(trace) print(c(lambda, temp))
    cat("lambda = ", formatC(lambda, digits = 4, width = 8, format = "f"), method, "= ", formatC(temp, digits = 5, width = 10, format = "f"), "\n")
    temp
  }
  
  lambda.min <- optimize(fun, interval=interval, X=X, y=y, alpha=alpha, family=family, intercept=intercept, 
                         weights=weights, offset=offset, unpenalized=unpenalized, control=control, k=k,
                         n=n, p=p, trace=trace)$minimum
  
  return(lambda.min)
}
