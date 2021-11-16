islasso.path  <-  function(formula, family = gaussian(), lambda = NULL, nlambda = 100, lambda.min.ratio = ifelse(nobs < nvars, 1E-3, 1E-05), 
                           alpha = 1, data, weights, subset, offset, contrasts = NULL, unpenalized, control = is.control()) {
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
  x <- if(!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else stop("Model matrix is empty!")
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  temp <- which(attr(mt, "dataClasses")[names(attr(mt, "dataClasses")) %in% attr(mt, "term.labels")] %in% c("factor", "character"))
  temp <- which(attr(x, "assign") %in% temp)
  if(ioff){
    noff <- match(unlist(lapply(off, as.character)), colnames(x))
    if(!all(is.na(noff))) x <- x[, -noff[which(noff!=0)]]
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
    x <- x[, -1, drop = FALSE]
  }
  attributes(x)$dataClasses <- temp
  if(missing(unpenalized)) unpenalized <- NULL
  
  fit <- islasso.path.fit(X = x, y = y, family = family, lambda = lambda, nlambda = nlambda, 
                          lambda.min.ratio = lambda.min.ratio, alpha = alpha, intercept = intercept, 
                          weights = weights, offset = offset, unpenalized = unpenalized, control = control)
  
  fit$call <- this.call
  fit$formula <- formula
  fit$model <- mf
  fit$terms <- mt
  fit$data <- data
  fit$contrasts <- contrasts
  fit$xlevels <- .getXlevels(mt, mf)
  fit$contrasts <- contrasts
  class(fit) <- 'islasso.path'
  return(fit)
}

islasso.path.fit <- function(X, y, family = gaussian(), lambda, nlambda, 
                             lambda.min.ratio, alpha = 1, intercept = FALSE, 
                             weights = NULL, offset = NULL, 
                             unpenalized = NULL, control = is.control()){
  this.call <- match.call()
  
  # call the general input checking function
  prep <- checkinput.islasso.path(X, y, family, lambda, nlambda, lambda.min.ratio, alpha, intercept, weights, offset, unpenalized, control)
  
  # call the starting point function
  start <- startpoint.islasso.path(prep$X, prep$y, prep$lambda, alpha, prep$weights, prep$offset, prep$mustart, prep$family, intercept, prep$setting)
  
  # Lambda[prep$unpenalized] <- 0
  
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
  # browser()
  # out <- islasso.path.fit.lm(prep, start, prep$lambda, fam, link)
  out <- islasso.path.fit.glm(prep, start, prep$lambda, fam, link)
  out$control <- prep$setting
  return(out)
}

checkinput.islasso.path <- function(X, y, family, lambda, nlambda, lambda.min.ratio, alpha, intercept, 
                                    weights, offset, unpenalized, control) {
  get_start <- function (x, y, weights, family, intercept, offset, exclude, alpha) {
    nobs <- nrow(x)
    nvars <- ncol(x)
    is.offset <- ifelse(all(offset == 0), FALSE, TRUE)
    if (intercept) {
      if (is.offset) {
        suppressWarnings(tempfit <- glm(y ~ 1, family = family, weights = weights, offset = offset))
        mu <- tempfit$fitted.values
      }
      else {
        mu <- rep(weighted.mean(y, weights), times = nobs)
      }
    }
    else {
      mu <- family$linkinv(offset)
    }
    nulldev <- sum(family$dev.resids(y, mu, weights))
    ju <- rep(1, nvars)
    ju[exclude] <- 0
    r <- y - mu
    eta <- family$linkfun(mu)
    v <- family$variance(mu)
    m.e <- family$mu.eta(eta)
    rv <- r/v * m.e * weights
    g <- abs(drop(t(rv) %*% x)) * ju
    lambda_max <- max(g)/max(alpha, 0.001)
    list(nulldev = nulldev, mu = mu, lambda_max = lambda_max)
  }
  
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
    if(any(unpenalized <= 0)) stop("some element of 'unpenalized' is smaller than zero")
    if(any(unpenalized > (nvars-1*intercept))) stop("some element of 'unpenalized' is greater than the number of columns of the matrix 'X'")
    temp <- rep(FALSE, (nvars-1*intercept))
    temp[unpenalized] <- TRUE
    if(intercept) temp <- c(TRUE, temp)
    unpenalized <- temp
  }
  
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
  if(setting$itmax < 0) stop("'itmax' should be a non-negative value")
  if(setting$c[1] > 1) stop("'mixing parameter' should be fixed in (0,1) or estimated using a negative value")
  if(setting$c[1] < 0) {
    setting$c <- rep(1, nvars)
    setting$fix.c <- FALSE
  } 
  else {
    if(length(setting$c) != nvars) setting$c <- rep(setting$c[1], nvars)
    setting$fix.c <- TRUE
  }
  if(setting$g < 0 | setting$g > 1) {
    warning("gamma parameter have to be set in (0,1). Default parameter is 0.5")
    setting$g <- .5
  }
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
    mustart <- y
  }
  if(tempFamily %in% c("binomial", "quasibinomial")){
    if (NCOL(y) == 1) {
      if (is.factor(y)) y <- y != levels(y)[1L]
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
    mustart <- y + 0.1
  }
  if(tempFamily == "Gamma"){
    if (any(y <= 0)) stop("non-positive values not allowed for the 'Gamma' family")
    mustart <- y
  }
  
  # eval(family$initialize)
  ##
  ## Just in case this was not done in initialize()
  y <- drop(y)  # we don't like matrix responses
  
  if(is.null(lambda)) {
    lmax <- get_start(X, y, weights, family, intercept, offset, unpenalized, alpha)$lambda_max
    # lmax <- max(abs(crossprod(x, y))) * 1.01 * n
    lmin <- lmax * lambda.min.ratio * 1.01
    lambda <- exp(seq(log(lmin), log(lmax), length.out = nlambda))
  } 
  else {
    nlambda <- length(lambda)
  }
  
  return(list(
    y = y,
    X = X,
    intercept = intercept,
    lambda = lambda,
    nlambda = nlambda,
    alpha = alpha,
    mustart = mustart,
    weights = weights,
    offset = offset,
    unpenalized = unpenalized,
    nobs = nobs,
    nvars = nvars,
    nms = nms,
    family = family,
    tempFamily = tempFamily,
    tempLink = tempLink,
    setting = setting
  ))
}

startpoint.islasso.path <- function(X, y, lambda, alpha, weights, offset, mustart, family, intercept, setting) {
  nobs <- nrow(X); nvars <- ncol(X)
  
  if(family$family %in% c("gaussian", "binomial", "quasibinomial", "poisson", "quasipoisson")){
    tempFamily <-  family$family
    if(tempFamily == "quasibinomial") tempFamily <- "binomial"
    if(tempFamily == "quasipoisson") tempFamily <- "poisson"
  }
  else tempFamily <- family
  
  if((nvars - intercept) < 2){
    est <- rep(0.1, nvars)
  }
  else{
    x <- as.matrix(X)
    y2 <- y
    weights2 <- weights
    if(family$family %in% c("binomial", "quasibinomial")) {
      y2 <- weights * cbind(1 - y, y)
      weights2 <- rep.int(1, nobs)
    }
    if(intercept) x <- x[, -1]
    obj <- suppressWarnings(glmnet(x = x, y = y2, family = tempFamily, alpha = alpha, weights = weights2,
                                   standardize = TRUE, intercept = intercept, offset = offset))
    est <- as.vector(coef(obj, s = min(lambda) / nobs))
    if(!intercept) est <- est[-1]
  }
  interval <- range(lambda)
  
  covar <- diag(.01, nrow = nvars, ncol = nvars)
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

islasso.path.fit.glm <- function(prep, start, lambda, fam, link){
  
  cov.unscaled <- start$covar
  se <- start$se
  b <- start$beta + .01
  
  setting <- prep$setting
  c <- setting$c
  fix.c <- setting$fix.c
  est.c <- !fix.c
  tol <- setting$tol
  maxIter <- setting$itmax
  trace <- setting$trace
  sigma2 <- setting$sigma2
  g <- setting$g
  nlambda <- length(lambda)
  alpha <- prep$alpha
  nvars <- prep$nvars
  nobs <- prep$nobs
  intercept <- prep$intercept
  x <- prep$X
  y <- prep$y
  weights <- prep$weights
  offset <- prep$offset
  unpenalized <- prep$unpenalized
  family <- prep$family
  
  storage.mode(x) <- "double"
  storage.mode(y) <- "double"
  storage.mode(intercept) <- "integer"
  storage.mode(nobs) <- "integer"
  storage.mode(nvars) <- "integer"
  storage.mode(b) <- "double"
  storage.mode(se) <- "double"
  storage.mode(cov.unscaled) <- "double"
  storage.mode(lambda) <- "double"
  storage.mode(alpha) <- "double"
  storage.mode(offset) <- "double"
  storage.mode(weights) <- "double"
  eta <- start$eta
  storage.mode(eta) <- "double"
  mu <- start$mu
  storage.mode(mu) <- "double"
  res <- start$residuals
  storage.mode(res) <- "double"
  
  storage.mode(c) <- "double"
  storage.mode(est.c) <- "integer"
  storage.mode(trace) <- "integer"
  adaptive <- setting$adaptive
  storage.mode(adaptive) <- "integer"
  stand <- setting$stand
  storage.mode(setting$stand) <- "integer"
  h <- as.double(1)
  storage.mode(maxIter) <- "integer"
  storage.mode(tol) <- "double"
  storage.mode(sigma2) <- "double"
  
  outputInfo <- matrix(NA, ncol = 7, nrow = nlambda)
  outputGoF <- matrix(NA, ncol = 6, nrow = nlambda)
  outputCoef <- matrix(NA, ncol = nvars, nrow = nlambda)
  outputSE <- matrix(NA, ncol = nvars, nrow = nlambda)
  outputWeight <- matrix(NA, ncol = nvars, nrow = nlambda)
  
  if (trace == 2L) cat('\r\n\r\n    Ind|Max \t      lambda \t    df \t\t     phi \t Iter \t Error')
  if(trace > 0) time0 <- Sys.time()
  
  for (rep in seq_len(nlambda)) {
    l <- lambda[rep]
    
    active <- as.vector(abs(b) > 1E-8)
    if(sum(active) == 0 | (intercept && all(which(active) == 1))) break
    
    #temp.fit <- islasso.path.fitter2(b[active], se[active], c[active], fix.c, l, alpha, x[, active, drop = FALSE], y, weights, sigma2,
    #                                cov.unscaled[active, active], nobs, intercept, tol, unpenalized[active], offset, family, maxIter)
    
    temp.fit <- if(prep$tempFamily == "gaussian"){
      .Fortran(C_islasso3, X = x[, active, drop = FALSE], y = y, n = nobs, p = sum(active), ntheta = b[active], se = se[active],
               cov = cov.unscaled[active, active], lambda = l * (!unpenalized[active]), alpha = alpha, pi = c[active], estpi = est.c,
               h = h, itmax = maxIter, tol = tol, phi = sigma2, trace = 0L, adaptive = adaptive, offset = offset, conv = integer(1),
               stand = stand, intercept = intercept, eta = eta, mu = mu, res = res, dev = double(1), weights = weights,
               hi = double(sum(active)), edf = double(1), grad = double(sum(active)))
    }else{
      .Fortran(C_islasso_glm2, X = x[, active, drop = FALSE], y = y, n = nobs, p = sum(active), ntheta = b[active], se = se[active],
               cov = cov.unscaled[active, active], lambda = l * (!unpenalized[active]), alpha = alpha, pi = c[active], estpi = est.c,
               h = h, itmax = maxIter, tol = tol, phi = sigma2, trace = 0L, adaptive = adaptive, offset = offset, conv = integer(1),
               stand = stand, intercept = intercept, eta = eta, mu = mu, dev = double(1), weights = weights,
               hi = double(sum(active)), edf = double(1), fam = as.integer(fam), link = as.integer(link), grad = double(sum(active)))
    }
                                    
    b[active] <- temp.fit$ntheta
    se[active] <- temp.fit$se
    c[active] <- temp.fit$pi
    cov.unscaled[active, active] <- temp.fit$cov
    iter <- temp.fit$itmax
    #iter2 <- temp.fit$iter2
    err <- temp.fit$tol
    #err2 <- temp.fit$err2
    #e <- temp.fit$dev
    s2 <- temp.fit$phi
    hi <- temp.fit$hi
    df <- temp.fit$edf
    conv <- temp.fit$conv
    
    # tau <- 1
    # ll <- n * log(2 * pi * tau) + e / tau
    eta <- temp.fit$eta
    mu <- temp.fit$mu
    ll <- dev <- sum(family$dev.resids(y, mu, weights))
    if ((nobs > nvars) & (prep$tempFamily %in% c("gaussian", "binomial", "poisson", "Gamma"))) ll <- family$aic(y, nobs, mu, weights, dev)
    
    aic <- ll  + 2 * df
    bic <- ll + log(nobs) * df
    aicc <- ll + 2 *  nobs * df * (df + 1) / (nobs - df - 1) + 2 * df
    ebic <- ll + (log(nobs) + 2 * g * log(nvars)) * df
    gcv <- ll / ((1 - df / nobs) ^ 2)
    gic <- ll + log(log(nobs)) * log(nvars) * df
    
    if (trace == 2L) {
      cat(
        '\r\n', formatC(rep, format = "d", width = 4, digits = 0), '|',
        formatC(nlambda, format = "d", width = 4, digits = 0), '\t ',
        formatC(l, format = "f", width = 10, digits = 4), '\t',
        formatC(df, format = "f", width = 10, digits = 4), '\t',
        formatC(s2, format = "f", width = 10, digits = 4), '\t',
        formatC(iter, format = "d", width = 4, digits = 0), '\t',
        formatC(err, format = "f", width = 5 + 3, digits = 7 + 1)
        #formatC(iter2, format = "d", width = 4, digits = 0), '\t',
        #formatC(err2, format = "f", width = 5 + 3, digits = 7 + 1)
      )
    } 
    else{
      if(trace == 1L) cat('\r ', rep, '|', nlambda, rep(' ', 20))
    }
    
    outputInfo[rep, ] <- c(l, df, s2, dev, -ll/2, iter, err)
    outputGoF[rep, ] <- c(aic, bic, aicc, ebic, gcv, gic)
    outputCoef[rep, ] <- c(b)
    outputSE[rep, ] <- c(se)
    outputWeight[rep, ] <- c(c)
  }
  
  colnames(outputInfo) <- c('lambda', "df", "phi", "deviance", "logLik", 
                            "iter", "err")
  colnames(outputGoF) <- c('AIC', 'BIC', 'AICc', "eBIC", 'GCV', "GIC")
  colnames(outputCoef) <- prep$nms
  colnames(outputSE) <- prep$nms
  colnames(outputWeight) <- prep$nms
  
  outputLinPred <- t(tcrossprod(x, na.omit(outputCoef)))
  outputFitted <- family$linkinv(outputLinPred)
  outputResid <- (rep(y, each = nrow(outputFitted)) - outputFitted) / family$mu.eta(outputFitted)
  
  if(trace > 0) cat('\n\n Executed in ', Sys.time() - time0, '(s) \n')
  output <- list(call = match.call(), Info = na.omit(outputInfo), GoF = na.omit(outputGoF), Coef = na.omit(outputCoef), 
                 SE = na.omit(outputSE), Weight = na.omit(outputWeight), Linear.predictors = outputLinPred,
                 Fitted.values = outputFitted, Residuals = outputResid,
                 Input = list(n = nobs, p = nvars, intercept = intercept, x = x, y = y, weights = weights,
                              offset = offset, family = family, alpha = alpha, unpenalized = unpenalized))
  
  return(output)
}

islasso.path.fitter2 <- function(b, se, c, fix.c, l, alpha, x, y, weights, sigma2, 
                                 cov.unscaled, n, intercept, tol, unpenalized, offset, 
                                 family, maxIter) {
  
  fn1 <- function(b, s = 1, c = .5, alpha, unpenalized) {
    lb <- length(b)
    bsc <- (b / s)
    r <- alpha * (c*(2 * pnorm(bsc, 0, 1) - 1) + (1 - c)*(2 * pnorm(bsc, 0, .00001) - 1)) / b + (1 - alpha)
    if(any(unpenalized)) r[unpenalized] <- 0
    if (length(r) > 1) {
      r <- diag(as.vector(r), nrow = lb, ncol = lb)
    }
    return(r)
  }
  fn2 <- function (b, s = 1, c = .5, alpha, unpenalized) {
    lb <- length(b)
    bsc <- (b / s)
    r <- 2 * alpha *(c * dnorm(bsc, 0, 1) + (1 - c) * dnorm(bsc, 0, .00001)) / s + (1 - alpha)
    r[is.nan(r)] <- 0
    if(any(unpenalized)) r[unpenalized] <- 0
    if (length(r) > 1) {
      r <- diag(as.vector(r), nrow = lb, ncol = lb)
    }
    return(r)
  }
  # fn0 <- function(x, y, b, s = 1, c = .5, lambda, alpha, unpenalized, family, offset, weights) {
  #   eta <- x %*% b + offset
  #   mu <- family$linkinv(eta)
  #   v <- family$variance(mu)
  #   m.e <- family$mu.eta(eta)
  #   r <- y - mu
  #   rv <- r / v * m.e * weights
  #   grad <- drop(t(rv) %*% x)
  # 
  #   bsc <- (b / s)
  #   # grad <- crossprod(x, y - drop(x %*% b))
  #   r <- alpha * (c*(2 * pnorm(bsc, 0, 1) - 1) + (1 - c)*(2 * pnorm(bsc, 0, 1E-5) - 1)) + (1 - alpha) * b
  #   if(any(unpenalized)) r[unpenalized] <- 0
  #   return(- grad + lambda * r)
  # }
  # armijo <- function(x, y, b, se, c, l, alpha, unpenalized, family, offset, weights, cov.unscaled, cov.new, h = 1, s2){
  #   grad0 <- drop(crossprod(fn0(x, y, b, se, c, l, alpha, unpenalized, family, offset, weights)))
  #   grad.new <- .Machine$double.xmax
  #   while(grad.new > grad0 & h > .1){
  #     covar <- cov.unscaled + h * (cov.new - cov.unscaled)
  #     cov.scaled <- s2 * covar
  #     se <- sqrt(diag(cov.scaled))
  #     grad.new <- drop(crossprod(fn0(x, y, b, se, c, l, alpha, unpenalized, family, offset, weights)))
  #     h <- .5 * h
  #   }
  #   return(list(cov.unscaled = covar, cov.scaled = cov.scaled, se = se, grad = grad.new, h = 2 * h))
  # }
  
  conv <- 0L
  iter <- iter2 <- 1
  err2 <- 1
  tmp <- b
  tmp2 <- se
  err.old <- 1
  s2 <- sigma2

  eta <- drop(x %*% b) + offset
  mu <- family$linkinv(eta)
  varmu <- family$variance(mu)
  mueta <- family$mu.eta(eta)
  res <- (y - mu) / mueta
  
  while(err2 > (tol*10)) {
    err <- 1

    if(!fix.c) c <- 2 * binomial()$linkinv(abs(b / se)) - 1
    z <- (eta - offset) + res
    w <- weights * (mueta^2) / varmu
    Xw <- x * w
    XtX <- crossprod(Xw, x)
    XtY <- crossprod(Xw, z)
    
    while (err > tol) {
      inV <- XtX + l * fn1(b = b, s = se, c = c, alpha = alpha, unpenalized = unpenalized)
      #inV <- try(chol2inv(chol(inV)), silent = TRUE)
      b <- qr.solve(inV, XtY)
      if(inherits(inV, "try-error")) {
        # b <- tmp 
        err <- err.old
        conv <- 1L
        warning(paste0("XtX + l*P' not invertible at lambda ", l, ". Safe exit!"))
        break()
      }
      #b <- inV %*% XtY
      
      err.old <- err <- max(abs(tmp - b))
      # if(err <= 10 ^ -(digit + 0)) print(data.frame(it = iter, err = err, t(b[1:10])))
      tmp <- b
      iter <- iter + 1
      if(iter > maxIter) break
    }
    if(conv == 1) break()
    
    eta <- drop(x %*% b) + offset
    mu <- family$linkinv(eta)
    varmu <- family$variance(mu)
    mueta <- family$mu.eta(eta)
    res <- (y - mu) / mueta
    e <- sum(weights * (res^2))
    # e <- sum(weights * ((y - mu)^2))
    inVg <- XtX + l * fn2(b = b, s = se, c = c, alpha = alpha, unpenalized = unpenalized)
    invH <- chol2inv(chol(inVg))
    tempMat <- invH %*% XtX
    hi <- diag(tempMat)
    df <- sum(hi)
    if(sigma2 <= 0) s2 <- e / (n - df)
    
    cov.new <- tempMat %*% invH
    
    # temp.covar <- armijo(x, y, b, se, c, l, alpha, unpenalized, family, offset, weights, cov.unscaled, cov.new, h = 1, s2)
    # cov.unscaled <- temp.covar$cov.unscale
    # cov.scaled <- temp.covar$cov.scale
    # se <- temp.covar$se
    
    # temp.covar <- list(h = .25, grad = drop(crossprod(fn0(x, y, b, se, c, l, alpha, unpenalized, family, offset, weights))),
    #                    n.frob = norm(cov.new - cov.unscaled, "F"))
    cov.unscaled <- cov.unscaled + .1 * (cov.new - cov.unscaled)
    cov.scaled <- s2 * cov.unscaled
    se <- sqrt(diag(cov.scaled))
    
    err2 <- max(abs(tmp2 - se))
    
    # print(data.frame(it = iter2, err = err, err2 = err2, s2 = s2, df = df, e = e))
    # print(data.frame(it = iter2, err = err, err2 = err2, s2 = s2, df = df, e = e,
    #                  h = temp.covar$h, grad = temp.covar$grad, n.frob = temp.covar$n.frob))
    tmp2 <- se
    iter2 <- iter2 + 1
    if(iter2 > maxIter) break
  }
  
  list(b = b, se = se, c = c, cov.unscaled = cov.unscaled, iter = iter, 
       err = err, iter2 = iter2, err2 = err2,  e = e, s2 = s2, hi = hi, 
       df = df, conv = conv, eta = eta, mu = mu, res = res)
}

GoF.islasso.path <- function(object, plot = TRUE, ...){
  lambda.seq <- object$Info[,"lambda"]
  nlambda <- length(lambda.seq)
  gof.name <- c("AIC", "BIC", "AICc", "eBIC", "GCV", "GIC")
  gof <- object$GoF
  id.min <- apply(gof, 2, which.min)
  lambda.min <- lambda.seq[id.min]
  names(lambda.min) <- gof.name
  
  if(plot) {
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mfrow = c(2, 3))
    for(i in seq_len(length(gof.name))) {
      plot(log(lambda.seq), gof[, i], lwd = 1.5, type = "l", xlab = "log Lambda", ylab = gof.name[i])
      abline(v = log(lambda.min[i]), lty = 3, lwd = 2)
    }
  }
  
  return(list(gof = gof, minimum = id.min, lambda.min = lambda.min))
}


# islasso.path.fit.lm <- function(prep, start, lambda, fam, link){
#   
#   cov.unscaled <- start$covar
#   tmp2 <- se <- start$se
#   tmp <- b <- start$beta + .01
#   
#   c <- prep$setting$c
#   fix.c <- prep$setting$fix.c
#   tol <- prep$setting$tol
#   itmax <- prep$setting$itmax
#   trace <- prep$setting$trace
#   sigma2 <- prep$setting$sigma2
#   nlambda <- length(lambda)
#   alpha <- prep$alpha
#   nvars <- prep$nvars
#   nobs <- prep$nobs
#   intercept <- prep$intercept
#   x <- prep$X
#   y <- prep$y
#   weights <- prep$weights
#   offset <- prep$offset
#   unpenalized <- prep$unpenalized
#   family <- prep$family
#   
#   xw <- x * weights
#   XtX <- crossprod(xw, x)
#   XtY <- crossprod(xw, y)
#   
#   outputInfo <- matrix(NA, ncol = 13, nrow = nlambda)
#   outputCoef <- matrix(NA, ncol = nvars, nrow = nlambda)
#   outputSE <- matrix(NA, ncol = nvars, nrow = nlambda)
#   outputWeight <- matrix(NA, ncol = nvars, nrow = nlambda)
#   
#   if (trace == 2L) cat('\r\n\r\n    Ind|Max \t      lambda \t    df \t\t     phi \t Iter1 \t ErrorB \t Iter2 \t Error SE')
#   time0 <- Sys.time()
#   
#   for (rep in seq_len(nlambda)) {
#     l <- lambda[rep]
#     
#     active <- as.vector(abs(b) > 1E-6)
#     if(sum(active) == 0 | (intercept && all(which(active) == 1))) break
#     
#     temp.fit <- islasso.path.fitter(b[active], se[active], c[active], fix.c, l, alpha, XtX[active, active], XtY[active,], 
#                                     xw[, active, drop = FALSE], y, cov.unscaled[active, active], nobs, intercept, tol, 
#                                     unpenalized, offset, weights, sigma2)
#     b[active] <- temp.fit$b
#     se[active] <- temp.fit$se
#     c[active] <- temp.fit$c
#     cov.unscaled[active, active] <- temp.fit$cov.unscaled
#     iter <- temp.fit$iter
#     iter2 <- temp.fit$iter2
#     err <- temp.fit$err
#     err2 <- temp.fit$err2
#     e <- temp.fit$e
#     s2 <- temp.fit$s2
#     hi <- temp.fit$hi
#     df <- temp.fit$df
#     conv <- temp.fit$conv
#     
#     tau <- 1
#     ll <- n * log(2 * pi * tau) + e / tau
#     aicc <- ll + 2 *  n * df * (df + 1) / (n - df - 1) + 2 * df
#     bic <- ll + log(n) * df
#     gcv <- ll / ((1 - df / n) ^ 2)
#     AIC <- ll  + 2 * df
#     
#     if (trace == 2L) {
#       cat(
#         '\r\n', formatC(rep, format = "d", width = 4, digits = 0), '|',
#         formatC(nlambda, format = "d", width = 4, digits = 0), '\t ',
#         formatC(l, format = "f", width = 10, digits = 4), '\t',
#         formatC(df, format = "f", width = 10, digits = 4), '\t',
#         formatC(s2, format = "f", width = 10, digits = 4), '\t',
#         formatC(iter, format = "d", width = 4, digits = 0), '\t',
#         formatC(err,, format = "f", width = 5 + 3, digits = 5 + 1), '\t',
#         formatC(iter2, format = "d", width = 4, digits = 0), '\t',
#         formatC(err2, format = "f", width = 5 + 3, digits = 5 + 1)
#       )
#     } 
#     else{
#       if(trace == 1L) cat('\r ', rep, '|', nlambda, rep(' ', 20))
#     }
#     
#     outputInfo[rep, ] <- c(l, df, s2, e, ll, iter, iter2, err, err2, aicc, AIC, bic, gcv)
#     outputCoef[rep, ] <- c(b)
#     outputSE[rep, ] <- c(se)
#     outputWeight[rep, ] <- c(c)
#   }
#   
#   colnames(outputInfo) <- c('lambda', "df", "phi", "deviance", "logLik", 
#                             "iterB", "iterSE", "errB", "errSE", 
#                             'AICc', 'AIC', 'BIC', 'GCV')
#   colnames(outputCoef) <- prep$nms
#   colnames(outputSE) <- prep$nms
#   colnames(outputWeight) <- prep$nms
#   
#   outputLinPred <- t(tcrossprod(x, na.omit(outputCoef)))
#   outputFitted <- outputLinPred
#   outputResid <- rep(y, each = nrow(outputFitted)) - outputFitted
#   
#   cat('\n\n Executed in ', Sys.time() - time0, '(s) \n')
#   output <- list(call = match.call(), Info = na.omit(outputInfo), Coef = na.omit(outputCoef), 
#                  SE = na.omit(outputSE), Weight = na.omit(outputWeight), Linear.predictors = outputLinPred,
#                  Fitted.values = outputFitted, Residuals = outputResid,
#                  Input = list(n = nobs, p = nvars, intercept = intercept, x = x, y = y, alpha = alpha, unpenalized = unpenalized))
#   
#   return(output)
# }
# 
# islasso.path.fitter <- function(b, se, c, fix.c, l, alpha, XtX, XtY, newx, y, 
#                                 cov.unscaled, n, intercept, tol, unpenalized, 
#                                 offset, weights, sigma2) {
#   fn1 <- function(b, s = 1, c = .5, alpha, unpenalized) {
#     lb <- length(b)
#     bsc <- (b / s)
#     r <- alpha * (c*(2 * pnorm(bsc, 0, 1) - 1) + (1 - c)*(2 * pnorm(bsc, 0, 1E-5) - 1)) / b + (1 - alpha)
#     if(any(unpenalized)) r[unpenalized] <- 0
#     if (length(r) > 1) {
#       r <- diag(as.vector(r), nrow = lb, ncol = lb)
#     }
#     return(r)
#   }
#   fn2 <- function (b, s = 1, c = .5, alpha, unpenalized) {
#     lb <- length(b)
#     bsc <- (b / s)
#     r <- alpha *(c*(2 * dnorm(bsc, 0, 1) / s) + (1 - c)*(2 * dnorm(bsc, 0, 1E-5) / s)) + (1 - alpha)
#     r[is.nan(r)] <- 0
#     if(any(unpenalized)) r[unpenalized] <- 0
#     if (length(r) > 1) {
#       r <- diag(as.vector(r), nrow = lb, ncol = lb)
#     }
#     return(r)
#   }
#   
#   conv <- 0L
#   iter <- iter2 <- 1
#   err2 <- 1
#   tmp <- b
#   tmp2 <- se
#   err.old <- 1
#   s2 <- sigma2
#   
#   while(err2 > tol) {
#     err <- 1
#     
#     if(!fix.c) c <- 2*binomial()$linkinv(abs(b / se)) - 1
#     
#     while (err > tol) {
#       inV <- XtX + l * fn1(b = b, s = se, c = c, alpha = alpha, unpenalized = unpenalized)
#       b <- temp.check <- try(solve(inV, XtY), silent = TRUE)
#       if(inherits(temp.check, "try-error")) {
#         b <- tmp 
#         err <- err.old
#         conv <- 1L
#         warning(paste0("XtX + l*P' not invertible at lambda ", l, ". Safe exit!"))
#         break()
#       }
#       
#       err.old <- err <- max(abs(tmp - b))
#       # if(err <= 10 ^ -(digit + 0)) print(data.frame(it = iter, err = err, t(b[1:10])))
#       tmp <- b
#       iter <- iter + 1
#     }
#     if(conv == 1) break()
#     eta <- newx %*% b + offset
#     res <- y - eta
#     e <- sum(weights * (res^2))
#     inVg <- XtX + l * fn2(b = b, s = se, c = c, alpha = alpha, unpenalized = unpenalized)
#     invH <- chol2inv(chol(inVg))
#     tempMat <- invH %*% XtX
#     hi <- diag(invH %*% XtX)
#     df <- sum(hi)
#     
#     cov.unscaled <- cov.unscaled + .25 * (tempMat %*% invH - cov.unscaled)
#     
#     if(sigma2 <= 0) s2 <- e / (n - df)
#     cov.scaled <- s2 * cov.unscaled
#     se <- sqrt(diag(cov.scaled))
#     err2 <- max(abs(tmp2 - se))
#     # print(data.frame(it = iter2, err2 = err2, s2 = s2, df = df, t(se[1:10])))
#     tmp2 <- se
#     iter2 <- iter2 + 1
#   }
#   
#   list(b = b, se = se, c = c, cov.unscaled = cov.unscaled, iter = iter, 
#        err = err, iter2 = iter2, err2 = err2,  e = e, s2 = s2, hi = hi, 
#        df = df, conv = conv, eta = eta, mu = eta, res = res)
# }
