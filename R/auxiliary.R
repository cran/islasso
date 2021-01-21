format.perc <- function (probs, digits){
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
