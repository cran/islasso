confint.islasso <- function (object, parm, level = 0.95, type.ci = "wald", trace = TRUE, ...) {
  type.ci <- match.arg(type.ci)
  if(type.ci != "wald") stop("Only Wald-type confidence intervals are implemented yet!")
  pnames <- names(B0 <- coef(object))
  if (missing(parm)) 
    parm <- seq_along(pnames)
  else if (is.character(parm)) 
    parm <- match(parm, pnames, nomatch = 0L)
  
  if(type.ci == "wald"){
    alpha <- (1 - level)/2
    a <- c(alpha, 1 - alpha)
    fac <- qnorm(a)
    pct <- paste(round(100 * a, 1), "%")
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(pnames[parm], pct))
    ses <- object$se
    ci[] <- B0[parm] + ses[parm] %o% fac
  }
  
  attr(ci, "coefficients") <- B0[parm]
  attr(ci, "level") <- level
  class(ci) <- "confint.islasso"
  ci
}

print.confint.islasso <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  print.default(format(round(x, 6), digits = digits), print.gap = 2, quote = FALSE)
  cat("\n")
  invisible(x)
}

plot.confint.islasso <- function(x, parm, ...){
  beta <- attr(x, "coefficients")
  if(missing(parm)) parm <- seq_along(beta)
  nms <- names(beta)
  p <- length(parm)
  L <- list(...)
  if(is.null(L$ylab)) L$ylab <- "Confidence intervals"
  if(is.null(L$xlab)) L$xlab <- ""
  if(is.null(L$srt)) L$srt <- 45
  if(is.null(L$lwd)) L$lwd <- 1
  if(is.null(L$col)) L$col <- "grey"
  if(is.null(L$cex)) L$cex <- 1
  if(is.null(L$text.cex)) L$text.cex <- 1
  if(is.null(L$pch)) L$pch <- 16
  
  plot(1:p, seq(min(x[parm,]), max(x[parm,]), l=p), type = "n", axes = FALSE, ylab = L$ylab, xlab=L$xlab)
  segments(1:p, x[parm,1], 1:p, x[parm,2], col=L$col, lwd=L$lwd)
  points(beta[parm], pch=16, cex=L$cex)
  abline(h=0, lty=3)
  axis(2); axis(2, at=c(-1000, 1000))
  axis(1, at=1:p, cex.axis=1, labels=FALSE); axis(1, at=c(-1000, 1000))
  if(is.null(L$text.y)) L$text.y <- par("usr")[3] * 1.1
  text(1:p, L$text.y, labels=nms[parm], srt=L$srt, pos=1, xpd=TRUE, cex=L$text.cex)
}

makeHyp <- function (cnames, hypothesis, rhs = NULL) {
  parseTerms <- function(terms) {
    component <- gsub("^[-\\ 0-9\\.]+", "", terms)
    component <- gsub(" ", "", component, fixed = TRUE)
    component
  }
  stripchars <- function(x) {
    x <- gsub("\\n", " ", x)
    x <- gsub("\\t", " ", x)
    x <- gsub(" ", "", x, fixed = TRUE)
    x <- gsub("*", "", x, fixed = TRUE)
    x <- gsub("-", "+-", x, fixed = TRUE)
    x <- strsplit(x, "+", fixed = TRUE)[[1]]
    x <- x[x != ""]
    x
  }
  char2num <- function(x) {
    x[x == ""] <- "1"
    x[x == "-"] <- "-1"
    as.numeric(x)
  }
  constants <- function(x, y) {
    with.coef <- unique(unlist(sapply(y, function(z) which(z == parseTerms(x)))))
    if (length(with.coef) > 0) x <- x[-with.coef]
    x <- if (is.null(x)) 0 else sum(as.numeric(x))
    if (any(is.na(x))) stop("The hypothesis \"", hypothesis, "\" is not well formed: contains bad coefficient/variable names.")
    x
  }
  coefvector <- function(x, y) {
    rv <- gsub(" ", "", x, fixed = TRUE) == parseTerms(y)
    if (!any(rv)) return(0)
    if (sum(rv) > 1) stop("The hypothesis \"", hypothesis, "\" is not well formed.")
    rv <- sum(char2num(unlist(strsplit(y[rv], x, fixed = TRUE))))
    if (is.na(rv)) stop("The hypothesis \"", hypothesis, "\" is not well formed: contains non-numeric coefficients.")
    rv
  }
  if (!is.null(rhs)) rhs <- rep(rhs, length.out = length(hypothesis))
  if (length(hypothesis) > 1) 
    return(rbind(Recall(cnames, hypothesis[1], rhs[1]), Recall(cnames, hypothesis[-1], rhs[-1])))
  cnames_symb <- sapply(c("@", "#", "~"), function(x) length(grep(x, cnames)) < 1)
  if (any(cnames_symb)) {
    cnames_symb <- head(c("@", "#", "~")[cnames_symb], 1)
    cnames_symb <- paste(cnames_symb, seq_along(cnames), cnames_symb, sep = "")
    hypothesis_symb <- hypothesis
    for (i in order(nchar(cnames), decreasing = TRUE)) 
      hypothesis_symb <- gsub(cnames[i], cnames_symb[i], hypothesis_symb, fixed = TRUE)
  }
  else {
    stop("The hypothesis \"", hypothesis, "\" is not well formed: contains non-standard coefficient names.")
  }
  lhs <- strsplit(hypothesis_symb, "=", fixed = TRUE)[[1]]
  if (is.null(rhs)) {
    if (length(lhs) < 2) 
      rhs <- "0"
    else if (length(lhs) == 2) {
      rhs <- lhs[2]
      lhs <- lhs[1]
    }
    else stop("The hypothesis \"", hypothesis, "\" is not well formed: contains more than one = sign.")
  }
  else {
    if (length(lhs) < 2) 
      as.character(rhs)
    else stop("The hypothesis \"", hypothesis, "\" is not well formed: contains a = sign although rhs was specified.")
  }
  lhs <- stripchars(lhs)
  rhs <- stripchars(rhs)
  rval <- sapply(cnames_symb, coefvector, y = lhs) - sapply(cnames_symb, coefvector, y = rhs)
  rval <- c(rval, constants(rhs, cnames_symb) - constants(lhs, cnames_symb))
  names(rval) <- c(cnames, "*rhs*")
  rval
}

printHyp <- function (L, b, nms) {
  nomi <- apply(L, 1, function(l){
    id <- which(l != 0)
    temp <- paste(l[id], nms[id], sep = "*", collapse = " + ")
    if(nchar(temp) > min(getOption("width"), 50)) temp <- paste(strtrim(temp, min(getOption("width"), 50)), "...")
    temp
  })
  paste0(nomi, " = ", b)
}

anova.islasso <- function (object, A, b = NULL, ci, ...) {
  # beta <- if(missing(ci)) coef(object) else attr(ci, "unbiased")
  beta <- coef(object)
  V <- vcov(object)
  nms <- names(beta)
  
  if (any(aliased <- is.na(beta))) stop("there are aliased coefficients in the model")
  beta <- beta[!aliased]
  if (is.null(beta)) 
    stop(paste("there is no coef() method for models of class", paste(class(object), collapse = ", ")))
  if (is.character(A)) {
    L <- makeHyp(nms, A, b)
    if (is.null(dim(L))) L <- t(L)
    b <- L[, NCOL(L)]
    L <- L[, -NCOL(L), drop = FALSE]
    rownames(L) <- A
  }
  else {
    L <- if (is.null(dim(A))) t(A) else A
    if (is.null(b)) b <- rep(0, nrow(L))
  }
  q <- NROW(L)
  
  value.hyp <- L %*% beta - b
  vcov.hyp <- L %*% V %*% t(L)
  
  rval <- matrix(NA, nrow = q + 1L, ncol = 4L)
  colnames(rval) <- c("Estimate", "Std. Error", "Chisq", "Pr(>Chisq)")
  rownames(rval) <- c(printHyp(L, b, nms), "Overall")
  rval[1:q, 1L] <- value.hyp
  rval[1:q, 2L] <- sqrt(diag(vcov.hyp))
  rval[1:q, 3L] <- value.hyp^2 / diag(vcov.hyp)
  rval[1:q, 4L] <- pchisq(rval[1:q, 3L], 1L, lower.tail = FALSE)
  
  if (q > 1) {
    statistic2 <- as.vector(t(value.hyp) %*% solve(vcov.hyp) %*% value.hyp)
    pval2 <- pchisq(statistic2, q, lower.tail = FALSE)
    rval[q + 1L, 3:4] <- c(statistic2, pval2)
  } 
  else rval <- rval[-c(q + 1L), , drop = FALSE]
  
  ic <- NULL
  if(!missing(ci)) {
    ic <- t(apply(L, 1, function(a) cislasso(object, a = a, ci = ci)))
    colnames(ic) <- colnames(ci)
    rownames(ic) <- rownames(rval)[1:q]
  }
  
  object$anova <- list(result = as.data.frame(rval), ci = ic)
  class(object) <- c("anova.islasso", class(object))
  return(object)
}

print.anova.islasso <- function(x, digits = max(3, getOption("digits") - 3), ...){
  cat("\n\t", "Simultaneous Tests for General Linear Combination\n\n")
  call <- x$call
  if (!is.null(call)) {
    cat("Fit: ")
    print(call)
    cat("\n")
  }
  pq <- x$anova
  
  cat("Linear Hypothesis:\n")
  printCoefmat(pq$result, digits = digits, has.Pvalue = TRUE, P.values = TRUE, 
               eps.Pvalue = .Machine$double.eps, na.print = "")
  cat("\n")
  if(!is.null(pq$ci)){
    cat("Confidence Interval Estimation for Linear Hypothesis:\n")
    print.default(format(round(pq$ci, 6), digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
  }
  invisible(x)
}

cislasso <- function(object, a, ci){
  V <- vcov(object)
  est <- coef(object)
  if(missing(a)) stop("vector of linear combination is missing")
  if(missing(ci)) stop("confidence intervals are missing")
  th <- sum(a*est)
  low <- ci[, 1]
  up <- ci[, 2]
  difmin <- a*est - pmin(a*low, a*up)
  difmax <- a*est - pmax(a*low, a*up)
  Cl <- Cu <- 0
  R <- cov2cor(V)
  r <- R[col(R) > row(R)]
  DMIN <- outer(difmin, difmin)
  dmin <- DMIN[col(R) > row(R)]
  Cl <- 2*sum(r*dmin)
  DMAX <- outer(difmax, difmax)
  dmax <- DMAX[col(R) > row(R)]
  Cu <- 2*sum(r*dmax)
  L <- th - sqrt(sum(difmin^2) + Cl)
  U <- th + sqrt(sum(difmax^2) + Cu)
  r <- c(L, U)
  r
}

ci.fitted.islasso <- function(object, newx, ci = NULL, type.ci = "wald", conf.level=.95, only.ci = FALSE){
  type.ci <- match.arg(type.ci)
  if(type.ci != "wald") stop("Only Wald-type confidence intervals are implemented yet!")
  if(missing(newx)) newx <- model.matrix(object)
  if(is.null(ci)) ci <- confint.islasso(object, level = conf.level, type.ci = type.ci, trace = FALSE)
  if(only.ci) return(ci)
  n <- nrow(newx)
  ris <- t(sapply(1:n, function(i) cislasso(object, newx[i,], ci)))
  colnames(ris) <- colnames(ci)
  rownames(ris) <- rownames(newx)
  ris
} 

