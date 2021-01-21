confint.islasso <- function(object, parm, level = 0.95, ...){
  cf <- coef(object)
  pnames <- names(cf)
  if(missing(parm)) parm <- seq_along(cf)
  
  alpha <- 1 - level
  a <- alpha/2
  a <- c(a, 1 - a)
  fac <- qnorm(a)
  pct <- format.perc(a, 3)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(pnames[parm], pct))
  ses <- object$se
  ci[] <- cf[parm] + ses[parm] %o% fac
 
  attr(ci, "coefficients") <- coef(object)[parm]
  attr(ci, "level") <- level
  class(ci) <- "confint.islasso"
  
  return(ci)
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

anova.islasso <- function(object, A, b, ci, ...){
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
  
  ic <- if(!missing(ci)) matrix(NA, k, 2, dimnames = list(NULL, colnames(ci))) else NULL
  
  if(mpc){
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
    if(!missing(ci)) ic[i, ] <- cislasso(object, a = A[i, ], ci = ci)
  }
  if(mpc) nomi <- c(nomi, "Overall")
  
  object$anova <- list(A=A, b=b, coefficients=beta_new, vcov=V_new, k=k, nms=nomi, mpc=mpc, 
                       tstat=sqrt(statistic), pvalues=pval, tstat2=statistic2, pvalues2=pval2, ci=ic)
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
  sig <- .Machine$double.eps
  if(pq$mpc) mtests <- rbind(mtests, c(0, 0, pq$tstat2, pq$pvalues2))
  colnames(mtests) <- c("Estimate", "Std. Error", "chi2 value", "Pr(>|chi2|)")
  rownames(mtests) <- paste(pq$nms, "=", format(pq$b, digits=3, trim=TRUE))
  if(pq$mpc) rownames(mtests)[nrow(mtests)] <- pq$nms[nrow(mtests)]
  cat("Linear Hypotheses:\n")
  if(pq$mpc){
    printCoefmat2(mtests, digits = digits, has.Pvalue = TRUE, P.values = TRUE, eps.Pvalue = sig)
  }else{
    printCoefmat(mtests, digits = digits, has.Pvalue = TRUE, P.values = TRUE, eps.Pvalue = sig)
  }
  cat("\n")
  if(!is.null(pq$ci)){
    cat("Confidence Interval Estimation for Linear Hypotheses:\n")
    rownames(pq$ci) <- pq$nms[1:pq$k]
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
  U <- th + sqrt(sum(difmin^2) + Cu)
  r <- c(L, U)
  r
}

ci.fitted.islasso <- function(object, newx, ci = NULL, conf.level=.95, only.ci = FALSE){
  if(missing(newx)) newx <- model.matrix(object)
  if(is.null(ci)) ci <- confint.islasso(object, level = conf.level)
  if(only.ci) return(ci)
  n <- nrow(newx)
  ris <- t(sapply(1:n, function(i) cislasso(object, newx[i,], ci)))
  colnames(ris) <- colnames(ci)
  rownames(ris) <- rownames(newx)
  ris
} 
