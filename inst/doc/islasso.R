## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.align = "center", message = FALSE, cache = TRUE, comment = ">", tidy = TRUE, warning = FALSE)
# the code in this chunk enables us to truncate the print output for each
# chunk using the `out.lines` option
# save the built-in output hook
hook_output <- knitr::knit_hooks$get("output")

# set a new output hook to truncate text output
knitr::knit_hooks$set(output = function(x, options) {
  if (!is.null(n <- options$out.lines)) {
    x <- xfun::split_lines(x)
    if (length(x) > n) {
        
      # truncate the output
      x <- c(head(x, n), "....\n")
    }
    x <- paste(x, collapse = "\n")
  }
  hook_output(x, options)
})

## ---- dont-eval, eval=FALSE---------------------------------------------------
#  islasso(formula, family = gaussian, lambda, alpha = 1, data, weights, subset,
#          offset, unpenalized, contrasts = NULL, control = is.control())

## ---- dont-eval1, eval=FALSE--------------------------------------------------
#  islasso.fit(X, y, family = gaussian, lambda, alpha = 1, intercept = FALSE,
#      weights = NULL, offset = NULL, unpenalized = NULL, control = is.control())

## -----------------------------------------------------------------------------
library(lars)
library(glmnet)

data("diabetes", package = "lars")

a1 <- with(diabetes, cv.glmnet(x2, y))
n <- nrow(diabetes)
a1$lambda.min * n

b <- drop(coef(a1, "lambda.min", exact = TRUE))
length(b[b != 0])

## ---- fold.source = TRUE------------------------------------------------------
names(b[b != 0])

## -----------------------------------------------------------------------------
tail(b[b != 0])

## -----------------------------------------------------------------------------
library(islasso)
out <- islasso(y ~ x2, data = diabetes, lambda = a1$lambda.min * n)

## -----------------------------------------------------------------------------
summary(out, pval = 0.10)

## -----------------------------------------------------------------------------
sum(out$internal$hi)

## -----------------------------------------------------------------------------
lmb.bic <- aic.islasso(out, method = "bic", interval = c(1, 100))
out1 <- update(out, lambda = lmb.bic)
summary(out1, pval = .05)

## -----------------------------------------------------------------------------
out2 <- update(out, alpha = .5)
lmb.bic <- aic.islasso(out2, method = "bic", interval = c(1, 100))
out3 <- update(out, lambda = lmb.bic)
summary(out3, pval = .05)

