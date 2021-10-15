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

## ---- dont-eval-1, eval=FALSE-------------------------------------------------
#  islasso.path(formula, family = gaussian, lambda = NULL, nlambda = 100,
#               lambda.min.ratio = ifelse(nobs < nvars, 1E-3, 1E-05), alpha = 1,
#               data, weights, subset, offset, unpenalized, contrasts = NULL, control = is.control())

## ---- dont-eval1, eval=FALSE--------------------------------------------------
#  islasso.fit(X, y, family = gaussian(), lambda, alpha = 1, intercept = FALSE,
#              weights = NULL, offset = NULL, unpenalized = NULL, control = is.control())

## ---- dont-eval1-1, eval=FALSE------------------------------------------------
#  islasso.path.fit(X, y, family = gaussian(), lambda, nlambda, lambda.min.ratio,
#                   alpha = 1, intercept = FALSE, weights = NULL, offset = NULL,
#                   unpenalized = NULL, control = is.control())

## -----------------------------------------------------------------------------
library(islasso)

data("diabetes", package = "islasso")

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
out <- islasso.path(y ~ x2, data = diabetes, nlambda = 30L)
out

## -----------------------------------------------------------------------------
lmb.best <- GoF.islasso.path(out)
lmb.best$lambda.min

## -----------------------------------------------------------------------------
par(mfrow = c(1, 3))
plot(out, yvar = "coefficients")
plot(out, yvar = "se")
plot(out, yvar = "gradient")

## -----------------------------------------------------------------------------
lambda.bic <- lmb.best$lambda.min["BIC"]
out2 <- islasso(y ~ x2, data = diabetes, lambda = lambda.bic)
out2

## -----------------------------------------------------------------------------
summary(out2, pval = 0.10)

## -----------------------------------------------------------------------------
sum(out2$internal$hi)

## -----------------------------------------------------------------------------
lambda.bic2 <- aic.islasso(out2, method = "BIC", interval = c(1, 100))
out3 <- update(out2, lambda = lambda.bic2)
summary(out3, pval = .10)

## -----------------------------------------------------------------------------
# update the islasso path to fit an elastic-net model
out4 <- update(out, alpha = .5)
out4

# some diagnostic plot
par(mfrow = c(1, 2))
plot(out4, yvar = "coefficients")
plot(out4, yvar = "se")

# select the best tuning parameter
lmb.best2 <- GoF.islasso.path(out4)
lmb.best2$lambda.min

# fit a new islasso model with elastic-net penalty
lambda.bic3 <- lmb.best2$lambda.min["BIC"]
out5 <- update(out2, alpha = .5, lambda = lambda.bic3)
summary(out5, pval = .10)

# or select the best tuning parameter using BIC with an islasso object
lambda.bic4 <- aic.islasso(out5, method = "BIC", interval = c(1, 100))
out6 <- update(out5, lambda = lambda.bic4)
summary(out6, pval = .10)

