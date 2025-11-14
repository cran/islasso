## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----dont-eval, eval=FALSE----------------------------------------------------
# islasso(formula, family = gaussian, lambda, alpha = 1, data, weights, subset,
#         offset, unpenalized, contrasts = NULL, control = is.control())

## ----dont-eval-1, eval=FALSE--------------------------------------------------
# islasso.path(formula, family = gaussian, lambda = NULL, nlambda = 100,
#              lambda.min.ratio = ifelse(nobs < nvars, 1E-3, 1E-05), alpha = 1,
#              data, weights, subset, offset, unpenalized, contrasts = NULL, control = is.control())

## ----dont-eval1, eval=FALSE---------------------------------------------------
# islasso.fit(X, y, family = gaussian(), lambda, alpha = 1, intercept = FALSE,
#             weights = NULL, offset = NULL, unpenalized = NULL, control = is.control())

## ----dont-eval1-1, eval=FALSE-------------------------------------------------
# islasso.path.fit(X, y, family = gaussian(), lambda, nlambda, lambda.min.ratio,
#                  alpha = 1, intercept = FALSE, weights = NULL, offset = NULL,
#                  unpenalized = NULL, control = is.control())

## -----------------------------------------------------------------------------
library(islasso)

data("Prostate", package = "islasso")
x <- model.matrix(lpsa ~ ., data = Prostate)[, -1]
y <- Prostate$lpsa
a1 <- cv.glmnet(x, y, family = "gaussian")
n <- nrow(Prostate)
a1$lambda.min * n

b <- drop(coef(a1, "lambda.min", exact = TRUE))
length(b[b != 0])

## ----fold.source = TRUE-------------------------------------------------------
names(b[b != 0])

## -----------------------------------------------------------------------------
tail(b[b != 0], n = 3)

## -----------------------------------------------------------------------------
out <- islasso.path(lpsa ~ ., data = Prostate, nlambda = 50L, family = gaussian())
out

## -----------------------------------------------------------------------------
lmb.best <- GoF.islasso.path(out)
lmb.best$lambda.min

## ----plot1, fig.cap="Regularization path for coefficients, standard errors and gradients."----
p1 <- plot(out, yvar = "coefficients")
p2 <- plot(out, yvar = "se")
p3 <- plot(out, yvar = "gradient")
gridExtra::grid.arrange(p1, p2, p3, ncol = 1L)

## -----------------------------------------------------------------------------
lambda.aic <- lmb.best$lambda.min["AIC"]
out2 <- islasso(lpsa ~ ., data = Prostate, lambda = lambda.aic, family = gaussian())
out2

## -----------------------------------------------------------------------------
summary(out2, pval = 0.10)

## -----------------------------------------------------------------------------
sum(out2$internal$hi)

## -----------------------------------------------------------------------------
lambda.aic2 <- aic.islasso(out2, method = "AIC", interval = c(.1, 50))
out3 <- update(out2, lambda = lambda.aic2)
summary(out3, pval = .10)

## -----------------------------------------------------------------------------
# update the islasso path to fit an elastic-net model
out4 <- update(out, alpha = .5)
out4

# some diagnostic plot
p4 <- plot(out4, yvar = "coefficients")
p5 <- plot(out4, yvar = "se")
gridExtra::grid.arrange(p4, p5, ncol = 1L)

# select the best tuning parameter
lmb.best2 <- GoF.islasso.path(out4)
lmb.best2$lambda.min

# fit a new islasso model with elastic-net penalty
lambda.aic3 <- lmb.best2$lambda.min["AIC"]
out5 <- update(out2, alpha = .5, lambda = lambda.aic3)
summary(out5, pval = .10)

# or select the best tuning parameter using AIC with an islasso object
lambda.aic4 <- aic.islasso(out5, method = "AIC", interval = c(.1, 100))
out6 <- update(out5, lambda = lambda.aic4)
summary(out6, pval = .10)

