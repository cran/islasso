\name{is.control}
\alias{is.control}
\title{Auxiliary for controlling islasso model fitting}
\usage{
is.control(sigma2 = -1, tol = 1e-05, itmax = 500, stand = FALSE,
  trace = 0, nfolds = 5, seed=NULL, adaptive = FALSE, 
  b0 = NULL, V0 = NULL, c = 0.5)
}
\arguments{
\item{sigma2}{optional. The fixed value of dispersion parameter. If -1 (default) it is estimated from the data}

\item{tol}{tollerance value to declare convergence, dafault to 1e-5}

\item{itmax}{maximum number of iterations, default to 500}

\item{stand}{if TRUE, the covariates are standardized prior to fitting the model. However the coefficients are always returned on the original scale.}

\item{trace}{Should the iterative procedure be printed? 0: no printing, 1 = compact printing, 2 = enlarged printing, 3 = compact printing including Fisher scoring information (only used in glm family). }

\item{nfolds}{if \code{lambda} is unspecified in \code{islasso}, the number of folds to be used to perform cross valdation. Default to 5, and \code{nfolds>2} is allowed. Although nfolds can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. \code{nfolds} is ignored if \code{lambda} is supplied.}

\item{seed}{optional, the seed to be used to split the dataframe and to perform cross validation. Useful to make reproducible the results. }

\item{adaptive}{experimental, if TRUE the adaptive LASSO is implemented.}

\item{b0}{optional, starting values for the regression coefficients. If NULL, the point estimates from \code{glmnet} are used.}

\item{V0}{optional, starting value for the estimates covariance matrix, If NULL, the identity matrix is used.}

\item{c}{the weight of the mixture in the induced smoothed lasso, the default is \code{c = 0.5}. \code{c = -1} means to compute it at each step of the iterative algorithm. }

}
\description{
Auxiliary function for controlling the islasso model fitting.
}
