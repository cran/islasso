\name{islasso.control}
\alias{islasso.control}
\title{Auxiliary for controlling islasso model fitting}
\usage{
islasso.control(sigma2 = -1, h = 0.15, tol = 1e-05, itmax = 500, 
  display = FALSE, use.starts = TRUE, debias = FALSE, standardize = TRUE,
  adaptive = FALSE, b0 = NULL, V0 = NULL, pai = 1, estpai = TRUE)
}
\arguments{
\item{sigma2}{if -1 (default) sigma is estimated from the data}

\item{h}{step halving, default is 0.15}

\item{tol}{tollerance for convergence, dafault is 1e-5}

\item{itmax}{maximum number of iterations, default 500}

\item{display}{if TRUE, shows each iteration}

\item{use.starts}{if TRUE, sets up glmnet coefficients as starting points}

\item{debias}{if TRUE, allows to correct estimates}

\item{standardize}{if TRUE, for X variable standardization prior to fitting the model, the coefficients are always returnd on the original scale}

\item{adaptive}{if TRUE, for adaptive LASSO}

\item{b0}{is for setting up starting coefficients}

\item{V0}{is for setting up starting variance and covariace matrix}

\item{pai}{the weight of the mixture in the induced smoothed lasso}

\item{estpai}{if TRUE, the probability of the mixtures are estimated}

\item{...}{further arguments passed to or from other methods.}
}
\description{
Auxiliary for controlling islasso model fitting
}
