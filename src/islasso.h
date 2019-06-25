#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/RS.h>

void F77_SUB(inv)(int *n, double *A, double *invA, int *info);
void F77_SUB(inv2)(int *n, double *A, double *invA, int *info);
void F77_SUB(inv3)(int *n, double *A, double *invA, int *info);

void F77_SUB(solve)(double *A, double *b, int *n, int *info);
void F77_SUB(solve2)(double *A, double *b, int *n, int *info);
void F77_SUB(solve3)(double *A, double *b, int *n, int *info);

void F77_SUB(identitylink)(double *x, int *n, double *mu);
void F77_SUB(identitylinkinv)(double *x, int *n, double *eta);
void F77_SUB(identitymueta)(double *x, int *n, double *eta);

void F77_SUB(loglink)(double *x, int *n, double *mu);
void F77_SUB(loglinkinv)(double *x, int *n, double *eta);
void F77_SUB(logmueta)(double *x, int *n, double *eta);

void F77_SUB(probitlink)(double *x, int *n, double *mu);
void F77_SUB(probitlinkinv)(double *x, int *n, double *eta);
void F77_SUB(probitmueta)(double *x, int *n, double *eta);

void F77_SUB(logitlink)(double *x, int *n, double *mu);
void F77_SUB(logitlinkinv)(double *x, int *n, double *eta);
void F77_SUB(logitmueta)(double *x, int *n, double *eta);

void F77_SUB(inverselink)(double *x, int *n, double *mu);
void F77_SUB(inverselinkinv)(double *x, int *n, double *eta);
void F77_SUB(inversemueta)(double *x, int *n, double *eta);

void F77_SUB(binomial_variance)(double *x, int *n, double *varmu);
void F77_SUB(poisson_variance)(double *x, int *n, double *varmu);
void F77_SUB(gamma_variance)(double *x, int *n, double *varmu);

void
F77_SUB(gradient)(double *theta, double *se, double *lambda, double *xtw,
            double *res, double *pi, int *n, int *p, double *grad, double *alpha);

void
F77_SUB(hessian)(double *theta, double *se, double *lambda, double *xtx,
                  double *pi, int *p, double *hess, double *alpha);

void
F77_SUB(islasso2)(double *X, double *y, int *n, int *p, double *theta, double *se,
                 double *cov, double *lambda, double *alpha, double *pi, int *estpi,
                 double *h, int *itmax, double *tol, double *sigma2, int *trace,
                 int *adaptive, double *offset, int *conv, int *stand, int *intercept,
                 double *eta, double *mu, double *res, double *dev, double *weights,
                 double *hi, double *edf, double *grad2);

void
F77_SUB(islasso_glm)(double *X, double *y, int *n, int *p, double *theta, double *se,
                  double *cov, double *lambda, double *alpha, double *pi, int *estpi,
                  double *h, int *itmax, double *tol, double *sigma2, int *trace,
                  int *adaptive, double *offset, int *conv, int *stand, int *intercept,
                  double *eta, double *mu, double *dev, double *weights, double *hi,
                  double *edf, int *fam, int *link, double *grad2);

void
F77_SUB(family)(int *fam, int *link, int *func, double *x, int *n, double *y);

void
F77_SUB(standardize)(double *X, double *xm, double *xse, int *n, int *p, int *intercept);

void
F77_SUB(check_out)(double *theta, double *cov, double *xm, double *xse, int *p, int *intercept);

void
F77_SUB(crossp)(double *x, double *xtx, int *n, int *p);

void
F77_SUB(tcrossp)(double *x, double *txx, int *n, int *p);

void
F77_SUB(prod1)(double *x, double *w, double *xtw, double *xtx, int *n, int *p);

void
F77_SUB(prod2)(double *xtx, double *tempMat, double *invH, double *cov1, double *hi, int *p);

void
F77_SUB(linear_predictor)(double *x, double *beta, double *eta, double *offset, int *n, int *p);

void
F77_SUB(setdiff)(int *p, int *iprofile, int *ind_noprofile);

void
F77_SUB(armijo)(double *beta, double *dir, double *dtg, double *f0,
                double *alpha, double *h, double *x,
                double *y, double *w, double *offset, int *n, int *p,
                double *lambda, double *eta, double *res, double *eps);

