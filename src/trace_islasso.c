#include <R.h>

void F77_SUB(islasso_trace1_1)(void) {
    Rprintf("\nGaussian iteration\n");
}

void F77_SUB(islasso_trace2_1)(void) {
    Rprintf("\nGLM iteration\n");
}

void F77_SUB(islasso_trace2_2)(double *varmu, double *mu_eta_val, double *res, double *z, double *w) {
    Rprintf("\nvarmu = %f\n", *varmu);
    Rprintf("\nmu_eta_val = %f\n", *mu_eta_val);
    Rprintf("\nres = %f\n", *res);
    Rprintf("\nz = %f\n", *z);
    Rprintf("\nw = %f\n", *w);
}

void F77_SUB(islasso_trace2_3)(double *theta, double *mu, double *eta) {
    Rprintf("\ntheta = %f\n", *theta);
    Rprintf("\nmu = %f\n", *mu);
    Rprintf("\neta = %f\n", *eta);
}

void F77_SUB(islasso_trace1_2)(double *eps, int *i, double *ind) {
    Rprintf("\n\tIS-lasso algorithm step = %d\n", *i);
    Rprintf("\n\tChecking convergence criterion (threshold = %f)\n", *eps);
    Rprintf("\tmax||(SEn - SEo)||_1 = %f\n", *ind);
    if(*ind <= *eps) Rprintf("\n\tConvergence criterion is met!\n");
}

void F77_SUB(islasso_trace1_3)(int *conv) {
    Rprintf("\n\tinversion problem = %d\n", *conv);
}

void F77_SUB(islasso_trace1_4)(int *conv) {
    Rprintf("\n\tsolve problem = %d\n", *conv);
}

void F77_SUB(islasso_trace1_3_2)(int *conv) {
    Rprintf("\n\tinversion problem (Inversion function) = %d\n", *conv);
}

void F77_SUB(islasso_trace1_4_2)(int *conv) {
    Rprintf("\n\tsolve problem (Solve function) = %d\n", *conv);
}


void F77_SUB(islasso_trace1_5)(double *theta, double *se, double *lambda, double *xtx, double *pi, int *p, double *hess) {
    Rprintf("\n\thessian arguments\n");
    Rprintf("\ntheta = %f\n", *theta);
    Rprintf("\nse = %f\n", *se);
    Rprintf("\nlambda = %f\n", *lambda);
    Rprintf("\nxtx = %f\n", *xtx);
    Rprintf("\npi = %f\n", *pi);
    Rprintf("\np = %d\n", *p);
    Rprintf("\nhess = %f\n", *hess);
}
void F77_SUB(islasso_trace1_6)(double *theta, double *se, double *lambda, double *xtw, double *res, double *pi, int *n, int *p, double *grad) {
    Rprintf("\n\tgradient arguments\n");
    Rprintf("\ntheta = %f\n", *theta);
    Rprintf("\nse = %f\n", *se);
    Rprintf("\nlambda = %f\n", *lambda);
    Rprintf("\nxtw = %f\n", *xtw);
    Rprintf("\nres = %f\n", *res);
    Rprintf("\npi = %f\n", *pi);
    Rprintf("\nn = %d\n", *n);
    Rprintf("\np = %d\n", *p);
    Rprintf("\ngrad = %f\n", *grad);
}
