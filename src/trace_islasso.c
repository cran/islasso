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

void F77_SUB(islasso_trace2_4)(int *i, double *h, double *dev, double *devold) {
    Rprintf("\n\nIter %d, h = %5.3f, dev = %5.5f (old = %5.5f)\n", *i, *h, *dev, *devold);
}

void F77_SUB(islasso_trace1_2)(double *eps, int *i, double *lmb, double *dev, double *redf, double *s2, double *h, double *ind, double *ind2) {
    Rprintf("\n==================================================================\n");
    Rprintf("IS-lasso algorithm step = %d\n", *i);
    Rprintf("  Choosen lambda value = %7.3f\n", *lmb);
    Rprintf("  Step halving = %1.8f\n", *h);
    Rprintf("  Residual deviance = %8.4f on %5.2f degrees of freedom\n", *dev, *redf);
    Rprintf("  Estimated dispersion parameter = %7.4f\n", *s2);
    Rprintf("  Checking convergence criterion (threshold = %g):\n", *eps);
    Rprintf("     max|(SEn - SEo)| = %2.8f\n", *ind);
    Rprintf("     max|(BETAn - BETAo)| = %2.8f", *ind2);
}

void F77_SUB(islasso_trace1_7)(double *eps, int *i, double *lmb, double *dev, double *redf, double *s2, double *h, double *ind, double *ind2) {
    if(*i == 1)Rprintf("\nIS-lasso algorithm (choosen lambda = %7.3f, threshold = %g)\n\n", *lmb, *eps);
    Rprintf("Step = %4d, h = %1.8f, deviance = %8.4f (%5.2f df), phi = %7.4f, max|(SEn - SEo)| = %2.8f\n", *i, *h, *dev, *redf, *s2, *ind, *eps);
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

void F77_SUB(islasso_trace1_8)() {  Rprintf("\n================================================================================================================\nConvergence criterion is met!\n");
}
