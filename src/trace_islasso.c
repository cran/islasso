#include <R.h>

void F77_SUB(islasso_trace1_7)(double *eps, int *i, double *lmb, double *f0, double *dev, double *redf, double *s2, double *h, double *h2, double *ind3) {
    if(*i == 1)Rprintf("\nIS-lasso algorithm (choosen lambda = %7.3f, threshold = %g)\n\n", *lmb, *eps);
    Rprintf("Step = %4d, h = %1.4f, h2 = %1.4f, DEV = %10.4f (%5.2f df), phi = %7.4f, ||grad||_2 = %12.6f (relative = %2.8f)\n", *i, *h, *h2, *dev, *redf, *s2, *f0, *ind3);
}

void F77_SUB(islasso_trace1_7_2)(double *eps, int *i, double *lmb, double *dev, double *redf, double *s2, double *ind, double *ind2) {
    if(*i == 1)Rprintf("\nIS-lasso algorithm (choosen lambda = %7.3f, threshold = %g)\n\n", *lmb, *eps);
    Rprintf("Step = %4d, DEV = %10.4f (%5.2f df), phi = %7.4f, ||(BETAn - BETAo)||_1 = %2.8f e ||(SEn - SEo)||_1 = %2.8f\n", *i, *dev, *redf, *s2, *ind2, *ind);
}

void F77_SUB(islasso_trace1_2)(double *eps, int *i, double *lmb, double *f0, double *dev, double *redf, double *s2, double *h, double *h2, double *ind, double *ind2, double *ind3, double *ind4) {
    Rprintf("\n=============================================================\n");
    Rprintf("IS-lasso algorithm step = %d\n", *i);
    Rprintf("  Choosen lambda value = %7.3f\n", *lmb);
    Rprintf("  Step halving (beta), h = %1.6f\n", *h);
    Rprintf("  Step halving (cov), h2 = %1.6f\n", *h2);
    Rprintf("  Residual deviance = %10.6f on %5.2f degrees of freedom\n", *dev, *redf);
    Rprintf("  Estimated dispersion parameter = %7.4f\n", *s2);
    Rprintf("  Checking convergence criterion (threshold = %g):\n", *eps);
    Rprintf("     ||(SEn - SEo)||_2 = %2.8f\n", *ind);
    Rprintf("     ||(BETAn - BETAo)||_2 = %2.8f\n", *ind2);
    Rprintf("     (||d(SE)||_2 + ||d(BETA)||_2)/2 = %2.8f\n", *ind4);
    Rprintf("     ||gradn||_2 = %10.6f (relative = %2.8f)", *f0, *ind3);
}

void F77_SUB(islasso_trace1_2_2)(double *eps, int *i, double *lmb, double *dev, double *redf, double *s2, double *ind, double *ind2) {
    Rprintf("\n=============================================================\n");
    Rprintf("IS-lasso algorithm step = %d\n", *i);
    Rprintf("  Choosen lambda value = %7.3f\n", *lmb);
    Rprintf("  Residual deviance = %10.6f on %5.2f degrees of freedom\n", *dev, *redf);
    Rprintf("  Estimated dispersion parameter = %7.4f\n", *s2);
    Rprintf("  Checking convergence criterion (threshold = %g):\n", *eps);
    Rprintf("     ||(BETAn - BETAo)||_1 = %2.8f\n", *ind2);
    Rprintf("     ||(SEn - SEo)||_1 = %2.8f\n", *ind);
}

void F77_SUB(islasso_trace1_8)() {
    Rprintf("\n===================================");
    Rprintf("\nConvergence criterion is met!\n\n");
}

void F77_SUB(islasso_trace2_5)(double *eps, double *lmb) {
    Rprintf("\nIS-lasso (GLM) algorithm (choosen lambda = %7.3f, threshold = %g)\n\n", *lmb, *eps);
}

void F77_SUB(islasso_trace2_7)(int *i, int *i2, double *dev, double *redf, double *s2, double *ind3, double *ind4) {
    Rprintf("Step = %4d, nit = %4d, deviance = %10.4f (%5.2f df), phi = %7.4f, (||dSE|| + ||dBETA||)/2 = %2.8f (relative = %2.8f)\n", *i, *i2, *dev, *redf, *s2, *ind4, *ind3);
}

void F77_SUB(islasso_trace2_7_2)(double *eps, int *i, double *lmb, double *dev, double *redf, double *s2, double *ind, double *ind2) {
    if(*i == 1)Rprintf("\nIS-lasso (GLM) algorithm (choosen lambda = %7.3f, threshold = %g)\n\n", *lmb, *eps);
    Rprintf("Step = %4d, DEV = %10.4f (%5.2f df), phi = %7.4f, ||(BETAn - BETAo)||_1 = %2.8f e ||(SEn - SEo)||_1 = %2.8f\n", *i, *dev, *redf, *s2, *ind2, *ind);
}

void F77_SUB(islasso_trace2_2)(double *eps, int *i, int *i2, double *lmb, double *f0, double *dev, double *redf, double *s2, double *ind, double *ind2, double *ind3, double *ind4) {
    Rprintf("\n=============================================================\n");
    Rprintf("IS-lasso (GLM) algorithm step = %d\n", *i);
    Rprintf("  Fisher scoring iteration = %d\n", *i2);
    Rprintf("  Choosen lambda value = %7.3f\n", *lmb);
    Rprintf("  Residual deviance = %10.6f on %5.2f degrees of freedom\n", *dev, *redf);
    Rprintf("  Estimated dispersion parameter = %7.4f\n", *s2);
    Rprintf("  Checking convergence criterion (threshold = %g):\n", *eps);
    Rprintf("     ||gradn||_2 = %10.6f\n", *f0);
    Rprintf("     ||(SEn - SEo)||_2 = %2.8f\n", *ind);
    Rprintf("     ||(BETAn - BETAo)||_2 = %2.8f\n", *ind2);
    Rprintf("     (||d(SE)||_2 + ||d(BETA)||_2)/2 = %2.8f (relative = %2.8f)\n", *ind4, *ind3);
}

void F77_SUB(islasso_trace2_2_2)(double *eps, int *i, double *lmb, double *dev, double *redf, double *s2, double *ind, double *ind2) {
    Rprintf("\n=============================================================\n");
    Rprintf("IS-lasso (GLM) algorithm step = %d\n", *i);
    Rprintf("  Choosen lambda value = %7.3f\n", *lmb);
    Rprintf("  Residual deviance = %10.6f on %5.2f degrees of freedom\n", *dev, *redf);
    Rprintf("  Estimated dispersion parameter = %7.4f\n", *s2);
    Rprintf("  Checking convergence criterion (threshold = %g):\n", *eps);
    Rprintf("     ||(BETAn - BETAo)||_1 = %2.8f\n", *ind2);
    Rprintf("     ||(SEn - SEo)||_1 = %2.8f\n", *ind);
}

void F77_SUB(islasso_trace2_3)(int *i) {
    if(*i == 1) Rprintf("\nFisher scoring step:\n");
    Rprintf(".");
    if((*i % 100) == 0) Rprintf(" %5d\n", *i);
}
void F77_SUB(islasso_trace2_6)(int *i) {
    Rprintf(" %5d. Convergence criterion is met!\n\n", *i);
}




void F77_SUB(islasso_trace1_3_2)(int *conv) {
    Rprintf("\n\tinversion problem (Inversion function) = %d\n", *conv);
}

void F77_SUB(islasso_trace1_4_2)(int *conv) {
    Rprintf("\n\tsolve problem (Solve function) = %d\n", *conv);
}




void F77_SUB(islasso_trace2_4)(int *i, double *h, double *dev, double *devold) {
    Rprintf("\n\nIter %d, h = %5.3f, dev = %5.5f (old = %5.5f)\n", *i, *h, *dev, *devold);
}

void F77_SUB(islasso_trace1_3)(int *conv) {
    Rprintf("\n\tinversion problem = %d\n", *conv);
}

void F77_SUB(islasso_trace1_4)(int *conv) {
    Rprintf("\n\tsolve problem = %d\n", *conv);
}


