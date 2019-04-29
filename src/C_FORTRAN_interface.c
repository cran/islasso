/* $Id: C_FORTRAN_interface.c 313 2015-09-16 20:20:04Z mmaechler $
*
*  wrapper for calling R's random number generator from
*  the original FORTRAN code
*
*/

#include "islasso.h"

double F77_SUB(pnm)(double *x){ return pnorm(*x, 0.0, 1.0, 1, 0); }
double F77_SUB(dnm)(double *x){ return dnorm(*x, 0.0, 1.0, 0); }
double F77_SUB(qnm)(double *x){ return qnorm(*x, 0.0, 1.0, 1, 0); }

double F77_SUB(pnm2)(double *x){ return pnorm(*x, 0.0, 0.001, 1, 0); }
double F77_SUB(dnm2)(double *x){ return dnorm(*x, 0.0, 0.001, 0); }
