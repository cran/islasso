#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include "islasso.h"

static const R_FortranMethodDef FortEntries[] = {
    {"solve", (DL_FUNC) &F77_SUB(solve), 4},
    {"inv", (DL_FUNC) &F77_SUB(inv), 4},
    {"identitylink", (DL_FUNC) &F77_SUB(identitylink), 3},
    {"identitylinkinv", (DL_FUNC) &F77_SUB(identitylinkinv), 3},
    {"identitymueta", (DL_FUNC) &F77_SUB(identitymueta), 3},
    {"loglink", (DL_FUNC) &F77_SUB(loglink), 3},
    {"loglinkinv", (DL_FUNC) &F77_SUB(loglinkinv), 3},
    {"logmueta", (DL_FUNC) &F77_SUB(logmueta), 3},
    {"probitlink", (DL_FUNC) &F77_SUB(probitlink), 3},
    {"probitlinkinv", (DL_FUNC) &F77_SUB(probitlinkinv), 3},
    {"probitmueta", (DL_FUNC) &F77_SUB(probitmueta), 3},
    {"logitlink", (DL_FUNC) &F77_SUB(logitlink), 3},
    {"logitlinkinv", (DL_FUNC) &F77_SUB(logitlinkinv), 3},
    {"logitmueta", (DL_FUNC) &F77_SUB(logitmueta), 3},
    {"inverselink", (DL_FUNC) &F77_SUB(inverselink), 3},
    {"inverselinkinv", (DL_FUNC) &F77_SUB(inverselinkinv), 3},
    {"inversemueta", (DL_FUNC) &F77_SUB(inversemueta), 3},
    {"binomial_variance", (DL_FUNC) &F77_SUB(binomial_variance), 3},
    {"poisson_variance", (DL_FUNC) &F77_SUB(poisson_variance), 3},
    {"gamma_variance", (DL_FUNC) &F77_SUB(gamma_variance), 3},
    {"family", (DL_FUNC) &F77_SUB(family), 6},
    {"gradient", (DL_FUNC) &F77_SUB(gradient), 10},
    {"hessian", (DL_FUNC) &F77_SUB(hessian), 8},
    {"islasso3", (DL_FUNC) &F77_SUB(islasso3), 29},
    {"islasso_glm2", (DL_FUNC) &F77_SUB(islasso_glm2), 30},
    {"standardize", (DL_FUNC) &F77_SUB(standardize), 6},
    {"check_out", (DL_FUNC) &F77_SUB(check_out), 6},
    {"crossp", (DL_FUNC) &F77_SUB(crossp), 4},
    {"tcrossp", (DL_FUNC) &F77_SUB(tcrossp), 4},
    {"prod1", (DL_FUNC) &F77_SUB(prod1), 6},
    {"prod2", (DL_FUNC) &F77_SUB(prod2), 6},
    {"linear_predictor", (DL_FUNC) &F77_SUB(linear_predictor), 6},
    {"setdiff", (DL_FUNC) &F77_SUB(setdiff), 3},
    {NULL, NULL, 0}
};

void attribute_visible R_init_islasso(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
