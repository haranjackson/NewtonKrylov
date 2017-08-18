#ifndef LGMRES_H
#define LGMRES_H

#include "types.h"
#include "newton_krylov.h"


Vec lgmres(VecFunc matvec, VecFunc psolve,
           Vecr b, Vec x, std::vector<Vec> & outer_v,
           const double tol, const int maxiter,
           const int inner_m, const unsigned int outer_k);


#endif // LGMRES_H
