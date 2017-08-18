#ifndef NEWTON_KRYLOV_H
#define NEWTON_KRYLOV_H

#include <limits>

#include "types.h"


const double INF = std::numeric_limits<double>::max();
const double mEPS = 2.2204460492503131e-16;


class KrylovJacobian
{
    VecFunc func;
    int maxiter;
    int inner_m;
    unsigned int outer_k;

    Vec x0;
    Vec f0;
    double rdiff;
    double omega;
    std::vector<Vec> outer_v;
    Mat M;

public:
    void update_diff_step();
    KrylovJacobian(Vecr x, Vecr f, VecFunc F);
    Vec matvec(Vec v);
    Vec psolve(Vec v);
    Vec solve(Vecr rhs, double tol);
    void update(Vecr x, Vecr f);
};


class TerminationCondition
{
    double f_tol;
    double f_rtol;
    double x_tol;
    double x_rtol;
    double f0_norm;
public:
    TerminationCondition(double ftol, double frtol, double xtol, double xrtol);
    int check(Vecr f, Vecr x, Vecr dx);
};


Vec nonlin_solve(VecFunc F, Vecr x,
                 double f_tol, double f_rtol, double x_tol, double x_rtol);


#endif // NEWTON_KRYLOV_H
