#include <cmath>

#include "types.h"
#include "newton_krylov.h"


DecQR qr_append(Matr Q, Matr R, Matr u)
{   // Returns QR factorization of Q*R with u appended as the last column
    Mat A = Q * R;
    A.conservativeResize(A.rows(),A.cols()+1);
    A.col(A.cols()-1) = u;
    DecQR ret;
    return ret.compute(A);
}

Vec lgmres(VecFunc matvec, VecFunc psolve,
           Vecr b, Vec x, std::vector<Vec> & outer_v,
           const double tol=1e-5, const int maxiter=1,
           const int inner_m=30, const unsigned int outer_k=10)
{
    int N = b.size();

    double b_norm = b.norm();
    if (b_norm == 0.)
        b_norm = 1.;

    for (int k_outer=0; k_outer<maxiter; k_outer++)
    {
        Vec r_outer = matvec(x) - b;

        double r_norm = r_outer.norm();
        if (r_norm <= tol * b_norm || r_norm <= tol)
            break;

        Vec vs0 = - psolve(r_outer);
        double inner_res_0 = vs0.norm();

        vs0 /= inner_res_0;
        std::vector<Vec> vs = {vs0};
        std::vector<Vec> ws;

        Mat Q = Mat::Ones(1,1);
        Mat R = Mat::Zero(1,0);

        unsigned int ind = 1 + inner_m + outer_v.size();
        unsigned int j;
        for (j=1; j<ind; j++)
        {
            Vec z(N);
            if (j < outer_v.size() + 1)
                z = outer_v[j-1];
            else if (j == outer_v.size() + 1)
                z = vs0;
            else
                z = vs.back();

            Vec v_new = psolve(matvec(z));
            double v_new_norm = v_new.norm();

            Vec hcur = Vec::Zero(j+1);
            for (unsigned int i=0; i<vs.size(); i++)
            {
                Vec v = vs[i];
                double alpha = v.dot(v_new);
                hcur(i) = alpha;
                v_new = v_new - alpha * v;
            }
            hcur(j) = v_new.norm();

            v_new /= hcur(j);
            vs.push_back(v_new);
            ws.push_back(z);

            Mat Q2 = Mat::Zero(j+1, j+1);
            Q2.topLeftCorner(j,j) = Q;
            Q2(j,j) = 1.;

            Mat R2 = Mat::Zero(j+1, j-1);
            R2.topRows(j) = R;

            DecQR QR = qr_append(Q2, R2, hcur);
            Q.conservativeResize(j+1, j+1);
            R.conservativeResize(j+1, j-1);
            Q = QR.householderQ();
            R = QR.matrixQR().triangularView<Eigen::Upper>();

            double inner_res = std::abs(Q(0,j)) * inner_res_0;

            if ((inner_res <= tol * inner_res_0) ||
                    (hcur(j) <= mEPS * v_new_norm))
            {
                break;
            }
        }
        Vec b = Q.topLeftCorner(1,j).transpose();
        Mat A = R.topLeftCorner(j,j);

        Vec y = A.jacobiSvd(Eigen::ComputeThinU|Eigen::ComputeThinV).solve(b);
        y *= inner_res_0;

        Vec dx = y(0) * ws[0];
        for (int i=1; i<y.size(); i++)
            dx += y(i) * ws[i];

        double nx = dx.norm();
        if (nx > 0.)
            outer_v.push_back(dx/nx);

        while (outer_v.size()>outer_k)
            outer_v.erase(outer_v.begin());

        x += dx;
    }
    return x;
}
