#include "kadath_polar.hpp"
#include "mpi.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>


using namespace Kadath;
using namespace std;

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // =================================================
    // 1. 空间设置：三维球对称 (r,theta)
    // =================================================

    const int dim = 2;
    int resol = 11;

    Dim_array res(dim);
    res.set(0) = resol;   // r
    res.set(1) = 1;       // 球对称

    Point center(dim);
    center.set(1) = 0;
    center.set(2) = 0;

    int ndom = 7;
    Array<double> bounds(ndom-1);
    bounds.set(0) = 1.0;
    bounds.set(1) = 2.0;
    bounds.set(2) = 4.0;
    bounds.set(3) = 8.0;
    bounds.set(4) = 16.0;
    bounds.set(5) = 48.0;

    Space_polar space(CHEB_TYPE, center, res, bounds);


    double omega  = 0.98;  
    double lambda = 0.0; 
    double phi_c  = 0.005; 

    double pi = M_PI;


    Scalar psi(space); 
    psi.annule_hard();
    psi.std_base();

    Scalar nu(space); 
    nu.annule_hard();
    nu.std_base();

    Scalar phi(space);
    phi.annule_hard();
    phi.std_base();


    Scalar r(space);
    r = 0;
    Scalar one(space);
    one = 1;
    one.std_base();

    for (int d = 0; d < ndom-1; ++d)
        r.set_domain(d) = space.get_domain(d)->mult_r(one(d));
    r.set_domain(ndom-1) = 1;

    for (int d = 0; d < ndom-1; ++d)
        phi.set_domain(d) = phi_c * exp(-r(d)*r(d));
    phi.set_domain(ndom-1) = 0;
    phi.std_base();

    System_of_eqs syst(space, 0, ndom-1);

    syst.add_cst("pi", pi);
    syst.add_cst("lambda", lambda);
    syst.add_cst("phi_c", phi_c);

    syst.add_var("psi", psi);
    syst.add_var("nu",  nu);
    syst.add_var("phi", phi);
    syst.add_var("ome", omega);

    // -------------------------------------------------
    // 定义辅助量
    // -------------------------------------------------

    syst.add_def("Psi = exp(psi)");
    syst.add_def("N   = exp(nu)");

    syst.add_def("V  = phi^2 + 0.5*lambda*phi^4");
    syst.add_def("dV = 1 + lambda*phi^2");   // dV / d|Phi|^2

    // =================================================
    // 5. PDE（严格对应你的三条方程）
    // =================================================

    // ---- (26) Psi 方程 ----
    syst.add_def(
        "eqPsi = lap2(psi) + scal(grad(psi),grad(psi))"
        " + pi*Psi^4*("
        "   (ome*phi/N)^2"
        " + scal(grad(phi),grad(phi))/Psi^4"
        " + V"
        " )"
    );

    // ---- (27) N 方程 ----
    syst.add_def(
        "eqN = lap2(nu) + scal(grad(nu),grad(nu))"
        " + 2*scal(grad(nu),grad(psi))"
        " - 4*pi*Psi^4*( 2*(ome^2)*phi^2/N^2 - V )"
    );

    // ---- (28) Klein–Gordon 方程 ----
    syst.add_def(
        "eqphi = lap2(phi)"
        " - Psi^4*( dV - ome^2/N^2 )*phi"
        " + scal(grad(phi),grad(nu))"
        " + 2*scal(grad(phi),grad(psi))"
    );

    // =================================================
    // 6. 注册为椭圆方程
    // =================================================

    space.add_eq(syst, "eqPsi = 0", "psi", "dn(psi)");
    space.add_eq(syst, "eqN   = 0", "nu",  "dn(nu)");
    space.add_eq(syst, "eqphi= 0", "phi", "dn(phi)");

    // =================================================
    // 7. 中心约束（绑定 omega）
    // =================================================

    Point C(2);
    C.set(1) = 0;

    space.add_eq_point(syst, C, "phi - 0.01");


    // =================================================
    // 8. 外边界条件（渐近平直 + 衰减）
    // =================================================

    syst.add_eq_bc(ndom-1, OUTER_BC, "psi = 0");
    syst.add_eq_bc(ndom-1, OUTER_BC, "nu  = 0");
    syst.add_eq_bc(ndom-1, OUTER_BC, "phi = 0");

    // =================================================
    // 9. Newton 求解
    // =================================================

    bool end = false;
    double conv = 1.0;
    int it = 0;

    while (!end) {
        end = syst.do_newton(1e-9, conv);
        if (rank == 0)
            cout << "Iter " << it
                 << "   conv=" << conv
                 << "   omega=" << omega << endl;
        it++;
        if (it > 30) break;
    }

    // =================================================
    // 10. 保存结果
    // =================================================

    if (rank == 0) 
    {
    std::ostringstream fname;
    fname << "boson_star_"
          << "lam" << std::fixed << std::setprecision(2) << lambda
          << "_om"  << std::fixed << std::setprecision(5) << omega
          << ".dat";

    FILE* f = fopen(fname.str().c_str(), "w");

    space.save(f);
    fwrite_be(&omega,  sizeof(double), 1, f);
    fwrite_be(&lambda, sizeof(double), 1, f);
    psi.save(f);
    nu.save(f);
    phi.save(f);
    fclose(f);

    std::cerr << "Saved solution to " << fname.str() << std::endl;
    }


    MPI_Finalize();
    return 0;
}
