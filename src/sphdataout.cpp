#include "kadath_polar.hpp"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace Kadath;

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " boson_star_lam0.00_om0.80526\n";
        return 1;
    }

    FILE* fin = fopen(argv[1], "r");
    if (!fin) {
        std::cerr << "Cannot open input file\n";
        return 1;
    }

    double omega;
    double lambda;
    Space_polar space(fin);
    fread_be(&omega, sizeof(double), 1, fin);
    fread_be(&lambda, sizeof(double), 1, fin);

    Scalar psi(space, fin);
    Scalar nu (space, fin);
    Scalar phi(space, fin);

    fclose(fin);

    // 物理量
    Scalar A(exp(psi));  A.std_base();
    Scalar B(exp(nu));   B.std_base();

    // 采样参数
    int Nr = 500;
    double rmax = 50.0;
    double theta = 0.0;   // 关键点：只取一个 theta

    std::ofstream fout("sphresult.txt");
    fout << "# r A B phi\n";

    Point P(2);

    for (int ir = 0; ir <= Nr; ++ir) {
        double r = rmax * ir / Nr;

        double x = r * sin(theta);
        double z = r * cos(theta);

        P.set(1) = x;
        P.set(2) = z;

        fout << r << " "
             << A.val_point(P)   << " "
             << B.val_point(P)   << " "
             << phi.val_point(P) << "\n";
    }

    fout.close();
    std::cerr << "Wrote sphresult.txt (theta = 0 only)\n";
    return 0;
}
