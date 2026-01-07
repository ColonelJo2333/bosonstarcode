#include "kadath_polar.hpp"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace Kadath ;

int main(int argc, char** argv) {

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << "bosl_1_0.798905_-1.000000\n";
        return 1;
    }

    const char* input_name = argv[1];
    FILE* fin = fopen(input_name, "r");
    if (!fin) {
        std::cerr << "Cannot open input file " << input_name << std::endl;
        return 1;
    }

    int kk;
    double omega;
    double lambda;

    Space_polar space(fin);
    fread_be(&kk, sizeof(int), 1, fin);
    fread_be(&omega, sizeof(double), 1, fin);
    fread_be(&lambda, sizeof(double), 1, fin);

    Scalar nu   (space, fin);
    Scalar incA (space, fin);
    Scalar incB (space, fin);
    Scalar incbt(space, fin);
    Scalar phi  (space, fin);

    fclose(fin);

    Param_tensor parameters;
    parameters.set_m_quant() = kk;
    phi.set_parameters() = parameters;

    Scalar ap (exp(nu));
    ap.std_base();

    Scalar bt (incbt.div_rsint());
    bt.std_base();

    Scalar bigB ((incB.div_rsint() + 1.0) / ap);
    bigB.std_base();

    Scalar bigA (exp(incA - nu));
    bigA.std_base();

    std::cerr << "Constructed metric quantities (ap, bigA, bigB, bt)." << std::endl;


    int Nr  = 500;
    int Nth = 500; 
    double rmax = 15.0;

    const char* output_name = "metricfields_2d.txt";
    std::ofstream fout(output_name);
    if (!fout) {
        std::cerr << "Cannot open output file " << output_name << std::endl;
        return 1;
    }


    fout << "# r theta ap A B bt phi\n";

    Point P(2);
    for (int it = 0; it <= Nth; ++it) {
        double theta = M_PI * it / Nth;
        for (int ir = 0; ir <= Nr; ++ir) {
            double r = rmax * ir / Nr;

            double x = r * sin(theta);
            double z = r * cos(theta);

            P.set(1) = x;
            P.set(2) = z;

            double ap_val  = ap.val_point(P);
            double A_val   = bigA.val_point(P);
            double B_val   = bigB.val_point(P);
            double bt_val  = bt.val_point(P);
            double phi_val = phi.val_point(P);

            fout << r << " " << theta << " "
                 << ap_val << " "
                 << A_val  << " "
                 << B_val  << " "
                 << bt_val << " "
                 << phi_val << "\n";
        }
        fout << "\n";
    }

    fout.close();
    std::cerr << "Wrote metricfields_2d.txt\n";
    return 0;
}