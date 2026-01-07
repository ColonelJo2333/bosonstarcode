// File: plot_initial.cpp
//
// Compiles with: g++ plot_initial.cpp -o plot_initial $(kadath_gflags)
// Runs with: ./plot_initial
//
// This program only computes the initial Gaussian wave packet for phi
// and saves the (z, x, phi) values to a text file "initial_phi_data.txt".

#include "kadath_polar.hpp"
#include "mpi.h"
#include <fstream> // For file output
#include <iostream> // For console output

using namespace Kadath;

int main(int argc, char** argv) {

    int rc = MPI_Init(&argc, &argv);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // --- Space and Parameter Setup (Identical to original code) ---
    int dim = 2;
    int type_coloc = CHEB_TYPE;
    int resol = 13; // Increased resolution for a better plot
    Dim_array res(dim);
    res.set(0) = resol; res.set(1) = resol -1 ;

    Point center(2);
    for (int i = 1; i <= dim; i++)
        center.set(i) = 0;
    int ndom = 8;
    Array<double> bounds(ndom - 1);
    bounds.set(0) = 1; bounds.set(1) = 2; bounds.set(2) = 4; bounds.set(3) = 8;
    bounds.set(4) = 16; bounds.set(5) = 32; bounds.set(6) = 48;
    Space_polar space(type_coloc, center, res, bounds);

    int kk = 1;
    Scalar one(space);
    one = 1;
    one.std_base();
    Scalar rsint(space);
    for (int d = 0; d < ndom - 1; d++)
        rsint.set_domain(d) = space.get_domain(d)->mult_r(space.get_domain(d)->mult_sin_theta(one(d)));
    rsint.set_domain(ndom - 1) = space.get_domain(ndom - 1)->mult_sin_theta(one(ndom - 1));

    double posmax = 2.5;
    double fmax = 0.05;
    double sigmax = 2 * posmax * posmax / double(kk);
    double sigmaz = sigmax / 4;
    double vmax = fmax / exp(0) / pow(posmax, kk) / exp(-posmax * posmax / sigmax);

    // --- Compute the initial phi field ---
    Scalar phi(space);
    for (int d = 0; d < ndom - 1; d++)
        phi.set_domain(d) = vmax * pow(rsint(d), kk) * exp(-pow(space.get_domain(d)->get_cart(1), 2) / sigmax) * exp(-pow(space.get_domain(d)->get_cart(2), 2) / sigmaz);
    phi.set_domain(ndom - 1) = 0; // Set to zero in the outer domain (compact support)
    phi.std_anti_base();


    // --- NEW PART: EXPORT DATA TO A FILE ---
    if (rank == 0) {
        std::cout << "Exporting initial 'phi' field data..." << std::endl;
        std::ofstream outfile("initial_phi_data.txt");
        if (!outfile.is_open()) {
            std::cerr << "Error: Could not open output file." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        FILE* fiche = fopen ("init.dat", "w") ;
        space.save(fiche) ;
        phi.save(fiche) ;
        fclose (fiche) ;
    }

    

    MPI_Finalize();
    return EXIT_SUCCESS;
}
