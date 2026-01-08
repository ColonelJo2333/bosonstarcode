#include "utils/io_commons.hpp"
#include <iostream>
#include <cmath>

using namespace Kadath;
using Io::SolutionKind;

static void usage(const char* prog) {
	std::cerr << "Usage: " << prog << " <solution.dat>\n";
}

int main(int argc, char** argv) {
	if (argc < 2) {
		usage(argv[0]);
		return 1;
	}

	const char* input = argv[1];

	try {
		int kk_guess = 0;
		SolutionKind kind = Io::detect_kind(input, kk_guess);

		if (kind == SolutionKind::Axisymmetric) {
			Io::load_axisymmetric(input, 0.0, [&](const Space_polar& space,
									  int kk,
									  double omega,
									  double lambda,
									  const Scalar& nu,
									  const Scalar& incA,
									  const Scalar& incB,
									  const Scalar& incbt,
									  const Scalar& phi,
									  bool /*has_lambda*/) {

			int ndom = space.get_nbr_domains();

			Scalar lapse(exp(nu));
			lapse.std_base();

			Scalar bigA(exp(incA - nu));
			bigA.std_base();

			Scalar bigB((incB.div_rsint() + 1) / lapse);
			bigB.std_base();

			int nbr = 200;
			double rmax = 100;
			double error_axe = 0;
			Point MM(2);
			for (int i = 0; i < nbr; i++) {
				double current = fabs(bigA.val_point(MM) - bigB.val_point(MM));
				if (current > error_axe)
					error_axe = current;
				MM.set(2) += rmax / nbr;
			}

			Scalar bt(incbt.div_rsint());

			Val_domain integadm(bigA(ndom - 1).der_r());
			double Madm = -space.get_domain(ndom - 1)->integ(integadm, OUTER_BC) / 4 / M_PI;

			Val_domain integkomar(lapse(ndom - 1).der_r());
			double Mkomar = space.get_domain(ndom - 1)->integ(integkomar, OUTER_BC) / 4 / M_PI;

			Val_domain auxiJs(bt(ndom - 1).der_r_rtwo());
			Val_domain integJs(auxiJs.mult_sin_theta().mult_sin_theta());
			double Js = -space.get_domain(ndom - 1)->integ(integJs, OUTER_BC) / 16 / M_PI;

			Scalar auxiJv(kk * (omega - bt * kk) * phi * phi / lapse * bigA * bigA * bigB);
			Scalar integJv(space);
			for (int d = 0; d < ndom; d++)
				integJv.set_domain(d) = auxiJv(d).mult_r().mult_sin_theta();
			double Jv = 0;
			for (int d = 0; d < ndom; d++)
				Jv += space.get_domain(d)->integ_volume(integJv(d));
			Jv *= 2 * M_PI;


			std::cout << "Axisymmetric boson star\n";
			std::cout << "k        = " << kk << "\n";
			std::cout << "omega    = " << omega << "\n";
			std::cout << "lambda   = " << lambda << "\n";
			std::cout << "Madm     = " << Madm << "\n";
			std::cout << "Mkomar   = " << Mkomar << "\n";
			std::cout << "Js       = " << Js << "\n";
			std::cout << "Jv       = " << Jv << "\n";
			std::cout << "diff Komar ADM = " << fabs(Madm - Mkomar) / fabs(Madm + Mkomar) << "\n";
			std::cout << "diff Js and Jv = " << fabs(Js - Jv) / fabs(Js + Jv) << "\n";
		});
		} else {
			Io::load_spherical(input, 0.0, [&](const Space_polar& space,
							   double omega,
							   double lambda,
							   const Scalar& psi,
							   const Scalar& nu,
							   const Scalar& phi,
							   bool has_lambda) {
				(void)space;
				(void)psi;
				std::cout << "Spherical solution\n";
				std::cout << "omega    = " << omega << "\n";
				std::cout << "lambda   = " << lambda << (has_lambda ? " (file)" : " (assumed)") << "\n";
				Point P(2);
				P.set(1) = 0;
				P.set(2) = 0;
				std::cout << "phi(0,0) = " << phi.val_point(P) << "\n";
			});
		}

	} catch (const std::exception& e) {
		std::cerr << "reader failed: " << e.what() << "\n";
		return 1;
	}

	return 0;
}




