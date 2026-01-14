#include "utils/io_commons.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

using namespace Kadath;
using Io::SolutionKind;

static void usage(const char* prog) {
	std::cerr << "Usage: " << prog << " <solution.dat> <tar> [output.txt]\n";
	std::cerr << "  tar=0  compute mode (auto-detect spherical/axisymmetric)\n";
	std::cerr << "  tar=1  export mode (auto-detect spherical/axisymmetric)\n";
	std::cerr << "  output.txt defaults to <solution>.txt when tar=1\n";
}

static std::string default_output_path(const std::string& input) {
	if (input.size() >= 4 && input.substr(input.size() - 4) == ".dat")
		return input.substr(0, input.size() - 4) + ".txt";
	return input + ".txt";
}

static void write_txt(const std::string& path,
					  const Space_polar& space,
					  double omega,
					  double lambda,
					  const std::vector<std::string>& names,
					  const std::vector<const Scalar*>& fields) {
	if (names.size() != fields.size())
		throw std::runtime_error("column names and fields size mismatch");

	std::ofstream out(path);
	if (!out)
		throw std::runtime_error("Cannot open output file: " + path);

	out << "# omega=" << omega << " lambda=" << lambda << "\n";
	out << "# ";
	for (size_t i = 0; i < names.size(); ++i) {
		out << names[i];
		if (i + 1 < names.size())
			out << " ";
	}
	out << "\n";
	out << std::setprecision(16);

	int ndom = space.get_nbr_domains();
	for (int d = 0; d < ndom; ++d) {
		const Domain* dom = space.get_domain(d);
		Dim_array dims = dom->get_nbr_points();
		Index idx(dims);
		for (const Scalar* field : fields)
			field->operator()(d).coef_i();
		bool more = true;
		while (more) {
			for (size_t i = 0; i < fields.size(); ++i) {
				double val = (*fields[i])(d)(idx);
				out << val;
				if (i + 1 < fields.size())
					out << " ";
			}
			out << "\n";
			more = idx.inc();
		}
	}
}

int main(int argc, char** argv) {
	if (argc < 3) {
		usage(argv[0]);
		return 1;
	}

	const char* input = argv[1];
	int tar = std::atoi(argv[2]);
	if (tar != 0 && tar != 1) {
		std::cerr << "Invalid tar: " << argv[2] << "\n";
		usage(argv[0]);
		return 1;
	}
	std::string output_path;
	if (tar == 1) {
		if (argc >= 4)
			output_path = argv[3];
		else
			output_path = default_output_path(input);
	}

	try {
		int kk_guess = 0;
		SolutionKind kind = Io::detect_kind(input, kk_guess);
		std::cout << "mode      = " << (tar == 0 ? "compute" : "export") << "\n";
		std::cout << "data_type = " << (kind == SolutionKind::Axisymmetric ? "axisymmetric" : "spherical") << "\n";

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

			Scalar ap(exp(nu));
			ap.std_base();

			Scalar A(exp(incA - nu));
			A.std_base();

			Scalar B((incB.div_rsint() + 1) / ap);
			B.std_base();

			Scalar bt(incbt.div_rsint());
			bt.std_base();

			if (tar == 0) {
				int ndom = space.get_nbr_domains();

				Val_domain integadm(A(ndom - 1).der_r());
				double Madm = -space.get_domain(ndom - 1)->integ(integadm, OUTER_BC) / 4 / M_PI;

				Val_domain integkomar(ap(ndom - 1).der_r());
				double Mkomar = space.get_domain(ndom - 1)->integ(integkomar, OUTER_BC) / 4 / M_PI;

				Val_domain auxiJs(bt(ndom - 1).der_r_rtwo());
				Val_domain integJs(auxiJs.mult_sin_theta().mult_sin_theta());
				double Js = -space.get_domain(ndom - 1)->integ(integJs, OUTER_BC) / 16 / M_PI;

				Scalar auxiJv(kk * (omega - bt * kk) * phi * phi / ap * A * A * B);
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
			} else {
				std::vector<std::string> names = {"ap", "A", "B", "bt", "phi"};
				std::vector<const Scalar*> fields = {&ap, &A, &B, &bt, &phi};
				write_txt(output_path, space, omega, lambda, names, fields);
				std::cout << "exported: " << output_path << "\n";
			}
		});
		} else {
			Io::load_spherical(input, 0.0, [&](const Space_polar& space,
							   double omega,
							   double lambda,
							   const Scalar& psi,
							   const Scalar& nu,
							   const Scalar& phi,
							   bool has_lambda) {
				Scalar Psi(exp(psi));
				Psi.std_base();

				Scalar N(exp(nu));
				N.std_base();

				if (tar == 0) {
					int ndom = space.get_nbr_domains();
					Scalar A(Psi * Psi * Psi *Psi);
					A.std_base();

					Val_domain integadm(A(ndom - 1).der_r());
					double Madm = -space.get_domain(ndom - 1)->integ(integadm, OUTER_BC) / 4 / M_PI;

					Val_domain integkomar(N(ndom - 1).der_r());
					double Mkomar = space.get_domain(ndom - 1)->integ(integkomar, OUTER_BC) / 4 / M_PI;

					std::cout << "Spherical solution\n";
					std::cout << "omega    = " << omega << "\n";
					std::cout << "lambda   = " << lambda << (has_lambda ? " (file)" : " (assumed)") << "\n";
					std::cout << "Madm     = " << Madm << "\n";
					std::cout << "Mkomar   = " << Mkomar << "\n";
					std::cout << "diff Komar ADM = " << fabs(Madm - Mkomar) / fabs(Madm + Mkomar) << "\n";
				} else {
					std::vector<std::string> names = {"Psi", "N", "phi"};
					std::vector<const Scalar*> fields = {&Psi, &N, &phi};
					write_txt(output_path, space, omega, lambda, names, fields);
					std::cout << "exported: " << output_path << "\n";
				}
			});
		}

	} catch (const std::exception& e) {
		std::cerr << "reader failed: " << e.what() << "\n";
		return 1;
	}

	return 0;
}


