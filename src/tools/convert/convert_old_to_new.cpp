#include "utils/io_commons.hpp"
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

using namespace Kadath;
using Io::SolutionKind;

int main(int argc, char** argv) {
    (void)argc;
    (void)argv;

    // 配置：在此修改输入/输出文件和可选的 lambda 覆盖值
    const char* input = "bosinit.dat";                // 旧格式输入文件
    std::string output = "tools/converted/bosinit.dat"; // 转换后输出文件
    bool override_lambda = false;                      // 设为 true 可强制覆盖 lambda_override
    double lambda_override = 0.0;                      // 覆盖的 lambda 值

    try {
        int kk_guess = 0;
        SolutionKind kind = Io::detect_kind(input, kk_guess);

        if (kind == SolutionKind::Axisymmetric) {
            Io::load_axisymmetric(input, lambda_override, [&](const Space_polar& space,
                                                             int kk,
                                                             double omega,
                                                             double lambda,
                                                             const Scalar& nu,
                                                             const Scalar& incA,
                                                             const Scalar& incB,
                                                             const Scalar& incbt,
                                                             const Scalar& phi,
                                                             bool has_lambda) {
                double lambda_out = override_lambda ? lambda_override : (has_lambda ? lambda : 0.0);
                Io::save_axisymmetric(output.c_str(), space, kk, omega, lambda_out, nu, incA, incB, incbt, phi);
                std::cerr << "Converted axisymmetric file: kk=" << kk
                          << " omega=" << omega
                          << " lambda=" << lambda_out
                          << " -> " << output << "\n";
            });
        } else {
            Io::load_spherical(input, lambda_override, [&](const Space_polar& space,
                                                           double omega,
                                                           double lambda,
                                                           const Scalar& psi,
                                                           const Scalar& nu,
                                                           const Scalar& phi,
                                                           bool has_lambda) {
                double lambda_out = override_lambda ? lambda_override : (has_lambda ? lambda : 0.0);
                Io::save_spherical(output.c_str(), space, omega, lambda_out, psi, nu, phi);
                std::cerr << "Converted spherical file: omega=" << omega
                          << " lambda=" << lambda_out
                          << " -> " << output << "\n";
            });
        }

    } catch (const std::exception& e) {
        std::cerr << "Conversion failed: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
