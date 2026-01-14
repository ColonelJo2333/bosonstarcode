#pragma once

#include "kadath_polar.hpp"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include <string>
#include <functional>

namespace Io {

static constexpr int KK_MAX_ABS = 32;

enum class SolutionKind {
    Axisymmetric,
    Spherical
};

inline SolutionKind detect_kind(const char* path, int& kk_guess) {
    FILE* f = fopen(path, "r");
    if (!f) throw std::runtime_error(std::string("Cannot open file: ") + path);
    Kadath::Space_polar space(f);
    (void)space;
    long pos = ftell(f);
    int kk = 0;
    size_t r = Kadath::fread_be(&kk, sizeof(int), 1, f);
    fclose(f);
    kk_guess = kk;
    if (r == 1 && std::abs(kk) <= KK_MAX_ABS) {
        return SolutionKind::Axisymmetric;
    }
    return SolutionKind::Spherical;
}

inline bool looks_like_lambda(double val) {
    return std::isfinite(val) && std::fabs(val) < 1e3;
}

template <class Fn>
inline void load_axisymmetric(const char* path, double lambda_default, Fn&& fn) {
    FILE* f = fopen(path, "r");
    if (!f) throw std::runtime_error(std::string("Cannot open file: ") + path);

    Kadath::Space_polar space(f);
    int kk = 0;
    double omega = 0.0;
    double lambda = lambda_default;

    if (Kadath::fread_be(&kk, sizeof(int), 1, f) != 1)
        throw std::runtime_error("Failed to read kk");
    if (Kadath::fread_be(&omega, sizeof(double), 1, f) != 1)
        throw std::runtime_error("Failed to read omega");

    long pos = ftell(f);
    double lambda_probe = lambda_default;
    bool has_lambda = false;
    if (Kadath::fread_be(&lambda_probe, sizeof(double), 1, f) == 1 && looks_like_lambda(lambda_probe)) {
        long after_lambda = ftell(f);
        int base_flag = 0;
        int ndim_flag = 0;
        bool peek_ok = Kadath::fread_be(&base_flag, sizeof(int), 1, f) == 1 &&
                       Kadath::fread_be(&ndim_flag, sizeof(int), 1, f) == 1;
        if (peek_ok && (base_flag == 0 || base_flag == 1) && ndim_flag == space.get_ndim()) {
            lambda = lambda_probe;
            has_lambda = true;
            fseek(f, after_lambda, SEEK_SET);
        } else {
            fseek(f, pos, SEEK_SET);
        }
    } else {
        fseek(f, pos, SEEK_SET);
    }

    Kadath::Scalar nu   (space, f);
    Kadath::Scalar incA (space, f);
    Kadath::Scalar incB (space, f);
    Kadath::Scalar incbt(space, f);
    Kadath::Scalar phi  (space, f);

    fclose(f);

    fn(space, kk, omega, lambda, nu, incA, incB, incbt, phi, has_lambda);
}

template <class Fn>
inline void load_spherical(const char* path, double lambda_default, Fn&& fn) {
    FILE* f = fopen(path, "r");
    if (!f) throw std::runtime_error(std::string("Cannot open file: ") + path);

    Kadath::Space_polar space(f);
    double omega = 0.0;
    double lambda = lambda_default;

    if (Kadath::fread_be(&omega, sizeof(double), 1, f) != 1)
        throw std::runtime_error("Failed to read omega");

    long pos = ftell(f);
    double lambda_probe = lambda_default;
    bool has_lambda = false;
    if (Kadath::fread_be(&lambda_probe, sizeof(double), 1, f) == 1 && looks_like_lambda(lambda_probe)) {
        long after_lambda = ftell(f);
        int base_flag = 0;
        int ndim_flag = 0;
        bool peek_ok = Kadath::fread_be(&base_flag, sizeof(int), 1, f) == 1 &&
                       Kadath::fread_be(&ndim_flag, sizeof(int), 1, f) == 1;
        if (peek_ok && (base_flag == 0 || base_flag == 1) && ndim_flag == space.get_ndim()) {
            lambda = lambda_probe;
            has_lambda = true;
            fseek(f, after_lambda, SEEK_SET);
        } else {
            fseek(f, pos, SEEK_SET);
        }
    } else {
        fseek(f, pos, SEEK_SET);
    }

    Kadath::Scalar psi(space, f);
    Kadath::Scalar nu (space, f);
    Kadath::Scalar phi(space, f);

    fclose(f);
    fn(space, omega, lambda, psi, nu, phi, has_lambda);
}

inline void save_axisymmetric(const char* path,
                              const Kadath::Space_polar& space,
                              int kk,
                              double omega,
                              double lambda,
                              const Kadath::Scalar& nu,
                              const Kadath::Scalar& incA,
                              const Kadath::Scalar& incB,
                              const Kadath::Scalar& incbt,
                              const Kadath::Scalar& phi) {
    FILE* f = fopen(path, "w");
    if (!f) throw std::runtime_error(std::string("Cannot open for write: ") + path);
    space.save(f);
    Kadath::fwrite_be(&kk, sizeof(int), 1, f);
    Kadath::fwrite_be(&omega, sizeof(double), 1, f);
    Kadath::fwrite_be(&lambda, sizeof(double), 1, f);
    nu.save(f);
    incA.save(f);
    incB.save(f);
    incbt.save(f);
    phi.save(f);
    fclose(f);
}

inline void save_spherical(const char* path,
                           const Kadath::Space_polar& space,
                           double omega,
                           double lambda,
                           const Kadath::Scalar& psi,
                           const Kadath::Scalar& nu,
                           const Kadath::Scalar& phi) {
    FILE* f = fopen(path, "w");
    if (!f) throw std::runtime_error(std::string("Cannot open for write: ") + path);
    space.save(f);
    Kadath::fwrite_be(&omega, sizeof(double), 1, f);
    Kadath::fwrite_be(&lambda, sizeof(double), 1, f);
    psi.save(f);
    nu.save(f);
    phi.save(f);
    fclose(f);
}

} // namespace Io
