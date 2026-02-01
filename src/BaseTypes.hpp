#ifndef BASETYPES_HPP
#define BASETYPES_HPP

/** \file BaseTypes.hpp
 * \brief Basic type definitions and global constants */

#include <armadillo>

/// Type aliases for Armadillo types
using Matrix = arma::mat;
using Vector = arma::vec;
using Vector_un = arma::Col<unsigned>;

/// Mathematical constants and utilities
namespace math {

/// Pi constant
inline constexpr double pi = 3.14159265358979323846;

/// 1/(2*pi) - used frequently in Biot-Savart calculations
inline constexpr double rpi2 = 1.0 / (2.0 * 3.14159265358979323846);

/// Square a value
template<typename T>
[[nodiscard]] constexpr T sqr(T x) noexcept {
    return x * x;
}

} // namespace math

#endif // BASETYPES_HPP
