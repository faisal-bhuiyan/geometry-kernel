#pragma once

#include "types.hpp"

namespace geometry_kernel::core {

/// @brief Numerical tolerance for floating-point comparisons.
template <ScalarType T>
constexpr T kTolerance{static_cast<T>(1e-9)};

/**
 * @brief Robust sign test for floating-point values.
 *
 * Classifies a value as positive (+1), negative (-1), or zero (0) with
 * tolerance to avoid spurious results from floating-point error. Values
 * inside the dead zone [-tolerance, +tolerance] snap to zero:
 *
 *      -1                  0                  +1
 *
 *  <-------|---------------|---------------|------->
 *       -kTol            0.0             +kTol
 *
 *  value < -kTol   ->  -1
 *  |value| <= kTol ->   0   (dead zone)
 *  value >  kTol   ->  +1
 *
 * @param value Scalar to classify.
 * @param tolerance Half-width of the dead zone (defaults to kTolerance).
 * @return -1, 0, or +1.
 *
 * @see Shewchuk, J.R. "Adaptive Precision Floating-Point Arithmetic and
 *      Fast Robust Geometric Predicates." Discrete & Computational Geometry,
 *      18(3):305–363, 1997.
 */
template <ScalarType T>
inline int RobustSign(T value, T tolerance = kTolerance<T>) {
    // Strict inequalities: values at exactly +/-tolerance round to zero
    if (value > tolerance) {
        return 1;
    }
    if (value < -tolerance) {
        return -1;
    }
    return 0;
}

}  // namespace geometry_kernel::core
