#pragma once

#include "types.hpp"

namespace geometry_kernel::core {

/// @brief Numerical tolerance for floating-point comparisons.
template <ScalarType T>
constexpr T kTolerance{static_cast<T>(1e-9)};

/**
 * @brief Robust sign test for floating-point values.
 *
 * Classifies a value as positive (+1), negative (-1), or zero (0) with tolerance to avoid spurious
 * results from floating-point error.
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
