#pragma once

#include "types.hpp"

namespace geometry_kernel::core {

/// @brief Numerical tolerance for floating-point comparisons.
template <ScalarType T>
constexpr T kTolerance{static_cast<T>(1e-9)};

/**
 * @brief Robust sign test for floating-point values.
 *
 * Classifies a value as positive (+1), negative (-1), or zero (0) with tolerance
 * to avoid spurious results from floating-point error. Values inside the "dead zone"
 * zone [-tolerance, +tolerance] snap to zero.
 *
 *              -1                    0                   +1
 *
 *   <-----------------------[========|========]----------------------->
 *                   -kTolerance     0.0     +kTolerance
 *
 *                           |<-- dead zone -->|
 *                            both ends included
 *
 * NOT TO SCALE -> the band defaults to 1e-9 wide while the outer regions run to +/-infinity.
 *
 * --------------------------------------------------------------------------------------
 * ** WHY A DEAD ZONE IS NECESSARY? **
 *
 * This kernel decides degeneracy by asking whether a determinant is zero. Three
 * exactly-collinear points have a doubled area of exactly zero in real arithmetic,
 * but in floating point the subtractions and multiplications leave a residue of a
 * few ULP, so the value arrives as a tiny nonzero number whose sign is arbitrary.
 * A bare `== 0` test would almost never fire, and that arbitrary sign would
 * surface as a confident CW or CCW answer from TriangleOrientation and as a
 * wrong side-of-edge vote in PointInTriangle. Snapping the band to zero turns
 * "close enough to degenerate" into a decidable, repeatable outcome.
 *
 * --------------------------------------------------------------------------------------
 *
 * @param value Scalar to classify.
 * @param tolerance Half-width of the dead zone (defaults to kTolerance). Pass 0
 *        to recover exact sign semantics, where only a true zero returns 0.
 * @return -1, 0, or +1. Exactly +/-tolerance returns 0 -> the band is closed.
 *
 * @note For T = float the default band is narrower than the type's own epsilon
 *       (1e-9 against roughly 1.2e-7), so at coordinates near 1 a single
 *       rounding step already escapes it and the dead zone does nothing. Use
 *       double, or supply a tolerance suited to float.
 *
 * @see Shewchuk, J.R. "Adaptive Precision Floating-Point Arithmetic and
 *      Fast Robust Geometric Predicates." Discrete & Computational Geometry,
 *      18(3):305-363, 1997.
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
