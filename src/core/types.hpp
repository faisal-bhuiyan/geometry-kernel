#pragma once

#include <concepts>

namespace geometry_kernel::core {

//---------------------------------------------------------------------------
// Type aliases and helpers
//---------------------------------------------------------------------------

/**
 * @brief Concept satisfied by any standard floating-point type.
 *
 * Constrains template parameters to `float`, `double`, or `long double`,
 * ensuring signed arithmetic, fractional representation, and well-defined
 * behavior for edge cases (infinity, NaN) as required by IEEE 754.
 *
 * @tparam T A type satisfying `std::floating_point<T>`.
 * @see std::floating_point
 */
template <typename T>
concept ScalarType = std::floating_point<T>;

/// @brief Enumeration for the orientation of ordered triple (a, b, c) in the plane.
enum class Orientation {
    kCollinear = 0,         ///< Collinear orientation -> points lie on a line
    kCounterClockwise = 1,  ///< Counter-clockwise orientation -> points form a CCW polygon
    kClockwise = -1         ///< Clockwise orientation -> points form a CW polygon
};

}  // namespace geometry_kernel::core
