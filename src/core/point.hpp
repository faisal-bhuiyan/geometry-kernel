#pragma once

#include <cmath>

#include "tolerance.hpp"
#include "types.hpp"

namespace geometry_kernel::core {

//---------------------------------------------------------------------------
// Type aliases and helpers
//---------------------------------------------------------------------------

/// @brief 3D point in space.
template <ScalarType T>
struct Point {
    T x{};
    T y{};
    T z{};
};

using Point2D = Point<double>;  ///< 2D point in the plane with double precision

/**
 * @brief Returns true if two points are within @p kTolerance of each other in both coordinates.
 *
 * Used to detect duplicate vertices without exact floating-point equality.
 *
 * @param a First point.
 * @param b Second point.
 * @return True if |a.x - b.x| <= kTolerance and |a.y - b.y| <= kTolerance and |a.z - b.z| <=
 * kTolerance.
 */
template <ScalarType T>
[[nodiscard]] inline bool NearlyEqualPoint(const Point<T>& a, const Point<T>& b) {
    return (std::abs(a.x - b.x) <= kTolerance<T>) && (std::abs(a.y - b.y) <= kTolerance<T>) &&
           (std::abs(a.z - b.z) <= kTolerance<T>);
}

}  // namespace geometry_kernel::core
