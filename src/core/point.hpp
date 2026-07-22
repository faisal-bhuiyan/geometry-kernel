#pragma once

#include <cmath>

#include "tolerance.hpp"
#include "types.hpp"

namespace geometry_kernel::core {

//---------------------------------------------------------------------------
// Point in the plane (2D)
//---------------------------------------------------------------------------

/// @brief 2D point in the plane
template <ScalarType T>
struct Point2 {
    T x{};
    T y{};
};

using Point2D = Point2<double>;  ///< 2D point in the plane with double precision

/**
 * @brief Returns true if two points are within @p kTolerance of each other in all coordinates.
 *
 * Used to detect duplicate vertices without exact floating-point equality.
 *
 * @param a First point.
 * @param b Second point.
 * @return True if |a.x - b.x| <= kTolerance and |a.y - b.y| <= kTolerance.
 */
template <ScalarType T>
[[nodiscard]] inline bool NearlyEqualPoint(const Point2<T>& a, const Point2<T>& b) {
    return (std::abs(a.x - b.x) <= kTolerance<T>) && (std::abs(a.y - b.y) <= kTolerance<T>);
}

}  // namespace geometry_kernel::core
