#pragma once

#include <cmath>

#include "tolerance.hpp"
#include "types.hpp"

namespace geometry_kernel::core {

//---------------------------------------------------------------------------
// Point in 3D space
//---------------------------------------------------------------------------

/// @brief 3D point in space
template <ScalarType T>
struct Point3 {
    T x{};
    T y{};
    T z{};
};

using Point3D = Point3<double>;  ///< 3D point in space with double precision

/**
 * @brief Returns true if two points are within @p kTolerance of each other in all coordinates.
 *
 * Used to detect duplicate vertices without exact floating-point equality.
 *
 * @param a First point.
 * @param b Second point.
 * @return True if |a.x - b.x| <= kTolerance and likewise for y, z.
 */
template <ScalarType T>
[[nodiscard]] inline bool NearlyEqualPoint(const Point3<T>& a, const Point3<T>& b) {
    return (std::abs(a.x - b.x) <= kTolerance<T>) && (std::abs(a.y - b.y) <= kTolerance<T>) &&
           (std::abs(a.z - b.z) <= kTolerance<T>);
}

}  // namespace geometry_kernel::core
