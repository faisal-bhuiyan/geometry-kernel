#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include "core/point.hpp"
#include "core/predicates.hpp"
#include "core/tolerance.hpp"
#include "core/types.hpp"
#include "queries/intersection.hpp"

namespace geometry_kernel::queries {

using namespace geometry_kernel::core;

//---------------------------------------------------------------------------
// Triangle containment (3D, in-plane)
//---------------------------------------------------------------------------

/**
 * @brief Tests if the in-plane projection of @p point lies inside (or on the boundary
 *        of) triangle (v1, v2, v3), given the triangle's own face @p normal.
 *
 * The 3D counterpart to PointInTriangle: runs the same three-edge sign-agreement test,
 * using EdgeSign3 (a half-space test oriented by @p normal) in place of
 * SignedTriangleArea2. Because EdgeSign3 only reads point's in-plane position (see its
 * doc comment), @p point need not actually lie in the triangle's plane -- this tests
 * whether its perpendicular projection onto that plane would fall inside the triangle.
 *
 * @param point Point to test (need not lie in the triangle's plane)
 * @param v1 First vertex of the triangle
 * @param v2 Second vertex of the triangle
 * @param v3 Third vertex of the triangle
 * @param normal The triangle's face normal, e.g. Cross(v2 - v1, v3 - v1). Passed in
 *        rather than recomputed, since callers (e.g. ClosestPointOnTriangle) typically
 *        already need it for the plane-projection step
 * @return True if point's in-plane projection is inside or on the boundary of the triangle
 *
 * @note Winding-order agnostic and boundary-inclusive, same as the 2D version.
 */
template <ScalarType T>
[[nodiscard]] inline bool InPlaneTriangleContainment(
    const Point3<T>& point, const Point3<T>& v1, const Point3<T>& v2, const Point3<T>& v3,
    const Vector3<T>& normal
) {
    const int s1{RobustSign(EdgeSign3(v1, v2, point, normal))};
    const int s2{RobustSign(EdgeSign3(v2, v3, point, normal))};
    const int s3{RobustSign(EdgeSign3(v3, v1, point, normal))};

    const bool has_negative{(s1 < 0) || (s2 < 0) || (s3 < 0)};
    const bool has_positive{(s1 > 0) || (s2 > 0) || (s3 > 0)};
    return !(has_negative && has_positive);
}

}  // namespace geometry_kernel::queries
