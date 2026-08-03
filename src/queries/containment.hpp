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
// Triangle containment
//---------------------------------------------------------------------------

/**
 * @brief Tests if a point lies inside or on the boundary of triangle (v1, v2, v3).
 *
 * For each directed edge of the triangle, computes which side the point falls on
 * using the signed area of the triangle formed by that edge and the point.
 * A point is inside if and only if all three signs agree (all positive, all negative,
 * or zero for boundary). Mixed signs mean the point is outside at least one edge.
 *
 * Winding-order agnostic -> works for both CCW and CW triangles.
 * Boundary-inclusive -> points exactly on an edge or vertex return true.
 *
 * Case 1: point strictly inside -> true (all three signs agree)
 *
 *          v3
 *          /\
 *         /  \
 *        / *  \      * = point
 *       /______\
 *     v1        v2
 *
 * Case 2: point on an edge (or vertex) -> true (at least one sign is zero;
 *         the other two do not disagree)
 *
 *          v3
 *          /\
 *         /  \
 *        /    \
 *       /__*___\
 *     v1        v2
 *
 * Case 3: point outside -> false (mixed signs: on the "inside" of some edges
 *         and the "outside" of at least one)
 *
 *          v3
 *          /\
 *         /  \
 *        /    \
 *       /______\
 *     v1        v2
 *           *
 *
 * @param point Point to test.
 * @param v1 First vertex of the triangle.
 * @param v2 Second vertex of the triangle.
 * @param v3 Third vertex of the triangle.
 * @return True if @p point is inside or on the boundary of triangle (v1, v2, v3).
 */
template <ScalarType T>
[[nodiscard]] inline bool PointInTriangle(
    const Point2<T>& point, const Point2<T>& v1, const Point2<T>& v2, const Point2<T>& v3
) {
    const int s1{RobustSign(SignedTriangleArea2(v1, v2, point))};
    const int s2{RobustSign(SignedTriangleArea2(v2, v3, point))};
    const int s3{RobustSign(SignedTriangleArea2(v3, v1, point))};

    const bool has_negative{(s1 < 0) || (s2 < 0) || (s3 < 0)};
    const bool has_positive{(s1 > 0) || (s2 > 0) || (s3 > 0)};
    return !(has_negative && has_positive);
}

}  // namespace geometry_kernel::queries
