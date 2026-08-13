#pragma once

#include <algorithm>

#include "core/point2.hpp"
#include "core/predicates2.hpp"
#include "core/tolerance.hpp"
#include "core/types.hpp"

namespace geometry_kernel::queries {

using namespace geometry_kernel::core;

//---------------------------------------------------------------------------
// Line segment containment
//---------------------------------------------------------------------------

/**
 * @brief Tests if @p point lies on the closed segment [a,b].
 *
 * Assumes collinearity has already been established. This function only checks
 * that @p point falls inside the axis-aligned bounding box of [a,b], expanded by
 * kTolerance (so endpoints and near-endpoint floating-point noise still pass).
 *
 * --------------------------------------------------------------------------------------
 *
 * Case 1: point between a and b -> true
 *
 *     a-----*-----b
 *           point
 *
 * Case 2: point at an endpoint -> true (closed segment)
 *
 *     *-----------b
 *     a (= point)
 *
 * Case 3: point collinear but outside [a,b] -> false
 *
 *     * - - a-----b
 *   point
 *   (on the infinite line through a,b, but outside the segment extent)
 *
 * --------------------------------------------------------------------------------------
 *
 * @param point Query point (must already be collinear with a and b).
 * @param a Segment start.
 * @param b Segment end.
 * @return True if @p point's coordinates lie within the padded AABB of [a,b].
 *
 * @pre @p point is already known to be collinear with @p a and @p b (e.g. via a prior
 *      RobustSign(SignedTriangleAreaTimes2(a, b, point)) == 0 check). This function only
 *      performs the bounding-box containment check; it does not re-verify collinearity,
 *      since every call site in this file has already established it.
 */
template <ScalarType T>
[[nodiscard]] inline bool PointOnSegment(
    const Point2<T>& point, const Point2<T>& a, const Point2<T>& b
) {
    return point.x >= std::min(a.x, b.x) - kTolerance<T> &&
           point.x <= std::max(a.x, b.x) + kTolerance<T> &&
           point.y >= std::min(a.y, b.y) - kTolerance<T> &&
           point.y <= std::max(a.y, b.y) + kTolerance<T>;
}

//---------------------------------------------------------------------------
// Line segment intersection
//---------------------------------------------------------------------------

/**
 * @brief Tests if closed segments [a,b] and [c,d] intersect.
 *
 * Uses the standard orientation test with a collinear fallback:
 *
 * 1. Proper (crossing) intersection: c and d fall on strictly opposite sides of
 *    line ab, AND a and b fall on strictly opposite sides of line cd
 *    (ab_c * ab_d < 0 && cd_a * cd_b < 0).
 *
 * 2. Collinear / touching: an endpoint has sign zero relative to the other
 *    segment's line -> fall back to PointOnSegment to confirm it lies within
 *    that segment's extent, not just on its infinite line.
 *
 * --------------------------------------------------------------------------------------
 *
 * Case 1: proper crossing -> true
 *
 *          c
 *          |
 *    a-----+-----b
 *          |
 *          d
 *
 * Case 2: endpoint of one segment lies on the other -> true
 *
 *    a-----c-----b
 *          |
 *          d
 *
 *    (c lies on [a,b]; ab_c == 0 and PointOnSegment(c, a, b))
 *
 * Case 3: shared endpoint -> true
 *
 *    a-----b
 *          |
 *          d
 *
 *    (b == c; closed segments include their endpoints)
 *
 * Case 4: collinear but disjoint -> false
 *
 *    a-----b   c-----d
 *
 *    (endpoints have sign zero, but PointOnSegment rejects each as outside)
 *
 * Case 5: non-crossing (same-side) -> false
 *
 *          c
 *         /
 *        d
 *
 *    a-----b
 *
 *    (c and d are both above ab, so ab_c * ab_d > 0; segments do not cross)
 *
 * --------------------------------------------------------------------------------------
 *
 * @param a First segment start.
 * @param b First segment end.
 * @param c Second segment start.
 * @param d Second segment end.
 * @return True if the closed segments [a,b] and [c,d] share any point.
 */
template <ScalarType T>
[[nodiscard]] inline bool SegmentsIntersect(
    const Point2<T>& a, const Point2<T>& b, const Point2<T>& c, const Point2<T>& d
) {
    // Orientations of each endpoint relative to the other segment's supporting line
    const int ab_c{RobustSign(SignedTriangleAreaTimes2(a, b, c))};
    const int ab_d{RobustSign(SignedTriangleAreaTimes2(a, b, d))};
    const int cd_a{RobustSign(SignedTriangleAreaTimes2(c, d, a))};
    const int cd_b{RobustSign(SignedTriangleAreaTimes2(c, d, b))};

    // Case 1: proper crossing —> endpoints on strictly opposite sides of each line
    if (ab_c * ab_d < 0 && cd_a * cd_b < 0) {
        return true;
    }

    // Cases 2–4: collinear / touching -> an endpoint lies on the other line
    // We need to confirm it also lies within that segment's extent
    if (ab_c == 0 && PointOnSegment(c, a, b)) {
        return true;
    }
    if (ab_d == 0 && PointOnSegment(d, a, b)) {
        return true;
    }
    if (cd_a == 0 && PointOnSegment(a, c, d)) {
        return true;
    }
    if (cd_b == 0 && PointOnSegment(b, c, d)) {
        return true;
    }

    // Case 5: disjoint or same-side non-crossing -> segments do not cross
    return false;
}

}  // namespace geometry_kernel::queries
