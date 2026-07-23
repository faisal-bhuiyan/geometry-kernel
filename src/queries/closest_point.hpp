#pragma once

#include <algorithm>
#include <cmath>

#include "core/point.hpp"
#include "core/predicates.hpp"
#include "core/tolerance.hpp"
#include "core/types.hpp"
#include "core/vector.hpp"

namespace geometry_kernel::queries {

using namespace geometry_kernel::core;

//---------------------------------------------------------------------------
// Point-to-line queries
//---------------------------------------------------------------------------

/**
 * @brief Projects a point onto the infinite line through @p a and @p b.
 *
 * Parametrize the line as q(t) = a + t * (b - a). The closest point is the
 * foot of the perpendicular from @p point, found by requiring the
 * displacement (point - q(t)) to be orthogonal to the line's direction:
 *
 *     (point - a - t * direction) · direction = 0
 *
 * Solving for t:
 *
 *     t = Dot(point - a, direction) / LengthSquared(direction)
 *
 * t is an affine parameter along the line, not an arc-length distance
 * (unless direction happens to be a unit vector): t = 0 at a, t = 1 at b,
 * and t outside [0, 1] is still a valid point on the infinite line.
 *
 *              point
 *                *
 *                |
 *                | (point - q) ⟂ direction
 *                v
 *   a ---------- q ---------------- b -------->  direction = b - a
 *   t=0          t                  t=1
 *
 * @param point Point to project.
 * @param a First point defining the line.
 * @param b Second point defining the line.
 * @return Closest point on the infinite line to @p point: q = a + t * (b - a).
 *
 * @note If @p a and @p b are coincident (within tolerance), returns @p a.
 *       The result may lie outside segment [a,b]; use ClosestPointOnSegment
 *       if you need it clamped to the segment.
 */
template <ScalarType T>
[[nodiscard]] inline Point2<T> ClosestPointOnLine(
    const Point2<T>& point, const Point2<T>& a, const Point2<T>& b
) {
    // Compute the direction vector from a to b
    const Vector2<T> direction{b - a};
    const T length_squared{LengthSquared(direction)};

    // Degenerate case -> a and b coincide -> return a
    if (RobustSign(length_squared) == 0) {
        return a;
    }

    // Compute the parameter t for the closest point on the line
    const T t{Dot(point - a, direction) / length_squared};
    return a + direction * t;  // a + t * (b - a)
}

/**
 * @brief Projects a point onto the closed segment [a,b], clamped to its endpoints.
 *
 * Same projection as ClosestPointOnLine, but the affine parameter t is
 * clamped to [0, 1] so the result cannot leave the segment:
 *
 *     t = clamp(Dot(point - a, direction) / LengthSquared(direction), 0, 1)
 *
 * When the unclamped t falls outside [0, 1], the perpendicular foot lies
 * beyond an endpoint, and the nearer endpoint becomes the answer instead.
 *
 * Case 1: foot lands inside [a, b] -> no clamping needed
 *
 *     point
 *       *
 *       |
 *       v
 * a-----q-----b---->  direction = b - a
 * t=0   t     t=1
 *
 * result = q = a + t * (b - a)
 *
 * Case 2: unclamped t < 0 (foot falls before a) -> clamped to a
 *
 *    point
 *      *
 *      |
 *      v
 * - - -q- - -a-----b---->  direction = b - a
 *      t<0   t=0   t=1
 *
 * result = a
 *
 * Case 3: unclamped t > 1 (foot falls beyond b) -> clamped to b
 *
 *           point
 *             *
 *             |
 *             v
 * a-----b- - -q---->  direction = b - a
 * t=0   t=1   t>1
 *
 * result = b
 *
 * @param point Point to project.
 * @param a Segment start.
 * @param b Segment end.
 * @return Closest point on segment [a,b] to @p point (may be an endpoint).
 *
 * @note If @p a and @p b are coincident (within tolerance), returns @p a.
 */
template <ScalarType T>
[[nodiscard]] inline Point2<T> ClosestPointOnSegment(
    const Point2<T>& point, const Point2<T>& a, const Point2<T>& b
) {
    const Vector2<T> direction{b - a};
    const T length_squared{LengthSquared(direction)};

    // Degenerate case -> a and b coincide -> return a
    if (RobustSign(length_squared) == 0) {
        return a;
    }

    // Compute the parameter t for the closest point on the segment
    T t{Dot(point - a, direction) / length_squared};
    t = std::clamp(t, T{0}, T{1});  // clamp t to [0, 1]
    return a + direction * t;       // a + t * (b - a) -> closest point on the segment
}

/**
 * @brief Squared distance from @p point to segment [a,b].
 *
 * Prefer this over DistancePointToSegment when only comparing relative distances
 * (e.g. finding the nearest of several segments), since it avoids the sqrt.
 */
template <ScalarType T>
[[nodiscard]] inline T PointToSegmentDistanceSquared(
    const Point2<T>& point, const Point2<T>& a, const Point2<T>& b
) {
    return LengthSquared(point - ClosestPointOnSegment(point, a, b));
}

/// @brief Distance from @p point to segment [a,b].
template <ScalarType T>
[[nodiscard]] inline T PointToSegmentDistance(
    const Point2<T>& point, const Point2<T>& a, const Point2<T>& b
) {
    return std::sqrt(PointToSegmentDistanceSquared(point, a, b));
}

//---------------------------------------------------------------------------
// Point-to-triangle queries
//---------------------------------------------------------------------------

/**
 * @brief Closest point on the filled triangle (v1, v2, v3) to @p point.
 *
 * Algorithm (O(1), no per-vertex Voronoi-region logic needed):
 *
 * 1. A triangle is convex. If @p point lies inside or on its boundary,
 *    the closest point in the filled region is @p point itself (distance zero).
 * 2. Otherwise, the closest point on any convex region's boundary always lies
 *    on one of its edges. So it suffices to project @p point onto each of the
 *    three edges via ClosestPointOnSegment and keep whichever projection has
 *    the smallest distance.
 * 3. Squared distances (PointToSegmentDistanceSquared) are compared instead of
 *    true distances to avoid three unnecessary sqrt calls -- squared distance
 *    preserves the ordering needed to pick the minimum.
 *
 * Note that clamped segment projections already handle the "closest to a
 * vertex" case correctly without extra logic: when @p point is nearest to a
 * shared vertex (e.g. v2), both adjacent edges (v1,v2) and (v2,v3) clamp their
 * own projection to v2, so the minimum-of-three naturally converges there.
 *
 * Case 1: point inside the triangle -> return point itself
 *
 *          v3
 *          /\
 *         /  \
 *        / *  \      point (*) is inside -> closest = point
 *       /______\
 *     v1        v2
 *
 * Case 2: point outside, closest to an edge interior
 *
 *          v3
 *          /\
 *         /  \
 *        /    \
 *       /______\
 *     v1        v2
 *           *            point outside, below edge v1-v2
 *           |
 *           v
 *        closest = interior point on edge v1-v2 (not a vertex)
 *
 * Case 3: point outside, closest to a shared vertex
 *
 *          v3
 *          /\
 *         /  \
 *        /    \
 *       /______\
 *     v1        v2  *    point outside, nearest to v2
 *
 *   Both edge (v1,v2) and edge (v2,v3) clamp their projection to v2,
 *   so either one yields the correct closest = v2.
 *
 * @param point Point to query.
 * @param v1 First vertex of the triangle.
 * @param v2 Second vertex of the triangle.
 * @param v3 Third vertex of the triangle.
 * @return Closest point on the triangle (interior or boundary) to @p point.
 *
 * @note Winding-order agnostic, since PointInTriangle and ClosestPointOnSegment both are.
 */
template <ScalarType T>
[[nodiscard]] inline Point2<T> ClosestPointOnTriangle(
    const Point2<T>& point, const Point2<T>& v1, const Point2<T>& v2, const Point2<T>& v3
) {
    // Interior (or boundary) case
    if (PointInTriangle(point, v1, v2, v3)) {
        return point;
    }

    // Edge v1-v2
    Point2<T> closest{ClosestPointOnSegment(point, v1, v2)};
    T best_distance_squared{PointToSegmentDistanceSquared(point, v1, v2)};

    // Edge v2-v3
    if (const T d23{PointToSegmentDistanceSquared(point, v2, v3)}; d23 < best_distance_squared) {
        closest = ClosestPointOnSegment(point, v2, v3);
        best_distance_squared = d23;
    }

    // Edge v3-v1
    if (const T d31{PointToSegmentDistanceSquared(point, v3, v1)}; d31 < best_distance_squared) {
        closest = ClosestPointOnSegment(point, v3, v1);
        best_distance_squared = d31;
    }

    return closest;
}

/// @brief Distance from @p point to the filled triangle (v1, v2, v3). Zero if point is inside.
template <ScalarType T>
[[nodiscard]] inline T PointToTriangleDistance(
    const Point2<T>& point, const Point2<T>& v1, const Point2<T>& v2, const Point2<T>& v3
) {
    return Length(point - ClosestPointOnTriangle(point, v1, v2, v3));
}

}  // namespace geometry_kernel::queries
