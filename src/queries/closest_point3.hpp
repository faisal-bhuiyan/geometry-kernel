#pragma once

#include <algorithm>
#include <cmath>

#include "core/point3.hpp"
#include "core/predicates3.hpp"
#include "core/tolerance.hpp"
#include "core/types.hpp"
#include "core/vector3.hpp"
#include "queries/containment3.hpp"

namespace geometry_kernel::queries {

using namespace geometry_kernel::core;

//---------------------------------------------------------------------------
// Point-to-segment queries (3D, in-plane)
//---------------------------------------------------------------------------

/**
 * @brief Projects a point onto the closed segment [a,b] in 3D, clamped to its endpoints.
 *
 * Same projection and clamping logic as the 2D ClosestPointOnSegment -> underlying
 * math (Dot, LengthSquared, affine parameter t clamped to [0,1]) is dimension-agnostic,
 * so this is a direct extension rather than a new derivation.
 *
 * @param point Point to project
 * @param a Segment start
 * @param b Segment end
 * @return Closest point on segment [a,b] to @p point (may be an endpoint)
 *
 * @note If @p a and @p b are coincident (within tolerance), returns @p a.
 */
template <ScalarType T>
[[nodiscard]] inline Point3<T> ClosestPointOnSegment(
    const Point3<T>& point, const Point3<T>& a, const Point3<T>& b
) {
    const Vector3<T> direction{b - a};
    const T length_squared{LengthSquared(direction)};
    if (RobustSign(length_squared) == 0) {
        return a;
    }
    T t{Dot(point - a, direction) / length_squared};
    t = std::clamp(t, T{0}, T{1});
    return a + direction * t;
}

/// @brief Squared distance from @p point to segment [a,b] in 3D
template <ScalarType T>
[[nodiscard]] inline T PointToSegmentDistanceSquared(
    const Point3<T>& point, const Point3<T>& a, const Point3<T>& b
) {
    return LengthSquared(point - ClosestPointOnSegment(point, a, b));
}

/// @brief Distance from @p point to segment [a,b] in 3D
template <ScalarType T>
[[nodiscard]] inline T PointToSegmentDistance(
    const Point3<T>& point, const Point3<T>& a, const Point3<T>& b
) {
    return std::sqrt(PointToSegmentDistanceSquared(point, a, b));
}

//---------------------------------------------------------------------------
// Point-to-triangle queries (3D, in-plane)
//---------------------------------------------------------------------------

/**
 * @brief Closest point on the filled triangle (v1, v2, v3) to @p point, in 3D.
 *
 * --------------------------------------------------------------------------------------
 *
 * Unlike the 2D version, @p point generally does not lie in the triangle's plane, so
 * "inside the triangle" first means "does point's perpendicular projection onto the
 * triangle's plane fall inside the triangle" (InPlaneTriangleContainment). If so, that
 * projection -- not point itself -- is the closest point (point is pulled straight down
 * onto the plane along the normal). Otherwise, same as 2D: the closest point on a
 * convex region's boundary lies on one of its edges, so fall back to the minimum over
 * ClosestPointOnSegment across the three edges, using the true 3D point (not the
 * projection).
 *
 * Each case below is drawn twice: a side view with the plane seen edge-on, and a top
 * view looking down the normal. The side view shows the drop; the top view shows the
 * containment test that chooses the branch. The side view cannot decide it on its own,
 * since landing within the triangle's silhouette in one direction does not put you
 * inside the triangle.
 *
 * --------------------------------------------------------------------------------------
 *
 * Case 1: the foot lands inside the triangle -> the foot is the answer
 *
 *   side view                                   top view (down the normal)
 *
 *                   * point                              * v3
 *                   :                                   / \
 *                   :  drop along the normal           /   \
 *                   :                                 /     \
 *                   v                                /   *   \
 *   ---[============*============]---               /    f    \
 *                   f                           v1 *-----------* v2
 *
 *   [===] = the triangle, seen edge-on            f is inside -> result = f
 *
 * Case 2: the foot lands outside -> clamp to the nearest edge
 *
 *   side view                                   top view (down the normal)
 *
 *                               * point                 * v3
 *                               :                      / \
 *                               :                     /   \
 *                               :                    /     \
 *                               v                   /       \
 *   ---[=================]------*------          v1 *-----*-----* v2
 *                       q       f                         q
 *                                                         :
 *   the drop overshoots the triangle                      * f
 *
 *                                                 f is outside -> result = q,
 *                                                 the nearest point on edge v1-v2
 *
 * The fallback measures from the true 3D @p point rather than from f, but the two agree:
 * (point - f) runs along the normal, which is perpendicular to every in-plane edge
 * direction, so it contributes nothing to the projection parameter along an edge. That
 * is why the outside branch never needs to compute f at all.
 *
 * --------------------------------------------------------------------------------------
 *
 * @param point Point to query
 * @param v1 First vertex of the triangle
 * @param v2 Second vertex of the triangle
 * @param v3 Third vertex of the triangle
 * @return Closest point on the triangle (interior or boundary) to @p point
 *
 * @note Winding-order agnostic, same as the 2D version.
 * @see Ericson, C. Real-Time Collision Detection. Morgan Kaufmann, 2005, §5.1–5.4.
 */
template <ScalarType T>
[[nodiscard]] inline Point3<T> ClosestPointOnTriangle(
    const Point3<T>& point, const Point3<T>& v1, const Point3<T>& v2, const Point3<T>& v3
) {
    const Vector3<T> normal{Cross(v2 - v1, v3 - v1)};

    if (InPlaneTriangleContainment(point, v1, v2, v3, normal)) {
        // Project point onto the triangle's plane along the normal:
        //   t = Dot(point - v1, normal) / LengthSquared(normal)
        // gives the signed distance (in units of |normal|) from the plane; subtracting
        // t * normal from point removes exactly the out-of-plane component.
        const T t{Dot(point - v1, normal) / LengthSquared(normal)};
        return point + normal * (-t);
    }

    // Outside the triangle's footprint -> closest point is on one of the three edges
    Point3<T> closest{ClosestPointOnSegment(point, v1, v2)};
    T best_distance_squared{PointToSegmentDistanceSquared(point, v1, v2)};

    if (const T d23{PointToSegmentDistanceSquared(point, v2, v3)}; d23 < best_distance_squared) {
        closest = ClosestPointOnSegment(point, v2, v3);
        best_distance_squared = d23;
    }

    if (const T d31{PointToSegmentDistanceSquared(point, v3, v1)}; d31 < best_distance_squared) {
        closest = ClosestPointOnSegment(point, v3, v1);
        best_distance_squared = d31;
    }

    return closest;
}

/// @brief Distance from @p point to the filled triangle (v1, v2, v3) in 3D
template <ScalarType T>
[[nodiscard]] inline T PointToTriangleDistance(
    const Point3<T>& point, const Point3<T>& v1, const Point3<T>& v2, const Point3<T>& v3
) {
    return Length(point - ClosestPointOnTriangle(point, v1, v2, v3));
}

}  // namespace geometry_kernel::queries
