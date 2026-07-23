#pragma once

#include <algorithm>
#include <cmath>

#include "core/point.hpp"
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

}  // namespace geometry_kernel::queries
