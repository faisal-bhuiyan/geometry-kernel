#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include "core/point2.hpp"
#include "core/predicates2.hpp"
#include "core/tolerance.hpp"
#include "core/types.hpp"
#include "queries/intersection.hpp"

namespace geometry_kernel::queries {

using namespace geometry_kernel::core;

//---------------------------------------------------------------------------
// Triangle containment (2D)
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
 * --------------------------------------------------------------------------------------
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
 * --------------------------------------------------------------------------------------
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
    const int s1{RobustSign(SignedTriangleAreaTimes2(v1, v2, point))};
    const int s2{RobustSign(SignedTriangleAreaTimes2(v2, v3, point))};
    const int s3{RobustSign(SignedTriangleAreaTimes2(v3, v1, point))};

    const bool has_negative{(s1 < 0) || (s2 < 0) || (s3 < 0)};
    const bool has_positive{(s1 > 0) || (s2 > 0) || (s3 > 0)};
    return !(has_negative && has_positive);
}

//---------------------------------------------------------------------------
// Polygon containment (2D)
//---------------------------------------------------------------------------

/**
 * @brief Tests if a point lies inside or on the boundary of a simple polygon.
 *
 * --------------------------------------------------------------------------------------
 *
 * Two-phase algorithm:
 *
 * 1. Boundary: for each edge, if the point is collinear with that edge and
 *    PointOnSegment reports it lies within the edge's extent, return true.
 * 2. Interior: crossing-number (ray casting) — cast a horizontal ray from the
 *    point in the +x direction and count edges that cross it. An odd count
 *    means inside. Correct for convex and concave (non-self-intersecting)
 *    polygons.
 *
 * Case 1: point inside (odd number of ray crossings) -> true
 *
 *     v4------------v3
 *      |             |
 *      |  *--------->|----  ray in +x: 1 crossing -> odd -> inside
 *      |             |
 *     v1------------v2
 *
 * Case 2: point outside (even number of ray crossings) -> false
 *
 *     v4------------v3
 *      |             |
 *      |             |  *--------->  ray in +x: 0 crossings -> even -> outside
 *      |             |
 *     v1------------v2
 *
 * Case 3: point on the boundary -> true (caught in phase 1 before ray casting)
 *
 *     v4------------v3
 *      |             |
 *      *             |      * lies on edge v1-v4
 *      |             |
 *     v1------------v2
 *
 * Case 4: concave polygon -> ray casting still works (odd/even)
 *
 *     v5----v4
 *      |     \
 *      |  *   v3----v2      * inside the concave pocket:
 *      |               \     ray crosses one edge -> odd -> inside
 *     v0----------------v1
 *
 * --------------------------------------------------------------------------------------
 *
 * @param point Point to test.
 * @param polygon Vertices of a simple polygon as an open ring (last vertex does
 *                not repeat the first).
 * @return True if @p point is inside or on the boundary of @p polygon.
 *
 * @pre polygon.size() >= 3; fewer vertices returns false.
 * @note Self-intersecting polygons produce undefined results.
 */
template <ScalarType T>
[[nodiscard]] inline bool PointInPolygon(
    const Point2<T>& point, std::span<const Point2<T>> polygon
) {
    const std::size_t num_vertices{polygon.size()};
    if (num_vertices < 3U) {
        return false;
    }

    // Phase 1: boundary -> collinear with an edge and within its extent
    for (std::size_t i = 0; i < num_vertices; ++i) {
        const Point2<T>& v0{polygon[i]};
        const Point2<T>& v1{polygon[(i + 1) % num_vertices]};
        if (RobustSign(SignedTriangleAreaTimes2(v0, v1, point)) == 0 &&
            PointOnSegment(point, v0, v1)) {
            return true;
        }
    }

    // Phase 2: crossing-number (ray cast in +x) -> toggle inside on each crossing
    bool inside{false};
    for (std::size_t i = 0, j = num_vertices - 1; i < num_vertices; j = i++) {
        const Point2<T>& vi{polygon[i]};
        const Point2<T>& vj{polygon[j]};
        const bool straddles_y{(vi.y > point.y) != (vj.y > point.y)};
        if (straddles_y) {
            const T x_at_ray{(vj.x - vi.x) * (point.y - vi.y) / (vj.y - vi.y) + vi.x};
            if (point.x < x_at_ray) {
                inside = !inside;
            }
        }
    }
    return inside;
}

/// @brief Convenience overload -> accepts std::vector<Point2<T>> directly
template <ScalarType T>
[[nodiscard]] inline bool PointInPolygon(
    const Point2<T>& point, const std::vector<Point2<T>>& polygon
) {
    return PointInPolygon(point, std::span<const Point2<T>>{polygon});
}

}  // namespace geometry_kernel::queries
