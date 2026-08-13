#pragma once

#include <span>
#include <vector>

#include "point2.hpp"
#include "types.hpp"
#include "vector2.hpp"

namespace geometry_kernel::core {

//---------------------------------------------------------------------------
// Triangle
//---------------------------------------------------------------------------

/**
 * @brief Signed 2D cross product: z-component of vector (b - a) x (c - a).
 *
 * Equals twice the signed area of triangle (a, b, c), or equivalently the
 * signed area of the parallelogram spanned by u = (b - a) and v = (c - a).
 * In 2D the cross product is the scalar:
 *
 *     u × v = u.x * v.y - u.y * v.x
 *
 * Why that scalar is an area: rotating u by 90° CCW gives perp(u) = (-u.y, u.x),
 * and Dot(perp(u), v) = u.x * v.y - u.y * v.x = u × v. So the cross product is
 * |u| times the component of v perpendicular to u -> base × height =
 * area of the parallelogram spanned by u and v.
 *
 *         c *---------------* c + u
 *          /:              /
 *         / :             /        h = |v| sin θ
 *    v   /  :            /             (perpendicular distance from c
 *       /   :           /               to the line through a and b)
 *      /θ   :          /
 *     *-----+---------*
 *     a         u      b
 *
 *     u × v = |u| · h = |u| |v| sin θ = parallelogram area
 *
 * Why "times 2": the diagonal b–c cuts that parallelogram into two congruent
 * triangles -> a 180° rotation about the midpoint of bc maps T1 onto T2 -> and
 * T1 is exactly triangle (a, b, c):
 *
 *         c *---------------* c + u
 *          /    ·    T2    /
 *     v   /   T1    ·     /
 *        *---------------*
 *        a       u        b
 *
 *     u × v = 2 × area of triangle (a, b, c)
 *
 * Sign follows the right-hand rule: positive when v is CCW from u (c lies left
 * of the directed line a -> b), negative when CW, and exactly zero when a, b, c
 * are collinear -> the degenerate parallelogram has no area.
 *
 * @param a First vertex.
 * @param b Second vertex.
 * @param c Third vertex.
 * @return Signed parallelogram area in the xy-plane (twice the signed triangle area).
 *
 * @see https://en.wikipedia.org/wiki/Cross_product#Computational_geometry
 */
template <ScalarType T>
[[nodiscard]] inline T SignedTriangleAreaTimes2(
    const Point2<T>& a, const Point2<T>& b, const Point2<T>& c
) {
    return Cross(b - a, c - a);
}

/// @brief Area of triangle (a, b, c) in the plane
template <ScalarType T>
[[nodiscard]] inline T TriangleArea(const Point2<T>& a, const Point2<T>& b, const Point2<T>& c) {
    return static_cast<T>(0.5) * std::abs(SignedTriangleAreaTimes2(a, b, c));
}

//---------------------------------------------------------------------------
// Polygon
//---------------------------------------------------------------------------

/**
 * @brief Computes the signed area of a simple polygon using the shoelace formula.
 *
 * Each term of the sum is p0.x * p1.y - p1.x * p0.y = Cross(p0, p1), which is
 * SignedTriangleAreaTimes2 evaluated with the origin as the first vertex:
 *
 *     Cross(vi, vi+1) = SignedTriangleAreaTimes2(O, vi, vi+1)
 *
 * So the shoelace sum is a fan of triangles from the origin O out to every edge,
 * each counted as twice its signed area -> which is what the trailing 0.5
 * factor undoes. The fan is correct even when O lies outside the polygon,
 * because the surplus region is swept once with each sign:
 *
 *        v4 *-------------* v3
 *           |             |
 *           |   polygon   |
 *        v1 *-------------* v2
 *            \  ·       · /
 *             \  ·     · /     edge (v1,v2) winds CW about O  -> negative
 *              \  ·   · /      edge (v3,v4) winds CCW about O -> positive
 *               \  · · /
 *                \ · ·/
 *                  O
 *
 * The strip between O and the polygon is covered twice with opposite signs and
 * cancels exactly, leaving the polygon's own area. Edges that are collinear
 * with O contribute zero.
 *
 * Reversing traversal negates every term, so the total sign encodes winding:
 *
 *   v1 -> v2 -> v3 -> v4 (CCW)         v4 -> v3 -> v2 -> v1 (CW)
 *
 *   v4 *------<------* v3              v4 *------>------* v3
 *      |             |                    |             |
 *      v     CCW     ^     > 0            ^     CW      v     < 0
 *      |             |                    |             |
 *   v1 *------>------* v2              v1 *------<------* v2
 *
 * - Positive -> counter-clockwise (CCW) order
 * - Negative -> clockwise (CW) order
 * - Zero     -> degenerate (fewer than 3 vertices, or all points collinear)
 *
 * @param polygon Vertices of the polygon as an open ring (last vertex does not repeat first).
 * @return Signed area -> absolute value equals the geometric area of the polygon.
 * @note Self-intersecting polygons produce undefined results.
 * @see https://en.wikipedia.org/wiki/Shoelace_formula
 */
template <ScalarType T>
[[nodiscard]] inline T SignedPolygonArea(std::span<const Point2<T>> polygon) {
    // Degenerate case: less than 3 vertices -> area is zero
    if (polygon.size() < 3U) {
        return T{};  // zero area
    }

    // Shoelace -> each term is Cross(p0, p1) = twice the signed area of triangle (O, p0, p1)
    T sum{};
    const std::size_t num_vertices{polygon.size()};
    for (std::size_t i = 0; i < num_vertices; ++i) {
        const auto& p0 = polygon[i];
        const auto& p1 = polygon[(i + 1) % num_vertices];
        sum += p0.x * p1.y - p1.x * p0.y;
    }
    return static_cast<T>(0.5) * sum;
}

/// @brief Convenience overload for SignedPolygonArea -> accepts `std::vector<Point2<T>>` directly
template <ScalarType T>
[[nodiscard]] inline T SignedPolygonArea(const std::vector<Point2<T>>& polygon) {
    return SignedPolygonArea(std::span<const Point2<T>>{polygon});
}

/// @brief Area of a polygon
template <ScalarType T>
[[nodiscard]] inline T PolygonArea(std::span<const Point2<T>> polygon) {
    return std::abs(SignedPolygonArea(polygon));
}

/// @brief Convenience overload for PolygonArea -> accepts `std::vector<Point2<T>>` directly
template <ScalarType T>
[[nodiscard]] inline T PolygonArea(const std::vector<Point2<T>>& polygon) {
    return PolygonArea(std::span<const Point2<T>>{polygon});
}

}  // namespace geometry_kernel::core
