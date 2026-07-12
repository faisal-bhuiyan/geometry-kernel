#pragma once

#include <span>
#include <vector>

#include "point.hpp"
#include "predicates.hpp"
#include "types.hpp"

namespace geometry_kernel::core {

//---------------------------------------------------------------------------
// Triangle
//---------------------------------------------------------------------------

/**
 * @brief Signed 2D cross product: z-component of vector (b - a) x (c - a).
 *
 * @return Signed parallelogram area in the xy-plane; sign follows the right-hand rule.
 *         Equals twice the signed area of triangle (a, b, c).
 */
template <ScalarType T>
[[nodiscard]] inline T SignedTriangleArea2(const Point<T>& a, const Point<T>& b, const Point<T>& c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

/// @brief Area of triangle (a, b, c).
template <ScalarType T>
[[nodiscard]] inline T TriangleArea(const Point<T>& a, const Point<T>& b, const Point<T>& c) {
    return static_cast<T>(0.5) * std::abs(SignedTriangleArea2(a, b, c));
}

//---------------------------------------------------------------------------
// Polygon
//---------------------------------------------------------------------------

/**
 * @brief Computes the signed area of a simple polygon using the shoelace formula.
 *
 * The sign encodes winding order:
 * - Positive -> vertices are in counter-clockwise (CCW) order
 * - Negative -> vertices are in clockwise (CW) order
 * - Zero     -> polygon is degenerate (fewer than 3 vertices, or all points collinear)
 *
 * @param polygon Vertices of the polygon as an open ring (last vertex does not repeat first).
 * @return Signed area - absolute value equals the geometric area of the polygon.
 * @note Self-intersecting polygons produce undefined results.
 * @see https://en.wikipedia.org/wiki/Shoelace_formula
 */
template <ScalarType T>
[[nodiscard]] inline T SignedPolygonArea(std::span<const Point<T>> polygon) {
    // Degenerate case: less than 3 vertices -> area is zero
    if (polygon.size() < 3U) {
        return T{};  // zero area
    }

    // Shoelace formula: https://en.wikipedia.org/wiki/Shoelace_formula
    T sum{};
    const std::size_t num_vertices{polygon.size()};
    for (std::size_t i = 0; i < num_vertices; ++i) {
        const auto& p0 = polygon[i];
        const auto& p1 = polygon[(i + 1) % num_vertices];
        sum += p0.x * p1.y - p1.x * p0.y;
    }
    return static_cast<T>(0.5) * sum;
}

/// @brief Convenience overload for SignedPolygonArea -> accepts `std::vector<Point<T>>` directly
template <ScalarType T>
[[nodiscard]] inline T SignedPolygonArea(const std::vector<Point<T>>& polygon) {
    return SignedPolygonArea(std::span<const Point<T>>{polygon});
}

/// @brief Area of a polygon.
template <ScalarType T>
[[nodiscard]] inline T PolygonArea(std::span<const Point<T>> polygon) {
    return std::abs(SignedPolygonArea(polygon));
}

/// @brief Convenience overload for PolygonArea -> accepts `std::vector<Point<T>>` directly
template <ScalarType T>
[[nodiscard]] inline T PolygonArea(const std::vector<Point<T>>& polygon) {
    return PolygonArea(std::span<const Point<T>>{polygon});
}

}  // namespace geometry_kernel::core
