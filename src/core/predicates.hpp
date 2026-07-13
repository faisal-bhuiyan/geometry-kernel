#pragma once

#include "point.hpp"
#include "tolerance.hpp"
#include "types.hpp"

namespace geometry_kernel::core {

//---------------------------------------------------------------------------
// Triangle predicates
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

/**
 * @brief Orientation of ordered triple (a, b, c) in the plane.
 */
template <ScalarType T>
[[nodiscard]] inline Orientation TriangleOrientation(
    const Point<T>& a, const Point<T>& b, const Point<T>& c
) {
    const auto winding = RobustSign(SignedTriangleArea2(a, b, c));
    if (winding > 0) {
        return Orientation::kCounterClockwise;
    }
    if (winding < 0) {
        return Orientation::kClockwise;
    }
    return Orientation::kCollinear;
}

/**
 * @brief Returns true if previous -> current -> next makes a strict left turn (i.e. positive cross
 * product).
 * @pre Polygon vertices are in CCW order.
 */
template <ScalarType T>
[[nodiscard]] inline bool IsLeftTurn(
    const Point<T>& previous, const Point<T>& current, const Point<T>& next
) {
    return TriangleOrientation(previous, current, next) == Orientation::kCounterClockwise;
}

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
 */
template <ScalarType T>
[[nodiscard]] inline bool PointInTriangle(
    const Point<T>& point, const Point<T>& v1, const Point<T>& v2, const Point<T>& v3
) {
    const int s1{RobustSign(SignedTriangleArea2(v1, v2, point))};
    const int s2{RobustSign(SignedTriangleArea2(v2, v3, point))};
    const int s3{RobustSign(SignedTriangleArea2(v3, v1, point))};

    const bool has_negative{(s1 < 0) || (s2 < 0) || (s3 < 0)};
    const bool has_positive{(s1 > 0) || (s2 > 0) || (s3 > 0)};
    return !(has_negative && has_positive);
}

}  // namespace geometry_kernel::core
