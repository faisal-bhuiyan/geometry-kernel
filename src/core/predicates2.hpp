#pragma once

#include "area.hpp"
#include "point2.hpp"
#include "tolerance.hpp"
#include "types.hpp"
#include "vector2.hpp"

namespace geometry_kernel::core {

//---------------------------------------------------------------------------
// Triangle predicates (2D)
//---------------------------------------------------------------------------

/**
 * @brief Orientation of ordered triple (a, b, c) in the plane.
 */
template <ScalarType T>
[[nodiscard]] inline Orientation TriangleOrientation(
    const Point2<T>& a, const Point2<T>& b, const Point2<T>& c
) {
    const auto winding = RobustSign(SignedTriangleAreaTimes2(a, b, c));
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
    const Point2<T>& previous, const Point2<T>& current, const Point2<T>& next
) {
    return TriangleOrientation(previous, current, next) == Orientation::kCounterClockwise;
}

}  // namespace geometry_kernel::core
