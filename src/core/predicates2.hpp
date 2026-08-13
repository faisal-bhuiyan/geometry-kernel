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
 *
 * Classifies the winding of (a, b, c) via RobustSign(SignedTriangleAreaTimes2).
 * The doubled signed area is positive exactly when the walk a -> b -> c turns
 * counter-clockwise, so an area computation doubles as a winding test.
 *
 * Case 1: CCW -> kCounterClockwise
 *
 *               c
 *               *
 *              / \
 *             /   \
 *            v     ^
 *           /  CCW  \
 *          /         \
 *         *----->-----*
 *         a           b
 *
 * Case 2: CW -> kClockwise
 *
 *               b
 *               *
 *              / \
 *             /   \
 *            ^     v
 *           /  CW   \
 *          /         \
 *         *-----<-----*
 *         a           c
 *
 * Case 3: collinear -> kCollinear
 *
 *         *-----------*-----------*
 *         a           b           c
 *
 *         no enclosed area -> no winding to report
 *
 * The same predicate read as a side-of-line test: the cases differ only in
 * which side of the directed line a -> b the third point falls on.
 *
 *                      * c            c left of  a -> b  ->  kCounterClockwise
 *
 *      a *---------------------> b    c on       a -> b  ->  kCollinear
 *
 *                      * c            c right of a -> b  ->  kClockwise
 *
 * That equivalence is what lets the containment and intersection queries reuse
 * this as "is the point on the inner side of that edge".
 *
 * Two invariants follow from the underlying determinant:
 * - Cyclic rotation preserves it: (a, b, c), (b, c, a), (c, a, b) all agree.
 * - Any single swap reverses it:  (a, b, c) and (a, c, b) always disagree.
 *
 * @param a First point of the ordered triple.
 * @param b Second point of the ordered triple.
 * @param c Third point of the ordered triple.
 * @return Orientation enum value for the winding of (a, b, c).
 *
 * @note Classification runs through RobustSign, so a triple whose doubled area
 *       lands inside +/-kTolerance reports kCollinear rather than picking up an
 *       arbitrary winding from floating-point noise.
 *
 * @see de Berg, M. et al. Computational Geometry: Algorithms and Applications,
 *      3rd ed., Springer, 2008, §1.1.
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
 * @brief Returns true if previous -> current -> next makes a strict left turn
 *        (i.e. positive cross product / CCW orientation).
 *
 * A left turn at current means next lies strictly to the left of the directed
 * line previous -> current. It is TriangleOrientation phrased as a walk along a
 * path rather than as a triangle's winding.
 *
 * Case 1: left turn -> true
 *
 *                               * next
 *                              ^
 *                             /
 *                            /
 *      prev *------->------* curr
 *
 *      heading right along prev -> curr, then swinging CCW to reach next
 *
 * Case 2: collinear -> false
 *
 *      prev *------->------* curr ------->------* next
 *
 *      straight on -> no turn, and a strict test rejects it
 *
 * Case 3: right turn -> false
 *
 *      prev *------->------* curr
 *                           \
 *                            \
 *                             v
 *                              * next
 *
 *      swinging CW at curr
 *
 * Why the winding of the input polygon matters: in a CCW-wound polygon the
 * interior lies to the left of every directed edge, so "left turn at current"
 * is exactly "the interior angle at current is below 180 deg" -- that is,
 * current is a convex vertex. Hand it CW input and every answer inverts:
 * convex corners report false and reflex corners report true. That is why
 * TriangulatePolygonWithEarClipping reverses CW input before constructing its
 * PolygonMesh, whose IsConvexVertex is this function's only caller.
 *
 * Strictness matters for the same consumer: collinear corners are never
 * classified convex, so they are never clipped as ears and no zero-area
 * triangle reaches the output.
 *
 * @param previous Previous vertex on the path.
 * @param current Current vertex (the corner being tested).
 * @param next Next vertex on the path.
 * @return True iff (previous, current, next) is strictly CCW.
 *
 * @note The predicate itself is total and valid for any three points. The
 *       convex-vertex reading above is what requires CCW winding.
 *
 * @see TriangleOrientation
 */
template <ScalarType T>
[[nodiscard]] inline bool IsLeftTurn(
    const Point2<T>& previous, const Point2<T>& current, const Point2<T>& next
) {
    return TriangleOrientation(previous, current, next) == Orientation::kCounterClockwise;
}

}  // namespace geometry_kernel::core
