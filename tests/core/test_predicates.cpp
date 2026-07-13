#include <gtest/gtest.h>

#include "core/predicates.hpp"
#include "test_fixtures.hpp"

namespace geometry_kernel::test {

using namespace geometry_kernel::core;

//------------------------------------------------------------------------------
// signed_triangle_area2 tests
//------------------------------------------------------------------------------

TEST(SignedTriangleArea2, UnitRightTriangleCCW) {
    // CCW -> positive area -> 2× area of 0.5 triangle = 1.0
    EXPECT_DOUBLE_EQ(SignedTriangleArea2(kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]), 1.);
}

TEST(SignedTriangleArea2, UnitRightTriangleCW) {
    // CW -> negative area -> 2× area of 0.5 triangle = -1.0
    EXPECT_DOUBLE_EQ(SignedTriangleArea2(kCwTriangle[0], kCwTriangle[1], kCwTriangle[2]), -1.);
}

TEST(SignedTriangleArea2, LargeCoordinates) {
    // Scales correctly -> 100×100 right triangle -> 2× area = 10000.0
    EXPECT_DOUBLE_EQ(
        SignedTriangleArea2(Point2D{0., 0.}, Point2D{100., 0.}, Point2D{0., 100.}), 10000.
    );
}

TEST(SignedTriangleArea2, kCollinearShapeHorizontal) {
    EXPECT_DOUBLE_EQ(SignedTriangleArea2(Point2D{0., 0.}, Point2D{1., 0.}, Point2D{2., 0.}), 0.);
}

TEST(SignedTriangleArea2, kCollinearShapeVertical) {
    EXPECT_DOUBLE_EQ(SignedTriangleArea2(Point2D{0., 0.}, Point2D{0., 1.}, Point2D{0., 2.}), 0.);
}

TEST(SignedTriangleArea2, AllCoincident) {
    EXPECT_DOUBLE_EQ(SignedTriangleArea2(Point2D{1., 1.}, Point2D{1., 1.}, Point2D{1., 1.}), 0.);
}

TEST(SignedTriangleArea2, AntisymmetryProperty) {
    // Swapping any two vertices negates the result
    const Point2D a{0., 0.}, b{3., 0.}, c{1., 2.};
    EXPECT_DOUBLE_EQ(SignedTriangleArea2(a, b, c), -SignedTriangleArea2(a, c, b));
    EXPECT_DOUBLE_EQ(SignedTriangleArea2(b, c, a), -SignedTriangleArea2(b, a, c));
    EXPECT_DOUBLE_EQ(SignedTriangleArea2(c, a, b), -SignedTriangleArea2(c, b, a));
}

//------------------------------------------------------------------------------
// triangle_orientation tests
//------------------------------------------------------------------------------

TEST(TriangleOrientation, CCWTriangle) {
    EXPECT_EQ(
        TriangleOrientation(kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]),
        Orientation::kCounterClockwise
    );
}

TEST(TriangleOrientation, CWTriangle) {
    EXPECT_EQ(
        TriangleOrientation(kCwTriangle[0], kCwTriangle[1], kCwTriangle[2]), Orientation::kClockwise
    );
}

TEST(TriangleOrientation, kCollinearShapePoints) {
    EXPECT_EQ(
        TriangleOrientation(Point2D{0., 0.}, Point2D{1., 0.}, Point2D{2., 0.}),
        Orientation::kCollinear
    );
}

TEST(TriangleOrientation, CoincidentPoints) {
    EXPECT_EQ(
        TriangleOrientation(Point2D{1., 1.}, Point2D{1., 1.}, Point2D{1., 1.}),
        Orientation::kCollinear
    );
}

TEST(TriangleOrientation, ReversingVertexOrderFlipsOrientation) {
    const Point2D a{0., 0.}, b{1., 0.}, c{0., 1.};
    EXPECT_NE(TriangleOrientation(a, b, c), TriangleOrientation(a, c, b));
    EXPECT_EQ(TriangleOrientation(a, b, c), Orientation::kCounterClockwise);
    EXPECT_EQ(TriangleOrientation(a, c, b), Orientation::kClockwise);
}

//------------------------------------------------------------------------------
// is_left_turn tests
//------------------------------------------------------------------------------

TEST(IsLeftTurn, ClearLeftTurn) {
    EXPECT_TRUE(IsLeftTurn(Point2D{0., 0.}, Point2D{1., 0.}, Point2D{1., 1.}));
}

TEST(IsLeftTurn, ClearRightTurn) {
    EXPECT_FALSE(IsLeftTurn(Point2D{0., 0.}, Point2D{1., 0.}, Point2D{1., -1.}));
}

TEST(IsLeftTurn, kCollinearShapeIsNotALeftTurn) {
    // Strict left turn -> kCollinearShape must return false
    EXPECT_FALSE(IsLeftTurn(Point2D{0., 0.}, Point2D{1., 0.}, Point2D{2., 0.}));
}

TEST(IsLeftTurn, UTurnIsNotALeftTurn) {
    EXPECT_FALSE(IsLeftTurn(Point2D{0., 0.}, Point2D{1., 0.}, Point2D{0., 0.}));
}

TEST(IsLeftTurn, ConvexSquareCornerCCW) {
    // Each corner of a CCW square is a left turn
    EXPECT_TRUE(IsLeftTurn(Point2D{0., 0.}, Point2D{1., 0.}, Point2D{1., 1.}));
    EXPECT_TRUE(IsLeftTurn(Point2D{1., 0.}, Point2D{1., 1.}, Point2D{0., 1.}));
    EXPECT_TRUE(IsLeftTurn(Point2D{1., 1.}, Point2D{0., 1.}, Point2D{0., 0.}));
    EXPECT_TRUE(IsLeftTurn(Point2D{0., 1.}, Point2D{0., 0.}, Point2D{1., 0.}));
}

TEST(IsLeftTurn, ReflexVertexIsNotALeftTurn) {
    // Reflex vertex of a CW polygon -> right turn
    EXPECT_FALSE(IsLeftTurn(Point2D{0., 0.}, Point2D{1., 0.}, Point2D{0., -1.}));
}

//------------------------------------------------------------------------------
// triangle_contains_point tests
//------------------------------------------------------------------------------

//--------------------------------
// Interior tests
//--------------------------------

TEST(TriangleContainsPoint, CentroidIsInside_CCW) {
    const Point2D centroid{1. / 3., 1. / 3.};  // centroid of CCW triangle
    EXPECT_TRUE(PointInTriangle(centroid, kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]));
}

TEST(TriangleContainsPoint, CentroidIsInside_CW) {
    // Winding-order agnostic -> same result for CW triangle
    const Point2D centroid{1. / 3., 1. / 3.};
    EXPECT_TRUE(PointInTriangle(centroid, kCwTriangle[0], kCwTriangle[1], kCwTriangle[2]));
}

TEST(TriangleContainsPoint, PointClearlyInside) {
    EXPECT_TRUE(
        PointInTriangle(Point2D{0.25, 0.25}, kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2])
    );
}

//--------------------------------
// Boundary vertices tests
//--------------------------------

TEST(TriangleContainsPoint, AtVertexV1) {
    EXPECT_TRUE(PointInTriangle(
        kCcwTriangle[0],  // vertex v1
        kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]
    ));
}

TEST(TriangleContainsPoint, AtVertexV2) {
    EXPECT_TRUE(PointInTriangle(
        kCcwTriangle[1],  // vertex v2
        kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]
    ));
}

TEST(TriangleContainsPoint, AtVertexV3) {
    EXPECT_TRUE(PointInTriangle(
        kCcwTriangle[2],  // vertex v3
        kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]
    ));
}

//--------------------------------
// Boundary edge midpoints tests
//--------------------------------

TEST(TriangleContainsPoint, MidpointOfEdgeV1V2) {
    EXPECT_TRUE(PointInTriangle(
        Point2D{0.5, 0.},  // midpoint of edge v1 -> v2
        kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]
    ));
}

TEST(TriangleContainsPoint, MidpointOfEdgeV2V3) {
    EXPECT_TRUE(PointInTriangle(
        Point2D{0.5, 0.5},  // midpoint of edge v2 -> v3
        kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]
    ));
}

TEST(TriangleContainsPoint, MidpointOfEdgeV3V1) {
    EXPECT_TRUE(PointInTriangle(
        Point2D{0., 0.5},  // midpoint of edge v3 -> v1
        kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]
    ));
}

//--------------------------------
// Outside tests
//--------------------------------

TEST(TriangleContainsPoint, ClearlyOutside_FarAway) {
    EXPECT_FALSE(PointInTriangle(
        Point2D{5., 5.},  // far away from triangle
        kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]
    ));
}

TEST(TriangleContainsPoint, JustBeyondEachEdge) {
    EXPECT_FALSE(PointInTriangle(
        Point2D{0.5, -0.01},  // just outside edge v1 -> v2 (below x-axis)
        kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]
    ));
    EXPECT_FALSE(PointInTriangle(
        Point2D{0.6, 0.6},  // just outside hypotenuse v2 -> v3
        kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]
    ));
    EXPECT_FALSE(PointInTriangle(
        Point2D{-0.01, 0.5},  // just outside edge v3 -> v1 (left of y-axis)
        kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]
    ));
}

//--------------------------------
// Winding-order agnosticism tests
//--------------------------------

TEST(TriangleContainsPoint, SameResultForCCWAndCW) {
    const Point2D p{0.25, 0.25};
    EXPECT_EQ(
        PointInTriangle(p, kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]),  // CCW
        PointInTriangle(p, kCwTriangle[0], kCwTriangle[1], kCwTriangle[2])      // CW
    );
}

//--------------------------------
// Degenerate triangle tests
//--------------------------------

TEST(TriangleContainsPoint, DegenerateAllCoincident_MatchingPoint) {
    // All three vertices at origin -> only the origin itself is "inside"
    EXPECT_TRUE(PointInTriangle(Point2D{0., 0.}, Point2D{0., 0.}, Point2D{0., 0.}, Point2D{0., 0.}));
}

TEST(TriangleContainsPoint, DegenerateAllCoincident_NonMatchingPoint) {
    // Known limitation: when all three vertices are coincident, all signed areas
    // are zero, so the mixed-sign test never fires. Every point returns true.
    // The caller is responsible for filtering degenerate triangles upstream.
    EXPECT_TRUE(
        PointInTriangle(Point2D{1., 0.}, Point2D{0., 0.}, Point2D{0., 0.}, Point2D{0., 0.})
    );  // degeneracy: all-zero signs -> always inside
}

TEST(TriangleContainsPoint, DegeneratekCollinearShape_PointOffLine) {
    EXPECT_FALSE(
        PointInTriangle(Point2D{1., 1.}, Point2D{0., 0.}, Point2D{2., 0.}, Point2D{4., 0.})
    );
}

}  // namespace geometry_kernel::test
