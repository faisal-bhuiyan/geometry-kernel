#include <gtest/gtest.h>

#include "core/predicates2.hpp"
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

}  // namespace geometry_kernel::test
