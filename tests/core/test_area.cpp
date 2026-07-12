#include <gtest/gtest.h>

#include "core/area.hpp"
#include "core/point.hpp"
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
// triangle_area tests
//------------------------------------------------------------------------------

TEST(TriangleArea, UnitRightTriangle) {
    EXPECT_DOUBLE_EQ(TriangleArea(kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]), 0.5);
}

TEST(TriangleArea, AlwaysNonNegativeForCW) {
    // CW winding -> area must still be positive
    EXPECT_DOUBLE_EQ(TriangleArea(kCwTriangle[0], kCwTriangle[1], kCwTriangle[2]), 0.5);
}

TEST(TriangleArea, CCWAndCWGiveSameArea) {
    EXPECT_DOUBLE_EQ(
        TriangleArea(kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]),  // CCW
        TriangleArea(kCwTriangle[0], kCwTriangle[1], kCwTriangle[2])      // CW
    );
}

TEST(TriangleArea, ThreeFourFiveTriangle) {
    // Right triangle with legs 3 and 4 -> area = 6.0
    EXPECT_DOUBLE_EQ(TriangleArea(Point2D{0., 0.}, Point2D{3., 0.}, Point2D{0., 4.}), 6.);
}

TEST(TriangleArea, DegeneratekCollinearShape) {
    EXPECT_DOUBLE_EQ(TriangleArea(Point2D{0., 0.}, Point2D{1., 0.}, Point2D{2., 0.}), 0.);
}

TEST(TriangleArea, ScaleByFactorOfTwo) {
    // Doubling coordinates scales area by 4
    EXPECT_DOUBLE_EQ(
        TriangleArea(Point2D{0., 0.}, Point2D{2., 0.}, Point2D{0., 2.}),
        4. * TriangleArea(Point2D{0., 0.}, Point2D{1., 0.}, Point2D{0., 1.})
    );
}

//------------------------------------------------------------------------------
// polygon_area_signed tests
//------------------------------------------------------------------------------

//--------------------------------
// Known areas and sign tests
//--------------------------------

TEST(SignedPolygonArea, CCWUnitSquare_PositiveArea) {
    EXPECT_DOUBLE_EQ(SignedPolygonArea(kCcwSquare), 1.);
}

TEST(SignedPolygonArea, CWUnitSquare_NegativeArea) {
    EXPECT_DOUBLE_EQ(SignedPolygonArea(kCwSquare), -1.);
}

TEST(SignedPolygonArea, CCWRectangle) {
    const std::vector<Point2D> rectangle = {{0., 0.}, {3., 0.}, {3., 4.}, {0., 4.}};
    EXPECT_DOUBLE_EQ(SignedPolygonArea(rectangle), 12.);
}

TEST(SignedPolygonArea, ReversedVertexOrderNegatesResult) {
    std::vector<Point2D> reversed(kCcwSquare.rbegin(), kCcwSquare.rend());
    EXPECT_DOUBLE_EQ(SignedPolygonArea(reversed), -SignedPolygonArea(kCcwSquare));
}

//--------------------------------
// Concave polygons tests
//--------------------------------

TEST(SignedPolygonArea, LShapeChallengePolygon) {
    EXPECT_NEAR(SignedPolygonArea(kSimpleConcavePolygon), -0.65, 1e-10);
}

//--------------------------------
// Degenerate inputs tests
//--------------------------------

TEST(SignedPolygonArea, EmptySpan) {
    EXPECT_DOUBLE_EQ(SignedPolygonArea(std::span<const Point2D>{}), 0.);
}

TEST(SignedPolygonArea, ThreekCollinearShapeVertices) {
    const std::vector<Point2D> kCollinearShape = {{0., 0.}, {1., 0.}, {2., 0.}};
    EXPECT_DOUBLE_EQ(SignedPolygonArea(kCollinearShape), 0.);
}

//------------------------------------------------------------------------------
// polygon_area tests
//------------------------------------------------------------------------------

TEST(PolygonArea, CCWSquare) {
    EXPECT_DOUBLE_EQ(PolygonArea(kCcwSquare), 1.);
}

TEST(PolygonArea, CWSquareSameAsCCW) {
    // Absolute value: winding order does not affect result
    EXPECT_DOUBLE_EQ(PolygonArea(kCwSquare), PolygonArea(kCcwSquare));
}

TEST(PolygonArea, LShapeChallengePolygon) {
    EXPECT_NEAR(PolygonArea(kSimpleConcavePolygon), 0.65, 1e-10);
}

TEST(PolygonArea, EmptySpan) {
    EXPECT_DOUBLE_EQ(PolygonArea(std::span<const Point2D>{}), 0.);
}

TEST(PolygonArea, AlwaysNonNegative) {
    EXPECT_GT(PolygonArea(kCcwSquare), 0.);
    EXPECT_GT(PolygonArea(kCwSquare), 0.);
}

//------------------------------------------------------------------------------
// Cross-cutting properties tests
//------------------------------------------------------------------------------

TEST(Properties, TriangleAreaMatchesPolygonAreaForThreeVertices) {
    // A triangle as a polygon should give the same area as triangle_area
    const std::vector<Point2D> triangle(kCcwTriangle.begin(), kCcwTriangle.end());
    EXPECT_DOUBLE_EQ(
        PolygonArea(triangle), TriangleArea(kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2])
    );
}

TEST(Properties, SquareAreaEqualsSumOfTwoTriangles) {
    // Manually triangulate the CCW unit square and verify areas sum correctly
    const double t1 = TriangleArea(Point2D{0., 0.}, Point2D{1., 0.}, Point2D{1., 1.});
    const double t2 = TriangleArea(Point2D{0., 0.}, Point2D{1., 1.}, Point2D{0., 1.});
    EXPECT_DOUBLE_EQ(t1 + t2, PolygonArea(kCcwSquare));
}

TEST(Properties, IsLeftTurnConsistentWithTriangleOrientation) {
    // is_left_turn must agree with triangle_orientation == kCounterClockwise
    const Point2D a{0., 0.}, b{1., 0.}, c{1., 1.};
    EXPECT_EQ(IsLeftTurn(a, b, c), TriangleOrientation(a, b, c) == Orientation::kCounterClockwise);
}

TEST(Properties, SignedArea2AntisymmetryUnderVertexSwap) {
    const Point2D a{0., 0.}, b{3., 0.}, c{1., 2.};
    EXPECT_DOUBLE_EQ(SignedTriangleArea2(a, b, c), -SignedTriangleArea2(a, c, b));
}

TEST(Properties, CentroidOfAnyCCWTriangleIsContained) {
    // Centroid = average of vertices; always strictly interior
    const Point2D a{0., 0.}, b{4., 0.}, c{2., 3.};
    const Point2D centroid{(a.x + b.x + c.x) / 3., (a.y + b.y + c.y) / 3.};
    EXPECT_TRUE(PointInTriangle(centroid, a, b, c));
}

TEST(Properties, ScaleInvarianceOfPolygonArea) {
    // Scaling all coordinates by k scales area by k^2
    const double k = 3.;
    std::vector<Point2D> scaled;
    for (const auto& p : kCcwSquare)
        scaled.push_back(Point2D{p.x * k, p.y * k});
    EXPECT_NEAR(PolygonArea(scaled), k * k * PolygonArea(kCcwSquare), 1e-10);
}

TEST(Properties, SignedPolygonAreaAndPolygonAreaRelationship) {
    // polygon_area == |polygon_area_signed| for both windings
    EXPECT_DOUBLE_EQ(PolygonArea(kCcwSquare), std::abs(SignedPolygonArea(kCcwSquare)));
    EXPECT_DOUBLE_EQ(PolygonArea(kCwSquare), std::abs(SignedPolygonArea(kCwSquare)));
}

}  // namespace geometry_kernel::test
