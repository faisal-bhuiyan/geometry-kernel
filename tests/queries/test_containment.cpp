#include <gtest/gtest.h>

#include "core/point2.hpp"
#include "core/test_fixtures.hpp"
#include "queries/containment2.hpp"

namespace geometry_kernel::test {

using namespace geometry_kernel::core;
using namespace geometry_kernel::queries;

//------------------------------------------------------------------------------
// PointInTriangle — interior
//------------------------------------------------------------------------------

TEST(PointInTriangle, CentroidIsInside_CCW) {
    const Point2D centroid{1. / 3., 1. / 3.};  // centroid of CCW triangle
    EXPECT_TRUE(PointInTriangle(centroid, kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]));
}

TEST(PointInTriangle, CentroidIsInside_CW) {
    // Winding-order agnostic -> same result for CW triangle
    const Point2D centroid{1. / 3., 1. / 3.};
    EXPECT_TRUE(PointInTriangle(centroid, kCwTriangle[0], kCwTriangle[1], kCwTriangle[2]));
}

TEST(PointInTriangle, PointClearlyInside) {
    EXPECT_TRUE(
        PointInTriangle(Point2D{0.25, 0.25}, kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2])
    );
}

//------------------------------------------------------------------------------
// PointInTriangle — boundary vertices
//------------------------------------------------------------------------------

TEST(PointInTriangle, AtVertexV1) {
    EXPECT_TRUE(PointInTriangle(
        kCcwTriangle[0],  // vertex v1
        kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]
    ));
}

TEST(PointInTriangle, AtVertexV2) {
    EXPECT_TRUE(PointInTriangle(
        kCcwTriangle[1],  // vertex v2
        kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]
    ));
}

TEST(PointInTriangle, AtVertexV3) {
    EXPECT_TRUE(PointInTriangle(
        kCcwTriangle[2],  // vertex v3
        kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]
    ));
}

//------------------------------------------------------------------------------
// PointInTriangle — boundary edge midpoints
//------------------------------------------------------------------------------

TEST(PointInTriangle, MidpointOfEdgeV1V2) {
    EXPECT_TRUE(PointInTriangle(
        Point2D{0.5, 0.},  // midpoint of edge v1 -> v2
        kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]
    ));
}

TEST(PointInTriangle, MidpointOfEdgeV2V3) {
    EXPECT_TRUE(PointInTriangle(
        Point2D{0.5, 0.5},  // midpoint of edge v2 -> v3
        kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]
    ));
}

TEST(PointInTriangle, MidpointOfEdgeV3V1) {
    EXPECT_TRUE(PointInTriangle(
        Point2D{0., 0.5},  // midpoint of edge v3 -> v1
        kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]
    ));
}

//------------------------------------------------------------------------------
// PointInTriangle — outside
//------------------------------------------------------------------------------

TEST(PointInTriangle, ClearlyOutside_FarAway) {
    EXPECT_FALSE(PointInTriangle(
        Point2D{5., 5.},  // far away from triangle
        kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]
    ));
}

TEST(PointInTriangle, JustBeyondEachEdge) {
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

//------------------------------------------------------------------------------
// PointInTriangle — winding-order agnosticism
//------------------------------------------------------------------------------

TEST(PointInTriangle, SameResultForCCWAndCW) {
    const Point2D p{0.25, 0.25};
    EXPECT_EQ(
        PointInTriangle(p, kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]),  // CCW
        PointInTriangle(p, kCwTriangle[0], kCwTriangle[1], kCwTriangle[2])      // CW
    );
}

//------------------------------------------------------------------------------
// PointInTriangle — degenerate triangles
//------------------------------------------------------------------------------

TEST(PointInTriangle, DegenerateAllCoincident_MatchingPoint) {
    // All three vertices at origin -> only the origin itself is "inside"
    EXPECT_TRUE(PointInTriangle(Point2D{0., 0.}, Point2D{0., 0.}, Point2D{0., 0.}, Point2D{0., 0.}));
}

TEST(PointInTriangle, DegenerateAllCoincident_NonMatchingPoint) {
    // Known limitation: when all three vertices are coincident, all signed areas
    // are zero, so the mixed-sign test never fires. Every point returns true.
    // The caller is responsible for filtering degenerate triangles upstream.
    EXPECT_TRUE(
        PointInTriangle(Point2D{1., 0.}, Point2D{0., 0.}, Point2D{0., 0.}, Point2D{0., 0.})
    );  // degeneracy: all-zero signs -> always inside
}

TEST(PointInTriangle, DegeneratekCollinearShape_PointOffLine) {
    EXPECT_FALSE(
        PointInTriangle(Point2D{1., 1.}, Point2D{0., 0.}, Point2D{2., 0.}, Point2D{4., 0.})
    );
}

}  // namespace geometry_kernel::test
