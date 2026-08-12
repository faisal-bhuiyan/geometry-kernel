#include <gtest/gtest.h>

#include "core/point2.hpp"
#include "queries/intersection.hpp"

namespace geometry_kernel::test {

using namespace geometry_kernel::core;
using namespace geometry_kernel::queries;

//------------------------------------------------------------------------------
// PointOnSegment — documented cases
//------------------------------------------------------------------------------

TEST(PointOnSegment, BetweenEndpoints) {
    EXPECT_TRUE(PointOnSegment(Point2D{0.5, 0.}, Point2D{0., 0.}, Point2D{1., 0.}));
}

TEST(PointOnSegment, AtEndpointA) {
    EXPECT_TRUE(PointOnSegment(Point2D{0., 0.}, Point2D{0., 0.}, Point2D{1., 0.}));
}

TEST(PointOnSegment, AtEndpointB) {
    EXPECT_TRUE(PointOnSegment(Point2D{1., 0.}, Point2D{0., 0.}, Point2D{1., 0.}));
}

TEST(PointOnSegment, CollinearButOutside) {
    // On the infinite line through [a,b], but outside the segment extent
    EXPECT_FALSE(PointOnSegment(Point2D{-1., 0.}, Point2D{0., 0.}, Point2D{1., 0.}));
    EXPECT_FALSE(PointOnSegment(Point2D{2., 0.}, Point2D{0., 0.}, Point2D{1., 0.}));
}

TEST(PointOnSegment, OrderIndependent) {
    const Point2D a{0., 0.}, b{2., 2.}, p{1., 1.};
    EXPECT_EQ(PointOnSegment(p, a, b), PointOnSegment(p, b, a));
}

TEST(PointOnSegment, DiagonalSegment) {
    EXPECT_TRUE(PointOnSegment(Point2D{1., 1.}, Point2D{0., 0.}, Point2D{2., 2.}));
    EXPECT_FALSE(PointOnSegment(Point2D{3., 3.}, Point2D{0., 0.}, Point2D{2., 2.}));
}

//------------------------------------------------------------------------------
// SegmentsIntersect — documented cases
//------------------------------------------------------------------------------

TEST(SegmentsIntersect, ProperCrossing) {
    // a-----b horizontal, c-----d vertical through the midpoint
    const Point2D a{0., 0.}, b{2., 0.}, c{1., -1.}, d{1., 1.};
    EXPECT_TRUE(SegmentsIntersect(a, b, c, d));
}

TEST(SegmentsIntersect, EndpointOfOneLiesOnOther) {
    // c lies on [a,b]
    const Point2D a{0., 0.}, b{2., 0.}, c{1., 0.}, d{1., 1.};
    EXPECT_TRUE(SegmentsIntersect(a, b, c, d));
}

TEST(SegmentsIntersect, SharedEndpoint) {
    // b == c
    const Point2D a{0., 0.}, b{1., 0.}, c{1., 0.}, d{1., 1.};
    EXPECT_TRUE(SegmentsIntersect(a, b, c, d));
}

TEST(SegmentsIntersect, CollinearButDisjoint) {
    const Point2D a{0., 0.}, b{1., 0.}, c{2., 0.}, d{3., 0.};
    EXPECT_FALSE(SegmentsIntersect(a, b, c, d));
}

TEST(SegmentsIntersect, SameSideNonCrossing) {
    // Both c and d above ab -> no crossing
    const Point2D a{0., 0.}, b{2., 0.}, c{0.5, 1.}, d{1.5, 1.};
    EXPECT_FALSE(SegmentsIntersect(a, b, c, d));
}

TEST(SegmentsIntersect, SymmetricInSegmentOrder) {
    const Point2D a{0., 0.}, b{2., 0.}, c{1., -1.}, d{1., 1.};
    EXPECT_EQ(SegmentsIntersect(a, b, c, d), SegmentsIntersect(c, d, a, b));
}

TEST(SegmentsIntersect, EndpointOrderIndependent) {
    const Point2D a{0., 0.}, b{2., 0.}, c{1., -1.}, d{1., 1.};
    EXPECT_EQ(SegmentsIntersect(a, b, c, d), SegmentsIntersect(b, a, c, d));
    EXPECT_EQ(SegmentsIntersect(a, b, c, d), SegmentsIntersect(a, b, d, c));
}

TEST(SegmentsIntersect, CollinearOverlapping) {
    // Partially overlapping collinear segments share points -> true
    const Point2D a{0., 0.}, b{2., 0.}, c{1., 0.}, d{3., 0.};
    EXPECT_TRUE(SegmentsIntersect(a, b, c, d));
}

TEST(SegmentsIntersect, TJunction) {
    // Vertical segment ending on the horizontal one (endpoint on other)
    const Point2D a{0., 0.}, b{2., 0.}, c{1., 0.}, d{1., 2.};
    EXPECT_TRUE(SegmentsIntersect(a, b, c, d));
}

}  // namespace geometry_kernel::test
