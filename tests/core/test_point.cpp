#include <limits>

#include <gtest/gtest.h>

#include "core/point.hpp"
#include "core/tolerance.hpp"

namespace geometry_kernel::test {

using namespace geometry_kernel::core;

//------------------------------------------------------------------------------
// Point2 construction
//------------------------------------------------------------------------------

TEST(Point2, DefaultInitIsOrigin) {
    const Point2D p{};
    EXPECT_DOUBLE_EQ(p.x, 0.);
    EXPECT_DOUBLE_EQ(p.y, 0.);
}

TEST(Point2, AggregateInitStoresCoordinates) {
    const Point2D p{1.5, -2.25};
    EXPECT_DOUBLE_EQ(p.x, 1.5);
    EXPECT_DOUBLE_EQ(p.y, -2.25);
}

//------------------------------------------------------------------------------
// NearlyEqualPoint
//------------------------------------------------------------------------------

TEST(NearlyEqualPoint, IdenticalPoints) {
    const Point2D a{1., 2.};
    EXPECT_TRUE(NearlyEqualPoint(a, a));
}

TEST(NearlyEqualPoint, JustInsideTolerance) {
    const Point2D a{0., 0.};
    const double delta = kTolerance<double> * 0.5;
    EXPECT_TRUE(NearlyEqualPoint(a, Point2D{delta, 0.}));
    EXPECT_TRUE(NearlyEqualPoint(a, Point2D{0., delta}));
}

TEST(NearlyEqualPoint, ExactlyAtToleranceBoundary) {
    // Comparison is <=, so exactly kTolerance is still equal
    const Point2D a{0., 0.};
    EXPECT_TRUE(NearlyEqualPoint(a, Point2D{kTolerance<double>, 0.}));
    EXPECT_TRUE(NearlyEqualPoint(a, Point2D{0., kTolerance<double>}));
    EXPECT_TRUE(NearlyEqualPoint(a, Point2D{-kTolerance<double>, 0.}));
}

TEST(NearlyEqualPoint, JustOutsideTolerance) {
    const Point2D a{0., 0.};
    const double eps = std::numeric_limits<double>::epsilon();
    EXPECT_FALSE(NearlyEqualPoint(a, Point2D{kTolerance<double> + eps, 0.}));
    EXPECT_FALSE(NearlyEqualPoint(a, Point2D{0., kTolerance<double> + eps}));
}

TEST(NearlyEqualPoint, OneAxisOutsideFails) {
    // Both axes must be within tolerance
    const Point2D a{0., 0.};
    const double outside = kTolerance<double> * 2.;
    EXPECT_FALSE(NearlyEqualPoint(a, Point2D{0., outside}));
    EXPECT_FALSE(NearlyEqualPoint(a, Point2D{outside, 0.}));
    EXPECT_FALSE(NearlyEqualPoint(a, Point2D{outside, outside}));
}

TEST(NearlyEqualPoint, NegativeCoordinates) {
    // Avoid adding kTolerance to large magnitudes (rounding can push the delta
    // just outside <= kTolerance). Use a clearly-inside half-tolerance offset.
    const Point2D a{-10., -20.};
    const double half = kTolerance<double> * 0.5;
    EXPECT_TRUE(NearlyEqualPoint(a, Point2D{-10., -20.}));
    EXPECT_TRUE(NearlyEqualPoint(a, Point2D{-10. + half, -20. - half}));
    EXPECT_FALSE(NearlyEqualPoint(a, Point2D{-10. + 1., -20.}));
}

TEST(NearlyEqualPoint, LargeCoordinates) {
    const Point2D a{1e6, -1e6};
    const double half = kTolerance<double> * 0.5;
    EXPECT_TRUE(NearlyEqualPoint(a, Point2D{1e6, -1e6}));
    EXPECT_TRUE(NearlyEqualPoint(a, Point2D{1e6 + half, -1e6 - half}));
    EXPECT_FALSE(NearlyEqualPoint(a, Point2D{1e6 + 1., -1e6}));
}

TEST(NearlyEqualPoint, Symmetric) {
    const Point2D a{1., 2.};
    const Point2D b{1. + kTolerance<double> * 0.5, 2.};
    EXPECT_EQ(NearlyEqualPoint(a, b), NearlyEqualPoint(b, a));

    const Point2D c{1. + kTolerance<double> * 2., 2.};
    EXPECT_EQ(NearlyEqualPoint(a, c), NearlyEqualPoint(c, a));
}

}  // namespace geometry_kernel::test
