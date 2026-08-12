#include <limits>

#include <gtest/gtest.h>

#include "core/point3.hpp"
#include "core/tolerance.hpp"

namespace geometry_kernel::test {

using namespace geometry_kernel::core;

//------------------------------------------------------------------------------
// Point3 construction
//------------------------------------------------------------------------------

TEST(Point3, DefaultInitIsOrigin) {
    const Point3D p{};
    EXPECT_DOUBLE_EQ(p.x, 0.);
    EXPECT_DOUBLE_EQ(p.y, 0.);
    EXPECT_DOUBLE_EQ(p.z, 0.);
}

TEST(Point3, AggregateInitStoresCoordinates) {
    const Point3D p{1.5, -2.25, 3.5};
    EXPECT_DOUBLE_EQ(p.x, 1.5);
    EXPECT_DOUBLE_EQ(p.y, -2.25);
    EXPECT_DOUBLE_EQ(p.z, 3.5);
}

//------------------------------------------------------------------------------
// NearlyEqualPoint (Point3)
//------------------------------------------------------------------------------

TEST(NearlyEqualPoint3, IdenticalPoints) {
    const Point3D a{1., 2., 3.};
    EXPECT_TRUE(NearlyEqualPoint(a, a));
}

TEST(NearlyEqualPoint3, JustInsideTolerance) {
    const Point3D a{0., 0., 0.};
    const double delta = kTolerance<double> * 0.5;
    EXPECT_TRUE(NearlyEqualPoint(a, Point3D{delta, 0., 0.}));
    EXPECT_TRUE(NearlyEqualPoint(a, Point3D{0., delta, 0.}));
    EXPECT_TRUE(NearlyEqualPoint(a, Point3D{0., 0., delta}));
}

TEST(NearlyEqualPoint3, ExactlyAtToleranceBoundary) {
    // Comparison is <=, so exactly kTolerance is still equal
    const Point3D a{0., 0., 0.};
    EXPECT_TRUE(NearlyEqualPoint(a, Point3D{kTolerance<double>, 0., 0.}));
    EXPECT_TRUE(NearlyEqualPoint(a, Point3D{0., kTolerance<double>, 0.}));
    EXPECT_TRUE(NearlyEqualPoint(a, Point3D{0., 0., kTolerance<double>}));
    EXPECT_TRUE(NearlyEqualPoint(a, Point3D{-kTolerance<double>, 0., 0.}));
}

TEST(NearlyEqualPoint3, JustOutsideTolerance) {
    const Point3D a{0., 0., 0.};
    const double eps = std::numeric_limits<double>::epsilon();
    EXPECT_FALSE(NearlyEqualPoint(a, Point3D{kTolerance<double> + eps, 0., 0.}));
    EXPECT_FALSE(NearlyEqualPoint(a, Point3D{0., kTolerance<double> + eps, 0.}));
    EXPECT_FALSE(NearlyEqualPoint(a, Point3D{0., 0., kTolerance<double> + eps}));
}

TEST(NearlyEqualPoint3, OneAxisOutsideFails) {
    // All three axes must be within tolerance
    const Point3D a{0., 0., 0.};
    const double outside = kTolerance<double> * 2.;
    EXPECT_FALSE(NearlyEqualPoint(a, Point3D{0., 0., outside}));
    EXPECT_FALSE(NearlyEqualPoint(a, Point3D{0., outside, 0.}));
    EXPECT_FALSE(NearlyEqualPoint(a, Point3D{outside, 0., 0.}));
    EXPECT_FALSE(NearlyEqualPoint(a, Point3D{outside, outside, outside}));
}

TEST(NearlyEqualPoint3, NegativeCoordinates) {
    const Point3D a{-10., -20., -30.};
    const double half = kTolerance<double> * 0.5;
    EXPECT_TRUE(NearlyEqualPoint(a, Point3D{-10., -20., -30.}));
    EXPECT_TRUE(NearlyEqualPoint(a, Point3D{-10. + half, -20. - half, -30. + half}));
    EXPECT_FALSE(NearlyEqualPoint(a, Point3D{-10. + 1., -20., -30.}));
}

TEST(NearlyEqualPoint3, LargeCoordinates) {
    const Point3D a{1e6, -1e6, 1e6};
    const double half = kTolerance<double> * 0.5;
    EXPECT_TRUE(NearlyEqualPoint(a, Point3D{1e6, -1e6, 1e6}));
    EXPECT_TRUE(NearlyEqualPoint(a, Point3D{1e6 + half, -1e6 - half, 1e6 + half}));
    EXPECT_FALSE(NearlyEqualPoint(a, Point3D{1e6 + 1., -1e6, 1e6}));
}

TEST(NearlyEqualPoint3, Symmetric) {
    const Point3D a{1., 2., 3.};
    const Point3D b{1. + kTolerance<double> * 0.5, 2., 3.};
    EXPECT_EQ(NearlyEqualPoint(a, b), NearlyEqualPoint(b, a));

    const Point3D c{1. + kTolerance<double> * 2., 2., 3.};
    EXPECT_EQ(NearlyEqualPoint(a, c), NearlyEqualPoint(c, a));
}

}  // namespace geometry_kernel::test
