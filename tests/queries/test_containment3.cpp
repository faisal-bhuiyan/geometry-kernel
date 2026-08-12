#include <gtest/gtest.h>

#include "core/point3.hpp"
#include "core/vector3.hpp"
#include "queries/containment3.hpp"

namespace geometry_kernel::test {

using namespace geometry_kernel::core;
using namespace geometry_kernel::queries;

namespace {

// Unit right triangle in the xy-plane
const Point3D kV1{0., 0., 0.};
const Point3D kV2{1., 0., 0.};
const Point3D kV3{0., 1., 0.};
const Vector3D kNormal{0., 0., 1.};  // Cross(v2-v1, v3-v1)

// Same triangle, CW winding (reversed v2/v3)
const Point3D kCwV1{0., 0., 0.};
const Point3D kCwV2{0., 1., 0.};
const Point3D kCwV3{1., 0., 0.};
const Vector3D kCwNormal{0., 0., -1.};

}  // namespace

//------------------------------------------------------------------------------
// InPlaneTriangleContainment — interior
//------------------------------------------------------------------------------

TEST(InPlaneTriangleContainment, CentroidIsInside) {
    const Point3D centroid{1. / 3., 1. / 3., 0.};
    EXPECT_TRUE(InPlaneTriangleContainment(centroid, kV1, kV2, kV3, kNormal));
}

TEST(InPlaneTriangleContainment, PointClearlyInside) {
    EXPECT_TRUE(InPlaneTriangleContainment(Point3D{0.25, 0.25, 0.}, kV1, kV2, kV3, kNormal));
}

//------------------------------------------------------------------------------
// InPlaneTriangleContainment — boundary
//------------------------------------------------------------------------------

TEST(InPlaneTriangleContainment, AtVertex) {
    EXPECT_TRUE(InPlaneTriangleContainment(kV1, kV1, kV2, kV3, kNormal));
    EXPECT_TRUE(InPlaneTriangleContainment(kV2, kV1, kV2, kV3, kNormal));
    EXPECT_TRUE(InPlaneTriangleContainment(kV3, kV1, kV2, kV3, kNormal));
}

TEST(InPlaneTriangleContainment, MidpointOfEdge) {
    EXPECT_TRUE(InPlaneTriangleContainment(Point3D{0.5, 0., 0.}, kV1, kV2, kV3, kNormal));
    EXPECT_TRUE(InPlaneTriangleContainment(Point3D{0.5, 0.5, 0.}, kV1, kV2, kV3, kNormal));
    EXPECT_TRUE(InPlaneTriangleContainment(Point3D{0., 0.5, 0.}, kV1, kV2, kV3, kNormal));
}

//------------------------------------------------------------------------------
// InPlaneTriangleContainment — outside
//------------------------------------------------------------------------------

TEST(InPlaneTriangleContainment, ClearlyOutside) {
    EXPECT_FALSE(InPlaneTriangleContainment(Point3D{5., 5., 0.}, kV1, kV2, kV3, kNormal));
}

TEST(InPlaneTriangleContainment, JustBeyondEachEdge) {
    EXPECT_FALSE(InPlaneTriangleContainment(Point3D{0.5, -0.01, 0.}, kV1, kV2, kV3, kNormal));
    EXPECT_FALSE(InPlaneTriangleContainment(Point3D{0.6, 0.6, 0.}, kV1, kV2, kV3, kNormal));
    EXPECT_FALSE(InPlaneTriangleContainment(Point3D{-0.01, 0.5, 0.}, kV1, kV2, kV3, kNormal));
}

//------------------------------------------------------------------------------
// InPlaneTriangleContainment — winding-order agnosticism
//------------------------------------------------------------------------------

TEST(InPlaneTriangleContainment, SameResultForCCWAndCW) {
    const Point3D p{0.25, 0.25, 0.};
    EXPECT_EQ(
        InPlaneTriangleContainment(p, kV1, kV2, kV3, kNormal),
        InPlaneTriangleContainment(p, kCwV1, kCwV2, kCwV3, kCwNormal)
    );
}

//------------------------------------------------------------------------------
// InPlaneTriangleContainment — height invariance (unique to 3D)
//------------------------------------------------------------------------------

TEST(InPlaneTriangleContainment, OffPlaneInteriorProjectionStillInside) {
    // Centroid lifted along the normal — still "inside" the footprint
    const Point3D lifted{1. / 3., 1. / 3., 5.};
    EXPECT_TRUE(InPlaneTriangleContainment(lifted, kV1, kV2, kV3, kNormal));
}

TEST(InPlaneTriangleContainment, OffPlaneExteriorProjectionStillOutside) {
    const Point3D lifted{0.5, -1., 5.};
    EXPECT_FALSE(InPlaneTriangleContainment(lifted, kV1, kV2, kV3, kNormal));
}

TEST(InPlaneTriangleContainment, OffPlaneAgreesWithInPlane) {
    const Point3D in_plane{0.25, 0.25, 0.};
    const Point3D lifted{0.25, 0.25, -3.};
    EXPECT_EQ(
        InPlaneTriangleContainment(in_plane, kV1, kV2, kV3, kNormal),
        InPlaneTriangleContainment(lifted, kV1, kV2, kV3, kNormal)
    );
}

//------------------------------------------------------------------------------
// InPlaneTriangleContainment — degenerate
//------------------------------------------------------------------------------

TEST(InPlaneTriangleContainment, DegenerateAllCoincident_MatchingPoint) {
    const Point3D v{0., 0., 0.};
    const Vector3D n{0., 0., 1.};
    EXPECT_TRUE(InPlaneTriangleContainment(v, v, v, v, n));
}

TEST(InPlaneTriangleContainment, DegenerateCollinear_PointOffLine) {
    const Point3D a{0., 0., 0.}, b{2., 0., 0.}, c{4., 0., 0.};
    const Vector3D n{0., 0., 1.};
    EXPECT_FALSE(InPlaneTriangleContainment(Point3D{1., 1., 0.}, a, b, c, n));
}

}  // namespace geometry_kernel::test
