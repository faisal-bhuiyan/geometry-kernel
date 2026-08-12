#include <cmath>

#include <gtest/gtest.h>

#include "core/point3.hpp"
#include "core/vector3.hpp"
#include "queries/closest_point3.hpp"

namespace geometry_kernel::test {

using namespace geometry_kernel::core;
using namespace geometry_kernel::queries;

namespace {

// Unit right triangle in the xy-plane
const Point3D kV1{0., 0., 0.};
const Point3D kV2{1., 0., 0.};
const Point3D kV3{0., 1., 0.};

}  // namespace

//------------------------------------------------------------------------------
// ClosestPointOnSegment (3D)
//------------------------------------------------------------------------------

TEST(ClosestPointOnSegment3, InteriorPerpendicularFoot) {
    // Segment on x-axis, point above it in z
    const Point3D a{0., 0., 0.}, b{4., 0., 0.}, p{2., 0., 3.};
    const Point3D q = ClosestPointOnSegment(p, a, b);
    EXPECT_DOUBLE_EQ(q.x, 2.);
    EXPECT_DOUBLE_EQ(q.y, 0.);
    EXPECT_DOUBLE_EQ(q.z, 0.);
}

TEST(ClosestPointOnSegment3, ClampsBeyondB) {
    const Point3D a{0., 0., 0.}, b{1., 0., 0.}, p{3., 1., 1.};
    const Point3D q = ClosestPointOnSegment(p, a, b);
    EXPECT_DOUBLE_EQ(q.x, b.x);
    EXPECT_DOUBLE_EQ(q.y, b.y);
    EXPECT_DOUBLE_EQ(q.z, b.z);
}

TEST(ClosestPointOnSegment3, ClampsBeforeA) {
    const Point3D a{0., 0., 0.}, b{1., 0., 0.}, p{-2., 1., 1.};
    const Point3D q = ClosestPointOnSegment(p, a, b);
    EXPECT_DOUBLE_EQ(q.x, a.x);
    EXPECT_DOUBLE_EQ(q.y, a.y);
    EXPECT_DOUBLE_EQ(q.z, a.z);
}

TEST(ClosestPointOnSegment3, PointEqualsA) {
    const Point3D a{1., 2., 3.}, b{5., 6., 7.};
    const Point3D q = ClosestPointOnSegment(a, a, b);
    EXPECT_DOUBLE_EQ(q.x, a.x);
    EXPECT_DOUBLE_EQ(q.y, a.y);
    EXPECT_DOUBLE_EQ(q.z, a.z);
}

TEST(ClosestPointOnSegment3, PointEqualsB) {
    const Point3D a{1., 2., 3.}, b{5., 6., 7.};
    const Point3D q = ClosestPointOnSegment(b, a, b);
    EXPECT_DOUBLE_EQ(q.x, b.x);
    EXPECT_DOUBLE_EQ(q.y, b.y);
    EXPECT_DOUBLE_EQ(q.z, b.z);
}

TEST(ClosestPointOnSegment3, DegenerateCoincidentEndpoints) {
    const Point3D a{1., 2., 3.};
    const Point3D q = ClosestPointOnSegment(Point3D{5., 5., 5.}, a, a);
    EXPECT_DOUBLE_EQ(q.x, a.x);
    EXPECT_DOUBLE_EQ(q.y, a.y);
    EXPECT_DOUBLE_EQ(q.z, a.z);
}

TEST(ClosestPointOnSegment3, OrderIndependent) {
    const Point3D a{0., 0., 0.}, b{4., 0., 0.}, p{2., 3., 1.};
    const Point3D q_ab = ClosestPointOnSegment(p, a, b);
    const Point3D q_ba = ClosestPointOnSegment(p, b, a);
    EXPECT_DOUBLE_EQ(q_ab.x, q_ba.x);
    EXPECT_DOUBLE_EQ(q_ab.y, q_ba.y);
    EXPECT_DOUBLE_EQ(q_ab.z, q_ba.z);
}

//------------------------------------------------------------------------------
// PointToSegmentDistance (3D)
//------------------------------------------------------------------------------

TEST(PointToSegmentDistance3, ThreeFourFiveInteriorFoot) {
    // Segment (0,0,0)-(5,0,0), point (0,3,4) -> foot at (0,0,0)? Wait, closest is
    // actually t=0 since projection of (0,3,4) onto x-axis is origin.
    // Better: point (2,3,4) onto (0,0,0)-(4,0,0) -> foot (2,0,0), dist = 5
    const Point3D a{0., 0., 0.}, b{4., 0., 0.}, p{2., 3., 4.};
    EXPECT_DOUBLE_EQ(PointToSegmentDistance(p, a, b), 5.);
}

TEST(PointToSegmentDistance3, PointOnSegmentIsZero) {
    const Point3D a{0., 0., 0.}, b{4., 0., 0.}, p{2., 0., 0.};
    EXPECT_DOUBLE_EQ(PointToSegmentDistance(p, a, b), 0.);
}

TEST(PointToSegmentDistance3, SqrtConsistency) {
    const Point3D a{0., 0., 0.}, b{4., 0., 0.}, p{2., 3., 4.};
    EXPECT_DOUBLE_EQ(
        PointToSegmentDistance(p, a, b), std::sqrt(PointToSegmentDistanceSquared(p, a, b))
    );
}

TEST(PointToSegmentDistance3, DefinedViaClosestPoint) {
    const Point3D a{0., 0., 0.}, b{4., 0., 0.}, p{2., 3., 1.};
    EXPECT_DOUBLE_EQ(
        PointToSegmentDistanceSquared(p, a, b),
        LengthSquared(p - ClosestPointOnSegment(p, a, b))
    );
}

TEST(PointToSegmentDistance3, BeyondEndpointEqualsDistanceToEndpoint) {
    const Point3D a{0., 0., 0.}, b{1., 0., 0.}, p{3., 4., 0.};
    EXPECT_DOUBLE_EQ(PointToSegmentDistance(p, a, b), Length(p - b));
}

TEST(PointToSegmentDistance3, NeverNegative) {
    EXPECT_GE(PointToSegmentDistance(Point3D{1., 1., 1.}, Point3D{0., 0., 0.}, Point3D{2., 0., 0.}), 0.);
}

//------------------------------------------------------------------------------
// ClosestPointOnTriangle (3D) — interior / boundary (coplanar)
//------------------------------------------------------------------------------

TEST(ClosestPointOnTriangle3, InteriorPointReturnsItself) {
    const Point3D centroid{1. / 3., 1. / 3., 0.};
    const Point3D q = ClosestPointOnTriangle(centroid, kV1, kV2, kV3);
    EXPECT_DOUBLE_EQ(q.x, centroid.x);
    EXPECT_DOUBLE_EQ(q.y, centroid.y);
    EXPECT_DOUBLE_EQ(q.z, centroid.z);
}

TEST(ClosestPointOnTriangle3, EdgeMidpointReturnsItself) {
    const Point3D mid{0.5, 0., 0.};
    const Point3D q = ClosestPointOnTriangle(mid, kV1, kV2, kV3);
    EXPECT_DOUBLE_EQ(q.x, mid.x);
    EXPECT_DOUBLE_EQ(q.y, mid.y);
    EXPECT_DOUBLE_EQ(q.z, mid.z);
}

TEST(ClosestPointOnTriangle3, VertexReturnsItself) {
    const Point3D q = ClosestPointOnTriangle(kV2, kV1, kV2, kV3);
    EXPECT_DOUBLE_EQ(q.x, kV2.x);
    EXPECT_DOUBLE_EQ(q.y, kV2.y);
    EXPECT_DOUBLE_EQ(q.z, kV2.z);
}

//------------------------------------------------------------------------------
// ClosestPointOnTriangle (3D) — exterior (coplanar)
//------------------------------------------------------------------------------

TEST(ClosestPointOnTriangle3, OutsideNearEdge) {
    const Point3D p{0.5, -1., 0.};
    const Point3D q = ClosestPointOnTriangle(p, kV1, kV2, kV3);
    const Point3D expected = ClosestPointOnSegment(p, kV1, kV2);
    EXPECT_DOUBLE_EQ(q.x, expected.x);
    EXPECT_DOUBLE_EQ(q.y, expected.y);
    EXPECT_DOUBLE_EQ(q.z, expected.z);
}

TEST(ClosestPointOnTriangle3, OutsideNearestToVertex) {
    const Point3D p{-1., -1., 0.};
    const Point3D q = ClosestPointOnTriangle(p, kV1, kV2, kV3);
    EXPECT_DOUBLE_EQ(q.x, kV1.x);
    EXPECT_DOUBLE_EQ(q.y, kV1.y);
    EXPECT_DOUBLE_EQ(q.z, kV1.z);
}

//------------------------------------------------------------------------------
// ClosestPointOnTriangle (3D) — off-plane projection (unique to 3D)
//------------------------------------------------------------------------------

TEST(ClosestPointOnTriangle3, OffPlaneInteriorProjectsOntoPlane) {
    // Centroid lifted along +z; closest point is the in-plane centroid
    const Point3D lifted{1. / 3., 1. / 3., 5.};
    const Point3D q = ClosestPointOnTriangle(lifted, kV1, kV2, kV3);
    EXPECT_NEAR(q.x, 1. / 3., 1e-12);
    EXPECT_NEAR(q.y, 1. / 3., 1e-12);
    EXPECT_NEAR(q.z, 0., 1e-12);
}

TEST(ClosestPointOnTriangle3, OffPlaneExteriorFallsBackToEdge) {
    // Point below the base edge and above the plane -> closest is on edge v1-v2
    const Point3D p{0.5, -1., 2.};
    const Point3D q = ClosestPointOnTriangle(p, kV1, kV2, kV3);
    const Point3D expected = ClosestPointOnSegment(p, kV1, kV2);
    EXPECT_NEAR(q.x, expected.x, 1e-12);
    EXPECT_NEAR(q.y, expected.y, 1e-12);
    EXPECT_NEAR(q.z, expected.z, 1e-12);
}

TEST(ClosestPointOnTriangle3, SameResultForCCWAndCW_OffPlane) {
    const Point3D lifted{1. / 3., 1. / 3., 5.};
    const Point3D q_ccw = ClosestPointOnTriangle(lifted, kV1, kV2, kV3);
    const Point3D q_cw = ClosestPointOnTriangle(lifted, kV1, kV3, kV2);
    EXPECT_NEAR(q_ccw.x, q_cw.x, 1e-12);
    EXPECT_NEAR(q_ccw.y, q_cw.y, 1e-12);
    EXPECT_NEAR(q_ccw.z, q_cw.z, 1e-12);
}

//------------------------------------------------------------------------------
// PointToTriangleDistance (3D)
//------------------------------------------------------------------------------

TEST(PointToTriangleDistance3, CoplanarInteriorIsZero) {
    const Point3D centroid{1. / 3., 1. / 3., 0.};
    EXPECT_DOUBLE_EQ(PointToTriangleDistance(centroid, kV1, kV2, kV3), 0.);
}

TEST(PointToTriangleDistance3, OffPlaneInteriorEqualsHeight) {
    // Lifted centroid: distance to plane is |z| = 5
    const Point3D lifted{1. / 3., 1. / 3., 5.};
    EXPECT_NEAR(PointToTriangleDistance(lifted, kV1, kV2, kV3), 5., 1e-12);
}

TEST(PointToTriangleDistance3, DefinedViaClosestPoint) {
    const Point3D p{2., 2., 2.};
    EXPECT_DOUBLE_EQ(
        PointToTriangleDistance(p, kV1, kV2, kV3), Length(p - ClosestPointOnTriangle(p, kV1, kV2, kV3))
    );
}

TEST(PointToTriangleDistance3, NeverNegative) {
    EXPECT_GE(PointToTriangleDistance(Point3D{0.25, 0.25, 0.}, kV1, kV2, kV3), 0.);
    EXPECT_GE(PointToTriangleDistance(Point3D{-1., -1., 3.}, kV1, kV2, kV3), 0.);
}

TEST(PointToTriangleDistance3, AtMostDistanceToAnyEdge) {
    const Point3D p{2., 2., 2.};
    const double d = PointToTriangleDistance(p, kV1, kV2, kV3);
    EXPECT_LE(d, PointToSegmentDistance(p, kV1, kV2) + 1e-12);
    EXPECT_LE(d, PointToSegmentDistance(p, kV2, kV3) + 1e-12);
    EXPECT_LE(d, PointToSegmentDistance(p, kV3, kV1) + 1e-12);
}

//------------------------------------------------------------------------------
// Cross-cutting properties
//------------------------------------------------------------------------------

TEST(ClosestPointProperties3, TriangleClosestPointIsIdempotent) {
    const Point3D p{2., -1., 3.};
    const Point3D q = ClosestPointOnTriangle(p, kV1, kV2, kV3);
    const Point3D q2 = ClosestPointOnTriangle(q, kV1, kV2, kV3);
    EXPECT_NEAR(q2.x, q.x, 1e-12);
    EXPECT_NEAR(q2.y, q.y, 1e-12);
    EXPECT_NEAR(q2.z, q.z, 1e-12);
}

TEST(ClosestPointProperties3, TriangleRotationInvariance) {
    const Point3D p{2., -1., 1.};
    const Point3D q123 = ClosestPointOnTriangle(p, kV1, kV2, kV3);
    const Point3D q231 = ClosestPointOnTriangle(p, kV2, kV3, kV1);
    const Point3D q312 = ClosestPointOnTriangle(p, kV3, kV1, kV2);
    EXPECT_NEAR(q123.x, q231.x, 1e-12);
    EXPECT_NEAR(q123.y, q231.y, 1e-12);
    EXPECT_NEAR(q123.z, q231.z, 1e-12);
    EXPECT_NEAR(q123.x, q312.x, 1e-12);
    EXPECT_NEAR(q123.y, q312.y, 1e-12);
    EXPECT_NEAR(q123.z, q312.z, 1e-12);
}

}  // namespace geometry_kernel::test
