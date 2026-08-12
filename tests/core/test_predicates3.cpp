#include <gtest/gtest.h>

#include "core/area.hpp"
#include "core/point2.hpp"
#include "core/point3.hpp"
#include "core/predicates3.hpp"
#include "core/tolerance.hpp"
#include "core/vector3.hpp"

namespace geometry_kernel::test {

using namespace geometry_kernel::core;

namespace {

// Unit right triangle in the xy-plane: (0,0,0)-(1,0,0)-(0,1,0), normal = +z
const Point3D kV1{0., 0., 0.};
const Point3D kV2{1., 0., 0.};
const Point3D kV3{0., 1., 0.};
const Vector3D kNormal{0., 0., 1.};

}  // namespace

//------------------------------------------------------------------------------
// EdgeSign3 — documented cases
//------------------------------------------------------------------------------

TEST(EdgeSign3, PointLeftOfDirectedEdge_Positive) {
    // Looking against +z, edge v1 -> v2 along +x; point above the edge (positive y) is left
    EXPECT_GT(EdgeSign3(kV1, kV2, Point3D{0.5, 0.5, 0.}, kNormal), 0.);
}

TEST(EdgeSign3, PointRightOfDirectedEdge_Negative) {
    // Point below the edge (negative y) is right
    EXPECT_LT(EdgeSign3(kV1, kV2, Point3D{0.5, -0.5, 0.}, kNormal), 0.);
}

TEST(EdgeSign3, PointCollinearWithEdge_Zero) {
    EXPECT_EQ(RobustSign(EdgeSign3(kV1, kV2, Point3D{0.5, 0., 0.}, kNormal)), 0);
}

//------------------------------------------------------------------------------
// EdgeSign3 — properties
//------------------------------------------------------------------------------

TEST(EdgeSign3, AntisymmetricUnderEdgeReversal) {
    const Point3D p{0.5, 0.5, 0.};
    EXPECT_DOUBLE_EQ(EdgeSign3(kV1, kV2, p, kNormal), -EdgeSign3(kV2, kV1, p, kNormal));
}

TEST(EdgeSign3, HeightAlongNormalDoesNotChangeSign) {
    // Same in-plane position, lifted along the normal
    const Point3D in_plane{0.5, 0.5, 0.};
    const Point3D lifted{0.5, 0.5, 10.};
    EXPECT_EQ(
        RobustSign(EdgeSign3(kV1, kV2, in_plane, kNormal)),
        RobustSign(EdgeSign3(kV1, kV2, lifted, kNormal))
    );
    EXPECT_DOUBLE_EQ(
        EdgeSign3(kV1, kV2, in_plane, kNormal), EdgeSign3(kV1, kV2, lifted, kNormal)
    );
}

TEST(EdgeSign3, AgreesWithSignedTriangleAreaTimes2WhenCoplanar) {
    // For z=0 and normal=(0,0,1), EdgeSign3(vi,vj,p,n) == SignedTriangleAreaTimes2(vi_2d, vj_2d, p_2d)
    // because Cross(edge, p-vi) = (0,0, Cross2D) and Dot with (0,0,1) extracts that z-component.
    const Point2D a{0., 0.}, b{1., 0.}, c{0.5, 0.5};
    const Point3D a3{0., 0., 0.}, b3{1., 0., 0.}, c3{0.5, 0.5, 0.};
    EXPECT_DOUBLE_EQ(
        EdgeSign3(a3, b3, c3, kNormal), SignedTriangleAreaTimes2(a, b, c)
    );
}

TEST(EdgeSign3, NegatedNormalFlipsSign) {
    const Point3D p{0.5, 0.5, 0.};
    const Vector3D neg_normal{0., 0., -1.};
    EXPECT_DOUBLE_EQ(EdgeSign3(kV1, kV2, p, kNormal), -EdgeSign3(kV1, kV2, p, neg_normal));
}

TEST(EdgeSign3, InteriorPointOfTriangleHasConsistentSigns) {
    // For CCW triangle with +z normal, all three EdgeSign3 values are positive at the centroid
    const Point3D centroid{1. / 3., 1. / 3., 0.};
    EXPECT_GT(EdgeSign3(kV1, kV2, centroid, kNormal), 0.);
    EXPECT_GT(EdgeSign3(kV2, kV3, centroid, kNormal), 0.);
    EXPECT_GT(EdgeSign3(kV3, kV1, centroid, kNormal), 0.);
}

TEST(EdgeSign3, OutsidePointHasMixedSigns) {
    const Point3D outside{0.5, -0.5, 0.};  // below base edge
    const int s1 = RobustSign(EdgeSign3(kV1, kV2, outside, kNormal));
    const int s2 = RobustSign(EdgeSign3(kV2, kV3, outside, kNormal));
    const int s3 = RobustSign(EdgeSign3(kV3, kV1, outside, kNormal));
    const bool has_neg = (s1 < 0) || (s2 < 0) || (s3 < 0);
    const bool has_pos = (s1 > 0) || (s2 > 0) || (s3 > 0);
    EXPECT_TRUE(has_neg && has_pos);
}

}  // namespace geometry_kernel::test
