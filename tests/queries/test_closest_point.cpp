#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include "core/point.hpp"
#include "core/test_fixtures.hpp"
#include "core/tolerance.hpp"
#include "core/vector.hpp"
#include "queries/closest_point.hpp"

namespace geometry_kernel::test {

using namespace geometry_kernel::core;
using namespace geometry_kernel::queries;

//------------------------------------------------------------------------------
// ClosestPointOnLine
//------------------------------------------------------------------------------

TEST(ClosestPointOnLine, HorizontalLine) {
    // Line (0,0)-(4,0), point (2,3) -> foot (2,0)
    const Point2D q = ClosestPointOnLine(Point2D{2., 3.}, Point2D{0., 0.}, Point2D{4., 0.});
    EXPECT_DOUBLE_EQ(q.x, 2.);
    EXPECT_DOUBLE_EQ(q.y, 0.);
}

TEST(ClosestPointOnLine, VerticalLine) {
    // Line (0,0)-(0,4), point (3,2) -> foot (0,2)
    const Point2D q = ClosestPointOnLine(Point2D{3., 2.}, Point2D{0., 0.}, Point2D{0., 4.});
    EXPECT_DOUBLE_EQ(q.x, 0.);
    EXPECT_DOUBLE_EQ(q.y, 2.);
}

TEST(ClosestPointOnLine, DiagonalFortyFiveDegrees) {
    // Line (0,0)-(1,1), point (1,0) -> foot (0.5, 0.5)
    const Point2D q = ClosestPointOnLine(Point2D{1., 0.}, Point2D{0., 0.}, Point2D{1., 1.});
    EXPECT_DOUBLE_EQ(q.x, 0.5);
    EXPECT_DOUBLE_EQ(q.y, 0.5);
}

TEST(ClosestPointOnLine, PointAlreadyOnLine) {
    const Point2D a{0., 0.}, b{4., 0.}, p{2., 0.};
    const Point2D q = ClosestPointOnLine(p, a, b);
    EXPECT_DOUBLE_EQ(q.x, p.x);
    EXPECT_DOUBLE_EQ(q.y, p.y);
}

TEST(ClosestPointOnLine, PointEqualsA) {
    const Point2D a{1., 2.}, b{5., 6.};
    const Point2D q = ClosestPointOnLine(a, a, b);
    EXPECT_DOUBLE_EQ(q.x, a.x);
    EXPECT_DOUBLE_EQ(q.y, a.y);
}

TEST(ClosestPointOnLine, PointEqualsB) {
    const Point2D a{1., 2.}, b{5., 6.};
    const Point2D q = ClosestPointOnLine(b, a, b);
    EXPECT_DOUBLE_EQ(q.x, b.x);
    EXPECT_DOUBLE_EQ(q.y, b.y);
}

TEST(ClosestPointOnLine, FootBeforeA_Unclamped) {
    // Point left of a on horizontal line; foot has t < 0 and stays off-segment
    const Point2D a{0., 0.}, b{4., 0.}, p{-2., 3.};
    const Point2D q = ClosestPointOnLine(p, a, b);
    EXPECT_DOUBLE_EQ(q.x, -2.);
    EXPECT_DOUBLE_EQ(q.y, 0.);
}

TEST(ClosestPointOnLine, FootBeyondB_Unclamped) {
    // Point right of b on horizontal line; foot has t > 1 and stays off-segment
    const Point2D a{0., 0.}, b{4., 0.}, p{6., 3.};
    const Point2D q = ClosestPointOnLine(p, a, b);
    EXPECT_DOUBLE_EQ(q.x, 6.);
    EXPECT_DOUBLE_EQ(q.y, 0.);
}

TEST(ClosestPointOnLine, DegenerateCoincidentEndpoints) {
    const Point2D a{1., 2.};
    const Point2D q = ClosestPointOnLine(Point2D{5., 5.}, a, a);
    EXPECT_DOUBLE_EQ(q.x, a.x);
    EXPECT_DOUBLE_EQ(q.y, a.y);
}

TEST(ClosestPointOnLine, DegenerateNearToleranceBoundary) {
    // LengthSquared(direction) == kTolerance is treated as zero by RobustSign
    const Point2D a{0., 0.};
    const Point2D b{std::sqrt(kTolerance<double>), 0.};
    EXPECT_EQ(RobustSign(LengthSquared(b - a)), 0);
    const Point2D q = ClosestPointOnLine(Point2D{1., 1.}, a, b);
    EXPECT_DOUBLE_EQ(q.x, a.x);
    EXPECT_DOUBLE_EQ(q.y, a.y);
}

TEST(ClosestPointOnLine, JustAboveToleranceIsNotDegenerate) {
    // LengthSquared just above kTolerance -> normal projection, not early return
    const double eps = std::numeric_limits<double>::epsilon();
    const Point2D a{0., 0.};
    const Point2D b{std::sqrt(kTolerance<double> + eps), 0.};
    EXPECT_EQ(RobustSign(LengthSquared(b - a)), 1);
    const Point2D q = ClosestPointOnLine(Point2D{0., 1.}, a, b);
    EXPECT_NEAR(q.x, 0., 1e-9);
    EXPECT_NEAR(q.y, 0., 1e-9);
}

TEST(ClosestPointOnLine, OrderIndependent) {
    const Point2D a{0., 0.}, b{4., 2.}, p{1., 3.};
    const Point2D q_ab = ClosestPointOnLine(p, a, b);
    const Point2D q_ba = ClosestPointOnLine(p, b, a);
    EXPECT_DOUBLE_EQ(q_ab.x, q_ba.x);
    EXPECT_DOUBLE_EQ(q_ab.y, q_ba.y);
}

TEST(ClosestPointOnLine, NegativeCoordinates) {
    const Point2D q = ClosestPointOnLine(Point2D{-2., -3.}, Point2D{-4., 0.}, Point2D{0., 0.});
    EXPECT_DOUBLE_EQ(q.x, -2.);
    EXPECT_DOUBLE_EQ(q.y, 0.);
}

TEST(ClosestPointOnLine, LargeCoordinates) {
    const Point2D a{1e6, 0.}, b{1e6 + 4., 0.}, p{1e6 + 2., 3.};
    const Point2D q = ClosestPointOnLine(p, a, b);
    EXPECT_NEAR(q.x, 1e6 + 2., 1e-6);
    EXPECT_NEAR(q.y, 0., 1e-6);
}

//------------------------------------------------------------------------------
// ClosestPointOnSegment
//------------------------------------------------------------------------------

TEST(ClosestPointOnSegment, AgreesWithLineWhenTInRange) {
    const Point2D a{0., 0.}, b{4., 0.}, p{2., 3.};
    const Point2D on_line = ClosestPointOnLine(p, a, b);
    const Point2D on_seg = ClosestPointOnSegment(p, a, b);
    EXPECT_DOUBLE_EQ(on_seg.x, on_line.x);
    EXPECT_DOUBLE_EQ(on_seg.y, on_line.y);
}

TEST(ClosestPointOnSegment, ClampsBeyondB) {
    const Point2D a{0., 0.}, b{4., 0.}, p{6., 3.};
    const Point2D q = ClosestPointOnSegment(p, a, b);
    EXPECT_DOUBLE_EQ(q.x, b.x);
    EXPECT_DOUBLE_EQ(q.y, b.y);
}

TEST(ClosestPointOnSegment, ClampsBeforeA) {
    const Point2D a{0., 0.}, b{4., 0.}, p{-2., 3.};
    const Point2D q = ClosestPointOnSegment(p, a, b);
    EXPECT_DOUBLE_EQ(q.x, a.x);
    EXPECT_DOUBLE_EQ(q.y, a.y);
}

TEST(ClosestPointOnSegment, PointEqualsA) {
    const Point2D a{1., 2.}, b{5., 6.};
    const Point2D q = ClosestPointOnSegment(a, a, b);
    EXPECT_DOUBLE_EQ(q.x, a.x);
    EXPECT_DOUBLE_EQ(q.y, a.y);
}

TEST(ClosestPointOnSegment, PointEqualsB) {
    const Point2D a{1., 2.}, b{5., 6.};
    const Point2D q = ClosestPointOnSegment(b, a, b);
    EXPECT_DOUBLE_EQ(q.x, b.x);
    EXPECT_DOUBLE_EQ(q.y, b.y);
}

TEST(ClosestPointOnSegment, InteriorPerpendicularFoot) {
    const Point2D a{0., 0.}, b{4., 0.}, p{2., 5.};
    const Point2D q = ClosestPointOnSegment(p, a, b);
    EXPECT_DOUBLE_EQ(q.x, 2.);
    EXPECT_DOUBLE_EQ(q.y, 0.);
}

TEST(ClosestPointOnSegment, DegenerateCoincidentEndpoints) {
    const Point2D a{1., 2.};
    const Point2D q = ClosestPointOnSegment(Point2D{5., 5.}, a, a);
    EXPECT_DOUBLE_EQ(q.x, a.x);
    EXPECT_DOUBLE_EQ(q.y, a.y);
}

TEST(ClosestPointOnSegment, OrderIndependent) {
    // Same segment geometrically regardless of endpoint order
    const Point2D a{0., 0.}, b{4., 0.}, p{6., 3.};
    const Point2D q_ab = ClosestPointOnSegment(p, a, b);
    const Point2D q_ba = ClosestPointOnSegment(p, b, a);
    EXPECT_DOUBLE_EQ(q_ab.x, q_ba.x);
    EXPECT_DOUBLE_EQ(q_ab.y, q_ba.y);
}

//------------------------------------------------------------------------------
// PointToSegmentDistanceSquared / PointToSegmentDistance
//------------------------------------------------------------------------------

TEST(PointToSegmentDistance, ThreeFourFiveInteriorFoot) {
    // Segment on x-axis, point at (2,3) -> perpendicular distance 3
    const Point2D a{0., 0.}, b{4., 0.}, p{2., 3.};
    EXPECT_DOUBLE_EQ(PointToSegmentDistance(p, a, b), 3.);
    EXPECT_DOUBLE_EQ(PointToSegmentDistanceSquared(p, a, b), 9.);
}

TEST(PointToSegmentDistance, BeyondEndpointEqualsDistanceToEndpoint) {
    const Point2D a{0., 0.}, b{4., 0.}, p{6., 3.};
    EXPECT_DOUBLE_EQ(PointToSegmentDistance(p, a, b), Length(p - b));
    EXPECT_DOUBLE_EQ(PointToSegmentDistanceSquared(p, a, b), LengthSquared(p - b));
}

TEST(PointToSegmentDistance, PointOnSegmentIsZero) {
    const Point2D a{0., 0.}, b{4., 0.}, p{2., 0.};
    EXPECT_DOUBLE_EQ(PointToSegmentDistance(p, a, b), 0.);
    EXPECT_DOUBLE_EQ(PointToSegmentDistanceSquared(p, a, b), 0.);
}

TEST(PointToSegmentDistance, SqrtConsistency) {
    const Point2D a{0., 0.}, b{4., 2.}, p{1., 5.};
    EXPECT_DOUBLE_EQ(
        PointToSegmentDistance(p, a, b), std::sqrt(PointToSegmentDistanceSquared(p, a, b))
    );
}

TEST(PointToSegmentDistanceSquared, DefinedViaClosestPoint) {
    const Point2D a{0., 0.}, b{4., 0.}, p{6., 3.};
    EXPECT_DOUBLE_EQ(
        PointToSegmentDistanceSquared(p, a, b), LengthSquared(p - ClosestPointOnSegment(p, a, b))
    );
}

TEST(PointToSegmentDistance, DegenerateEqualsDistanceToA) {
    const Point2D a{1., 2.}, p{4., 6.};
    EXPECT_DOUBLE_EQ(PointToSegmentDistance(p, a, a), Length(p - a));
}

TEST(PointToSegmentDistance, NeverNegative) {
    const Point2D a{0., 0.}, b{4., 0.};
    EXPECT_GE(PointToSegmentDistance(Point2D{2., 3.}, a, b), 0.);
    EXPECT_GE(PointToSegmentDistance(Point2D{-1., 0.}, a, b), 0.);
    EXPECT_GE(PointToSegmentDistanceSquared(Point2D{6., 3.}, a, b), 0.);
}

TEST(PointToSegmentDistance, AtMostDistanceToEitherEndpoint) {
    const Point2D a{0., 0.}, b{4., 0.}, p{2., 3.};
    const double d = PointToSegmentDistance(p, a, b);
    EXPECT_LE(d, Length(p - a));
    EXPECT_LE(d, Length(p - b));
}

//------------------------------------------------------------------------------
// Cross-cutting properties
//------------------------------------------------------------------------------

TEST(ClosestPointProperties, LineProjectionIsCollinearWithAB) {
    const Point2D a{1., 2.}, b{5., -1.}, p{-3., 7.};
    const Point2D q = ClosestPointOnLine(p, a, b);
    EXPECT_NEAR(Cross(q - a, b - a), 0., 1e-12);
}

TEST(ClosestPointProperties, SegmentResultParameterInUnitInterval) {
    const Point2D a{0., 0.}, b{4., 2.};
    for (const Point2D& p : {Point2D{-2., 3.}, Point2D{2., 5.}, Point2D{10., -1.}}) {
        const Point2D q = ClosestPointOnSegment(p, a, b);
        const double t = Dot(q - a, b - a) / LengthSquared(b - a);
        EXPECT_GE(t, -1e-12);
        EXPECT_LE(t, 1. + 1e-12);
    }
}

TEST(ClosestPointProperties, DistanceOrderIndependent) {
    const Point2D a{0., 0.}, b{4., 2.}, p{1., 5.};
    EXPECT_DOUBLE_EQ(PointToSegmentDistance(p, a, b), PointToSegmentDistance(p, b, a));
}

TEST(ClosestPointProperties, SegmentDistanceAtLeastLineDistance) {
    // Clamping can only increase distance relative to the unconstrained line foot
    const Point2D a{0., 0.}, b{4., 0.};
    for (const Point2D& p : {Point2D{2., 3.}, Point2D{-2., 3.}, Point2D{6., 3.}}) {
        const double seg_sq = PointToSegmentDistanceSquared(p, a, b);
        const double line_sq = LengthSquared(p - ClosestPointOnLine(p, a, b));
        EXPECT_GE(seg_sq, line_sq - 1e-12);
    }
}

//------------------------------------------------------------------------------
// ClosestPointOnTriangle — interior / boundary
//------------------------------------------------------------------------------

TEST(ClosestPointOnTriangle, InteriorPointReturnsItself) {
    const Point2D& v1 = kCcwTriangle[0];
    const Point2D& v2 = kCcwTriangle[1];
    const Point2D& v3 = kCcwTriangle[2];
    const Point2D centroid{1. / 3., 1. / 3.};
    const Point2D q = ClosestPointOnTriangle(centroid, v1, v2, v3);
    EXPECT_DOUBLE_EQ(q.x, centroid.x);
    EXPECT_DOUBLE_EQ(q.y, centroid.y);
}

TEST(ClosestPointOnTriangle, EdgeMidpointReturnsItself) {
    const Point2D& v1 = kCcwTriangle[0];
    const Point2D& v2 = kCcwTriangle[1];
    const Point2D& v3 = kCcwTriangle[2];
    const Point2D mid{0.5, 0.};  // midpoint of edge v1-v2
    const Point2D q = ClosestPointOnTriangle(mid, v1, v2, v3);
    EXPECT_DOUBLE_EQ(q.x, mid.x);
    EXPECT_DOUBLE_EQ(q.y, mid.y);
}

TEST(ClosestPointOnTriangle, VertexReturnsItself) {
    const Point2D& v1 = kCcwTriangle[0];
    const Point2D& v2 = kCcwTriangle[1];
    const Point2D& v3 = kCcwTriangle[2];
    const Point2D q = ClosestPointOnTriangle(v2, v1, v2, v3);
    EXPECT_DOUBLE_EQ(q.x, v2.x);
    EXPECT_DOUBLE_EQ(q.y, v2.y);
}

//------------------------------------------------------------------------------
// ClosestPointOnTriangle — exterior edge-interior regions (one per branch)
//------------------------------------------------------------------------------

TEST(ClosestPointOnTriangle, OutsideNearEdgeV1V2) {
    // Unit right triangle; point below base edge v1-v2
    const Point2D& v1 = kCcwTriangle[0];
    const Point2D& v2 = kCcwTriangle[1];
    const Point2D& v3 = kCcwTriangle[2];
    const Point2D p{0.5, -1.};
    const Point2D q = ClosestPointOnTriangle(p, v1, v2, v3);
    const Point2D expected = ClosestPointOnSegment(p, v1, v2);
    EXPECT_DOUBLE_EQ(q.x, expected.x);
    EXPECT_DOUBLE_EQ(q.y, expected.y);
}

TEST(ClosestPointOnTriangle, OutsideNearEdgeV2V3) {
    // Point outside near the hypotenuse v2(1,0)-v3(0,1); midpoint is (0.5,0.5)
    // Direction along edge is (-1,1); outward normal (toward exterior) is (1,1)/sqrt(2)
    const Point2D& v1 = kCcwTriangle[0];
    const Point2D& v2 = kCcwTriangle[1];
    const Point2D& v3 = kCcwTriangle[2];
    const Point2D p{1., 1.};  // outside, past the hypotenuse
    const Point2D q = ClosestPointOnTriangle(p, v1, v2, v3);
    const Point2D expected = ClosestPointOnSegment(p, v2, v3);
    EXPECT_DOUBLE_EQ(q.x, expected.x);
    EXPECT_DOUBLE_EQ(q.y, expected.y);
}

TEST(ClosestPointOnTriangle, OutsideNearEdgeV3V1) {
    // Point left of vertical-ish edge v3(0,1)-v1(0,0) which is the y-axis segment
    const Point2D& v1 = kCcwTriangle[0];
    const Point2D& v2 = kCcwTriangle[1];
    const Point2D& v3 = kCcwTriangle[2];
    const Point2D p{-1., 0.5};
    const Point2D q = ClosestPointOnTriangle(p, v1, v2, v3);
    const Point2D expected = ClosestPointOnSegment(p, v3, v1);
    EXPECT_DOUBLE_EQ(q.x, expected.x);
    EXPECT_DOUBLE_EQ(q.y, expected.y);
}

//------------------------------------------------------------------------------
// ClosestPointOnTriangle — exterior vertex regions
//------------------------------------------------------------------------------

TEST(ClosestPointOnTriangle, OutsideNearestToV1) {
    const Point2D& v1 = kCcwTriangle[0];
    const Point2D& v2 = kCcwTriangle[1];
    const Point2D& v3 = kCcwTriangle[2];
    const Point2D p{-1., -1.};
    const Point2D q = ClosestPointOnTriangle(p, v1, v2, v3);
    EXPECT_DOUBLE_EQ(q.x, v1.x);
    EXPECT_DOUBLE_EQ(q.y, v1.y);
}

TEST(ClosestPointOnTriangle, OutsideNearestToV2) {
    const Point2D& v1 = kCcwTriangle[0];
    const Point2D& v2 = kCcwTriangle[1];
    const Point2D& v3 = kCcwTriangle[2];
    const Point2D p{2., -1.};
    const Point2D q = ClosestPointOnTriangle(p, v1, v2, v3);
    EXPECT_DOUBLE_EQ(q.x, v2.x);
    EXPECT_DOUBLE_EQ(q.y, v2.y);
}

TEST(ClosestPointOnTriangle, OutsideNearestToV3) {
    const Point2D& v1 = kCcwTriangle[0];
    const Point2D& v2 = kCcwTriangle[1];
    const Point2D& v3 = kCcwTriangle[2];
    const Point2D p{-1., 2.};
    const Point2D q = ClosestPointOnTriangle(p, v1, v2, v3);
    EXPECT_DOUBLE_EQ(q.x, v3.x);
    EXPECT_DOUBLE_EQ(q.y, v3.y);
}

//------------------------------------------------------------------------------
// ClosestPointOnTriangle — degenerate triangles
//------------------------------------------------------------------------------

TEST(ClosestPointOnTriangle, DegenerateAllCoincident_ReturnsQueryPoint) {
    // Known limitation inherited from PointInTriangle: when all vertices coincide,
    // PointInTriangle always returns true, so ClosestPointOnTriangle early-returns
    // the query point itself (implying zero distance to a single coincident point).
    const Point2D v{0., 0.};
    const Point2D p{1., 0.};
    const Point2D q = ClosestPointOnTriangle(p, v, v, v);
    EXPECT_DOUBLE_EQ(q.x, p.x);
    EXPECT_DOUBLE_EQ(q.y, p.y);
}

TEST(ClosestPointOnTriangle, DegenerateCollinearVertices) {
    // Three distinct collinear points; no crash/NaN, result lies on the line
    const Point2D v1{0., 0.}, v2{2., 0.}, v3{4., 0.};
    const Point2D p{2., 3.};
    const Point2D q = ClosestPointOnTriangle(p, v1, v2, v3);
    EXPECT_NEAR(q.y, 0., 1e-12);  // on the x-axis
    EXPECT_GE(q.x, 0. - 1e-12);
    EXPECT_LE(q.x, 4. + 1e-12);
    // Closest should be the foot onto the middle of the span
    EXPECT_NEAR(q.x, 2., 1e-12);
}

//------------------------------------------------------------------------------
// ClosestPointOnTriangle — winding-order agnosticism
//------------------------------------------------------------------------------

TEST(ClosestPointOnTriangle, SameResultForCCWAndCW_Interior) {
    const Point2D centroid{1. / 3., 1. / 3.};
    const Point2D q_ccw =
        ClosestPointOnTriangle(centroid, kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]);
    const Point2D q_cw =
        ClosestPointOnTriangle(centroid, kCwTriangle[0], kCwTriangle[1], kCwTriangle[2]);
    EXPECT_DOUBLE_EQ(q_ccw.x, q_cw.x);
    EXPECT_DOUBLE_EQ(q_ccw.y, q_cw.y);
}

TEST(ClosestPointOnTriangle, SameResultForCCWAndCW_Exterior) {
    const Point2D p{0.5, -1.};
    const Point2D q_ccw =
        ClosestPointOnTriangle(p, kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]);
    const Point2D q_cw = ClosestPointOnTriangle(p, kCwTriangle[0], kCwTriangle[1], kCwTriangle[2]);
    EXPECT_DOUBLE_EQ(q_ccw.x, q_cw.x);
    EXPECT_DOUBLE_EQ(q_ccw.y, q_cw.y);
}

//------------------------------------------------------------------------------
// PointToTriangleDistance
//------------------------------------------------------------------------------

TEST(PointToTriangleDistance, InteriorIsZero) {
    const Point2D centroid{1. / 3., 1. / 3.};
    EXPECT_DOUBLE_EQ(
        PointToTriangleDistance(centroid, kCcwTriangle[0], kCcwTriangle[1], kCcwTriangle[2]), 0.
    );
}

TEST(PointToTriangleDistance, DefinedViaClosestPoint) {
    const Point2D p{2., 2.};
    const Point2D& v1 = kCcwTriangle[0];
    const Point2D& v2 = kCcwTriangle[1];
    const Point2D& v3 = kCcwTriangle[2];
    EXPECT_DOUBLE_EQ(
        PointToTriangleDistance(p, v1, v2, v3), Length(p - ClosestPointOnTriangle(p, v1, v2, v3))
    );
}

TEST(PointToTriangleDistance, KnownNumericCase) {
    // Right triangle (0,0)-(4,0)-(0,3); point (2,-3) is 3 units below the base
    const Point2D v1{0., 0.}, v2{4., 0.}, v3{0., 3.};
    const Point2D p{2., -3.};
    EXPECT_DOUBLE_EQ(PointToTriangleDistance(p, v1, v2, v3), 3.);
}

TEST(PointToTriangleDistance, NeverNegative) {
    const Point2D& v1 = kCcwTriangle[0];
    const Point2D& v2 = kCcwTriangle[1];
    const Point2D& v3 = kCcwTriangle[2];
    EXPECT_GE(PointToTriangleDistance(Point2D{0.25, 0.25}, v1, v2, v3), 0.);
    EXPECT_GE(PointToTriangleDistance(Point2D{-1., -1.}, v1, v2, v3), 0.);
    EXPECT_GE(PointToTriangleDistance(Point2D{5., 5.}, v1, v2, v3), 0.);
}

TEST(PointToTriangleDistance, AtMostDistanceToAnyEdge) {
    const Point2D& v1 = kCcwTriangle[0];
    const Point2D& v2 = kCcwTriangle[1];
    const Point2D& v3 = kCcwTriangle[2];
    const Point2D p{2., 2.};
    const double d = PointToTriangleDistance(p, v1, v2, v3);
    EXPECT_LE(d, PointToSegmentDistance(p, v1, v2) + 1e-12);
    EXPECT_LE(d, PointToSegmentDistance(p, v2, v3) + 1e-12);
    EXPECT_LE(d, PointToSegmentDistance(p, v3, v1) + 1e-12);
}

//------------------------------------------------------------------------------
// ClosestPointOnTriangle — cross-cutting properties
//------------------------------------------------------------------------------

TEST(ClosestPointProperties, TriangleRotationInvariance) {
    // Same triangle, cyclic vertex permutations must agree
    const Point2D& v1 = kCcwTriangle[0];
    const Point2D& v2 = kCcwTriangle[1];
    const Point2D& v3 = kCcwTriangle[2];
    const Point2D p{2., -1.};
    const Point2D q123 = ClosestPointOnTriangle(p, v1, v2, v3);
    const Point2D q231 = ClosestPointOnTriangle(p, v2, v3, v1);
    const Point2D q312 = ClosestPointOnTriangle(p, v3, v1, v2);
    EXPECT_DOUBLE_EQ(q123.x, q231.x);
    EXPECT_DOUBLE_EQ(q123.y, q231.y);
    EXPECT_DOUBLE_EQ(q123.x, q312.x);
    EXPECT_DOUBLE_EQ(q123.y, q312.y);
}

TEST(ClosestPointProperties, TriangleClosestPointIsIdempotent) {
    const Point2D& v1 = kCcwTriangle[0];
    const Point2D& v2 = kCcwTriangle[1];
    const Point2D& v3 = kCcwTriangle[2];
    const Point2D p{2., -1.};
    const Point2D q = ClosestPointOnTriangle(p, v1, v2, v3);
    const Point2D q2 = ClosestPointOnTriangle(q, v1, v2, v3);
    EXPECT_DOUBLE_EQ(q2.x, q.x);
    EXPECT_DOUBLE_EQ(q2.y, q.y);
}

}  // namespace geometry_kernel::test
