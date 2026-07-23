#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include "core/point.hpp"
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

}  // namespace geometry_kernel::test
