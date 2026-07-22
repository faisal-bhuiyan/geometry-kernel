#include <cmath>

#include <gtest/gtest.h>

#include "core/point.hpp"
#include "core/predicates.hpp"
#include "core/vector.hpp"

namespace geometry_kernel::test {

using namespace geometry_kernel::core;

//------------------------------------------------------------------------------
// Vector2 construction
//------------------------------------------------------------------------------

TEST(Vector2, DefaultInitIsZero) {
    const Vector2D v{};
    EXPECT_DOUBLE_EQ(v.x, 0.);
    EXPECT_DOUBLE_EQ(v.y, 0.);
}

TEST(Vector2, AggregateInitStoresComponents) {
    const Vector2D v{3., -4.};
    EXPECT_DOUBLE_EQ(v.x, 3.);
    EXPECT_DOUBLE_EQ(v.y, -4.);
}

//------------------------------------------------------------------------------
// Dot
//------------------------------------------------------------------------------

TEST(Dot, OrthogonalVectors) {
    EXPECT_DOUBLE_EQ(Dot(Vector2D{1., 0.}, Vector2D{0., 1.}), 0.);
}

TEST(Dot, ParallelSameDirection) {
    EXPECT_DOUBLE_EQ(Dot(Vector2D{2., 0.}, Vector2D{3., 0.}), 6.);
}

TEST(Dot, ParallelOppositeDirection) {
    EXPECT_DOUBLE_EQ(Dot(Vector2D{2., 0.}, Vector2D{-3., 0.}), -6.);
}

TEST(Dot, KnownValue) {
    // (3,4) · (4,3) = 12 + 12 = 24
    EXPECT_DOUBLE_EQ(Dot(Vector2D{3., 4.}, Vector2D{4., 3.}), 24.);
}

TEST(Dot, Commutative) {
    const Vector2D a{1.5, -2.};
    const Vector2D b{-3., 0.5};
    EXPECT_DOUBLE_EQ(Dot(a, b), Dot(b, a));
}

TEST(Dot, SelfEqualsLengthSquared) {
    const Vector2D v{3., 4.};
    EXPECT_DOUBLE_EQ(Dot(v, v), LengthSquared(v));
}

//------------------------------------------------------------------------------
// Cross
//------------------------------------------------------------------------------

TEST(Cross, KnownValue) {
    // (1,0) × (0,1) = 1
    EXPECT_DOUBLE_EQ(Cross(Vector2D{1., 0.}, Vector2D{0., 1.}), 1.);
}

TEST(Cross, SelfIsZero) {
    const Vector2D v{3., 4.};
    EXPECT_DOUBLE_EQ(Cross(v, v), 0.);
}

TEST(Cross, Antisymmetric) {
    const Vector2D a{2., 3.};
    const Vector2D b{-1., 4.};
    EXPECT_DOUBLE_EQ(Cross(a, b), -Cross(b, a));
}

TEST(Cross, ParallelVectorsAreZero) {
    EXPECT_DOUBLE_EQ(Cross(Vector2D{2., 4.}, Vector2D{3., 6.}), 0.);
}

TEST(Cross, MagnitudeIsParallelogramArea) {
    // Unit square spanned by (1,0) and (0,1) -> area 1
    EXPECT_DOUBLE_EQ(std::abs(Cross(Vector2D{1., 0.}, Vector2D{0., 1.})), 1.);
    // Rectangle 3×4 -> area 12
    EXPECT_DOUBLE_EQ(std::abs(Cross(Vector2D{3., 0.}, Vector2D{0., 4.})), 12.);
}

//------------------------------------------------------------------------------
// LengthSquared / Length
//------------------------------------------------------------------------------

TEST(LengthSquared, ThreeFourFive) {
    EXPECT_DOUBLE_EQ(LengthSquared(Vector2D{3., 4.}), 25.);
}

TEST(LengthSquared, ZeroVector) {
    EXPECT_DOUBLE_EQ(LengthSquared(Vector2D{}), 0.);
}

TEST(Length, ThreeFourFive) {
    EXPECT_DOUBLE_EQ(Length(Vector2D{3., 4.}), 5.);
}

TEST(Length, ZeroVector) {
    EXPECT_DOUBLE_EQ(Length(Vector2D{}), 0.);
}

//------------------------------------------------------------------------------
// Normalize
//------------------------------------------------------------------------------

TEST(Normalize, ThreeFourFive) {
    const Vector2D u = Normalize(Vector2D{3., 4.});
    EXPECT_DOUBLE_EQ(u.x, 0.6);
    EXPECT_DOUBLE_EQ(u.y, 0.8);
}

TEST(Normalize, ResultHasUnitLength) {
    const Vector2D u = Normalize(Vector2D{3., 4.});
    EXPECT_DOUBLE_EQ(Length(u), 1.);
}

TEST(Normalize, ZeroVectorReturnsZero) {
    // Guard against divide-by-zero: returns {0,0} rather than NaN
    const Vector2D u = Normalize(Vector2D{});
    EXPECT_DOUBLE_EQ(u.x, 0.);
    EXPECT_DOUBLE_EQ(u.y, 0.);
}

TEST(Normalize, AlreadyUnitIsNoOp) {
    const Vector2D v{0.6, 0.8};
    const Vector2D u = Normalize(v);
    EXPECT_DOUBLE_EQ(u.x, v.x);
    EXPECT_DOUBLE_EQ(u.y, v.y);
}

//------------------------------------------------------------------------------
// operator-(Point2, Point2)
//------------------------------------------------------------------------------

TEST(PointMinusPoint, MatchesComponentSubtraction) {
    const Point2D a{5., 7.};
    const Point2D b{2., 3.};
    const Vector2D v = a - b;
    EXPECT_DOUBLE_EQ(v.x, 3.);
    EXPECT_DOUBLE_EQ(v.y, 4.);
}

TEST(PointMinusPoint, NegatesWhenOperandsSwapped) {
    const Point2D a{5., 7.};
    const Point2D b{2., 3.};
    const Vector2D ab = a - b;
    const Vector2D ba = b - a;
    EXPECT_DOUBLE_EQ(ab.x, -ba.x);
    EXPECT_DOUBLE_EQ(ab.y, -ba.y);
}

//------------------------------------------------------------------------------
// operator+(Point2, Vector2) / operator+(Vector2, Point2)
//------------------------------------------------------------------------------

TEST(PointPlusVector, TranslatesPoint) {
    const Point2D p{1., 2.};
    const Vector2D v{3., -4.};
    const Point2D q = p + v;
    EXPECT_DOUBLE_EQ(q.x, 4.);
    EXPECT_DOUBLE_EQ(q.y, -2.);
}

TEST(VectorPlusPoint, CommutativeWithPointPlusVector) {
    const Point2D p{1., 2.};
    const Vector2D v{3., -4.};
    const Point2D left = p + v;
    const Point2D right = v + p;
    EXPECT_DOUBLE_EQ(left.x, right.x);
    EXPECT_DOUBLE_EQ(left.y, right.y);
}

TEST(PointPlusVector, ZeroVectorIsNoOp) {
    const Point2D p{1.5, -2.5};
    const Point2D q = p + Vector2D{};
    EXPECT_DOUBLE_EQ(q.x, p.x);
    EXPECT_DOUBLE_EQ(q.y, p.y);
}

TEST(PointPlusVector, RoundTripViaDisplacement) {
    // p + (a - p) == a
    const Point2D p{1., 2.};
    const Point2D a{5., 7.};
    const Point2D recovered = p + (a - p);
    EXPECT_DOUBLE_EQ(recovered.x, a.x);
    EXPECT_DOUBLE_EQ(recovered.y, a.y);
}

//------------------------------------------------------------------------------
// operator+(Vector2, Vector2) / operator-(Vector2, Vector2)
//------------------------------------------------------------------------------

TEST(VectorPlusVector, ComponentWise) {
    const Vector2D sum = Vector2D{1., 2.} + Vector2D{3., 4.};
    EXPECT_DOUBLE_EQ(sum.x, 4.);
    EXPECT_DOUBLE_EQ(sum.y, 6.);
}

TEST(VectorMinusVector, ComponentWise) {
    const Vector2D diff = Vector2D{5., 7.} - Vector2D{2., 3.};
    EXPECT_DOUBLE_EQ(diff.x, 3.);
    EXPECT_DOUBLE_EQ(diff.y, 4.);
}

//------------------------------------------------------------------------------
// operator*(Vector2, T) / operator*(T, Vector2)
//------------------------------------------------------------------------------

TEST(VectorTimesScalar, ScalesComponents) {
    const Vector2D v = Vector2D{3., 4.} * 2.;
    EXPECT_DOUBLE_EQ(v.x, 6.);
    EXPECT_DOUBLE_EQ(v.y, 8.);
}

TEST(ScalarTimesVector, CommutativeWithVectorTimesScalar) {
    const Vector2D a = Vector2D{3., 4.} * 2.5;
    const Vector2D b = 2.5 * Vector2D{3., 4.};
    EXPECT_DOUBLE_EQ(a.x, b.x);
    EXPECT_DOUBLE_EQ(a.y, b.y);
}

TEST(VectorTimesScalar, ScaleByZero) {
    const Vector2D v = Vector2D{3., 4.} * 0.;
    EXPECT_DOUBLE_EQ(v.x, 0.);
    EXPECT_DOUBLE_EQ(v.y, 0.);
}

TEST(VectorTimesScalar, ScaleByNegativeOne) {
    const Vector2D v = Vector2D{3., 4.} * -1.;
    EXPECT_DOUBLE_EQ(v.x, -3.);
    EXPECT_DOUBLE_EQ(v.y, -4.);
}

TEST(VectorTimesScalar, ScaleByOneIsIdentity) {
    const Vector2D v{3., 4.};
    const Vector2D u = v * 1.;
    EXPECT_DOUBLE_EQ(u.x, v.x);
    EXPECT_DOUBLE_EQ(u.y, v.y);
}

//------------------------------------------------------------------------------
// Cross-cutting properties
//------------------------------------------------------------------------------

TEST(Properties, PerpendicularConstructionHasZeroDot) {
    // (x, y) perp (-y, x)
    const Vector2D v{3., 4.};
    const Vector2D perp{-v.y, v.x};
    EXPECT_DOUBLE_EQ(Dot(v, perp), 0.);
}

TEST(Properties, SignedTriangleArea2MatchesCrossOfEdgeVectors) {
    // SignedTriangleArea2 is implemented as Cross(b - a, c - a)
    const Point2D a{0., 0.}, b{3., 0.}, c{1., 2.};
    EXPECT_DOUBLE_EQ(SignedTriangleArea2(a, b, c), Cross(b - a, c - a));
}

TEST(Properties, NormalizedLengthIsOneWhenNonZero) {
    const Vector2D v{1.5, -2.5};
    EXPECT_NEAR(Length(Normalize(v)), 1., 1e-12);
}

TEST(Properties, LengthSquaredIsDotWithSelf) {
    const Vector2D v{-7., 2.};
    EXPECT_DOUBLE_EQ(LengthSquared(v), Dot(v, v));
}

}  // namespace geometry_kernel::test
