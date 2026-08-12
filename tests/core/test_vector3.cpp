#include <cmath>

#include <gtest/gtest.h>

#include "core/point3.hpp"
#include "core/vector3.hpp"

namespace geometry_kernel::test {

using namespace geometry_kernel::core;

//------------------------------------------------------------------------------
// Vector3 construction
//------------------------------------------------------------------------------

TEST(Vector3, DefaultInitIsZero) {
    const Vector3D v{};
    EXPECT_DOUBLE_EQ(v.x, 0.);
    EXPECT_DOUBLE_EQ(v.y, 0.);
    EXPECT_DOUBLE_EQ(v.z, 0.);
}

TEST(Vector3, AggregateInitStoresComponents) {
    const Vector3D v{3., -4., 5.};
    EXPECT_DOUBLE_EQ(v.x, 3.);
    EXPECT_DOUBLE_EQ(v.y, -4.);
    EXPECT_DOUBLE_EQ(v.z, 5.);
}

//------------------------------------------------------------------------------
// Dot (Vector3)
//------------------------------------------------------------------------------

TEST(Dot3, OrthogonalVectors) {
    EXPECT_DOUBLE_EQ(Dot(Vector3D{1., 0., 0.}, Vector3D{0., 1., 0.}), 0.);
    EXPECT_DOUBLE_EQ(Dot(Vector3D{1., 0., 0.}, Vector3D{0., 0., 1.}), 0.);
}

TEST(Dot3, ParallelSameDirection) {
    EXPECT_DOUBLE_EQ(Dot(Vector3D{2., 0., 0.}, Vector3D{3., 0., 0.}), 6.);
}

TEST(Dot3, ParallelOppositeDirection) {
    EXPECT_DOUBLE_EQ(Dot(Vector3D{2., 0., 0.}, Vector3D{-3., 0., 0.}), -6.);
}

TEST(Dot3, KnownValue) {
    // (1,2,3) · (4,5,6) = 4 + 10 + 18 = 32
    EXPECT_DOUBLE_EQ(Dot(Vector3D{1., 2., 3.}, Vector3D{4., 5., 6.}), 32.);
}

TEST(Dot3, Commutative) {
    const Vector3D a{1.5, -2., 0.5};
    const Vector3D b{-3., 0.5, 2.};
    EXPECT_DOUBLE_EQ(Dot(a, b), Dot(b, a));
}

TEST(Dot3, SelfEqualsLengthSquared) {
    const Vector3D v{3., 4., 12.};
    EXPECT_DOUBLE_EQ(Dot(v, v), LengthSquared(v));
}

//------------------------------------------------------------------------------
// Cross (Vector3)
//------------------------------------------------------------------------------

TEST(Cross3, KnownValue) {
    // i × j = k
    const Vector3D k = Cross(Vector3D{1., 0., 0.}, Vector3D{0., 1., 0.});
    EXPECT_DOUBLE_EQ(k.x, 0.);
    EXPECT_DOUBLE_EQ(k.y, 0.);
    EXPECT_DOUBLE_EQ(k.z, 1.);
}

TEST(Cross3, SelfIsZeroVector) {
    const Vector3D v{3., 4., 5.};
    const Vector3D z = Cross(v, v);
    EXPECT_DOUBLE_EQ(z.x, 0.);
    EXPECT_DOUBLE_EQ(z.y, 0.);
    EXPECT_DOUBLE_EQ(z.z, 0.);
}

TEST(Cross3, Antisymmetric) {
    const Vector3D a{2., 3., 1.};
    const Vector3D b{-1., 4., 2.};
    const Vector3D ab = Cross(a, b);
    const Vector3D ba = Cross(b, a);
    EXPECT_DOUBLE_EQ(ab.x, -ba.x);
    EXPECT_DOUBLE_EQ(ab.y, -ba.y);
    EXPECT_DOUBLE_EQ(ab.z, -ba.z);
}

TEST(Cross3, ParallelVectorsAreZero) {
    const Vector3D z = Cross(Vector3D{2., 4., 6.}, Vector3D{3., 6., 9.});
    EXPECT_DOUBLE_EQ(z.x, 0.);
    EXPECT_DOUBLE_EQ(z.y, 0.);
    EXPECT_DOUBLE_EQ(z.z, 0.);
}

TEST(Cross3, MagnitudeIsParallelogramArea) {
    // Unit square in xy-plane spanned by (1,0,0) and (0,1,0) -> area 1
    EXPECT_DOUBLE_EQ(Length(Cross(Vector3D{1., 0., 0.}, Vector3D{0., 1., 0.})), 1.);
    // Rectangle 3×4 in xy-plane -> area 12
    EXPECT_DOUBLE_EQ(Length(Cross(Vector3D{3., 0., 0.}, Vector3D{0., 4., 0.})), 12.);
}

TEST(Cross3, ResultPerpendicularToBothOperands) {
    const Vector3D a{1., 2., 3.};
    const Vector3D b{4., -1., 2.};
    const Vector3D n = Cross(a, b);
    EXPECT_NEAR(Dot(n, a), 0., 1e-12);
    EXPECT_NEAR(Dot(n, b), 0., 1e-12);
}

//------------------------------------------------------------------------------
// LengthSquared / Length (Vector3)
//------------------------------------------------------------------------------

TEST(LengthSquared3, ThreeFourTwelveThirteen) {
    // 3-4-12 right triangle in 3D: ||(3,4,12)||^2 = 9+16+144 = 169
    EXPECT_DOUBLE_EQ(LengthSquared(Vector3D{3., 4., 12.}), 169.);
}

TEST(LengthSquared3, ZeroVector) {
    EXPECT_DOUBLE_EQ(LengthSquared(Vector3D{}), 0.);
}

TEST(Length3, ThreeFourTwelveThirteen) {
    EXPECT_DOUBLE_EQ(Length(Vector3D{3., 4., 12.}), 13.);
}

TEST(Length3, ZeroVector) {
    EXPECT_DOUBLE_EQ(Length(Vector3D{}), 0.);
}

//------------------------------------------------------------------------------
// Normalize (Vector3)
//------------------------------------------------------------------------------

TEST(Normalize3, ThreeFourTwelveThirteen) {
    const Vector3D u = Normalize(Vector3D{3., 4., 12.});
    EXPECT_DOUBLE_EQ(u.x, 3. / 13.);
    EXPECT_DOUBLE_EQ(u.y, 4. / 13.);
    EXPECT_DOUBLE_EQ(u.z, 12. / 13.);
}

TEST(Normalize3, ResultHasUnitLength) {
    const Vector3D u = Normalize(Vector3D{3., 4., 12.});
    EXPECT_DOUBLE_EQ(Length(u), 1.);
}

TEST(Normalize3, ZeroVectorReturnsZero) {
    const Vector3D u = Normalize(Vector3D{});
    EXPECT_DOUBLE_EQ(u.x, 0.);
    EXPECT_DOUBLE_EQ(u.y, 0.);
    EXPECT_DOUBLE_EQ(u.z, 0.);
}

TEST(Normalize3, AlreadyUnitIsNoOp) {
    const Vector3D v{0., 0., 1.};
    const Vector3D u = Normalize(v);
    EXPECT_DOUBLE_EQ(u.x, v.x);
    EXPECT_DOUBLE_EQ(u.y, v.y);
    EXPECT_DOUBLE_EQ(u.z, v.z);
}

//------------------------------------------------------------------------------
// operator-(Point3, Point3)
//------------------------------------------------------------------------------

TEST(Point3MinusPoint3, MatchesComponentSubtraction) {
    const Point3D a{5., 7., 9.};
    const Point3D b{2., 3., 4.};
    const Vector3D v = a - b;
    EXPECT_DOUBLE_EQ(v.x, 3.);
    EXPECT_DOUBLE_EQ(v.y, 4.);
    EXPECT_DOUBLE_EQ(v.z, 5.);
}

TEST(Point3MinusPoint3, NegatesWhenOperandsSwapped) {
    const Point3D a{5., 7., 9.};
    const Point3D b{2., 3., 4.};
    const Vector3D ab = a - b;
    const Vector3D ba = b - a;
    EXPECT_DOUBLE_EQ(ab.x, -ba.x);
    EXPECT_DOUBLE_EQ(ab.y, -ba.y);
    EXPECT_DOUBLE_EQ(ab.z, -ba.z);
}

//------------------------------------------------------------------------------
// operator+(Point3, Vector3) / operator+(Vector3, Point3)
//------------------------------------------------------------------------------

TEST(Point3PlusVector3, TranslatesPoint) {
    const Point3D p{1., 2., 3.};
    const Vector3D v{3., -4., 5.};
    const Point3D q = p + v;
    EXPECT_DOUBLE_EQ(q.x, 4.);
    EXPECT_DOUBLE_EQ(q.y, -2.);
    EXPECT_DOUBLE_EQ(q.z, 8.);
}

TEST(Vector3PlusPoint3, CommutativeWithPointPlusVector) {
    const Point3D p{1., 2., 3.};
    const Vector3D v{3., -4., 5.};
    const Point3D left = p + v;
    const Point3D right = v + p;
    EXPECT_DOUBLE_EQ(left.x, right.x);
    EXPECT_DOUBLE_EQ(left.y, right.y);
    EXPECT_DOUBLE_EQ(left.z, right.z);
}

TEST(Point3PlusVector3, ZeroVectorIsNoOp) {
    const Point3D p{1.5, -2.5, 3.5};
    const Point3D q = p + Vector3D{};
    EXPECT_DOUBLE_EQ(q.x, p.x);
    EXPECT_DOUBLE_EQ(q.y, p.y);
    EXPECT_DOUBLE_EQ(q.z, p.z);
}

TEST(Point3PlusVector3, RoundTripViaDisplacement) {
    // p + (a - p) == a
    const Point3D p{1., 2., 3.};
    const Point3D a{5., 7., 9.};
    const Point3D recovered = p + (a - p);
    EXPECT_DOUBLE_EQ(recovered.x, a.x);
    EXPECT_DOUBLE_EQ(recovered.y, a.y);
    EXPECT_DOUBLE_EQ(recovered.z, a.z);
}

//------------------------------------------------------------------------------
// operator+(Vector3, Vector3) / operator-(Vector3, Vector3)
//------------------------------------------------------------------------------

TEST(Vector3PlusVector3, ComponentWise) {
    const Vector3D sum = Vector3D{1., 2., 3.} + Vector3D{4., 5., 6.};
    EXPECT_DOUBLE_EQ(sum.x, 5.);
    EXPECT_DOUBLE_EQ(sum.y, 7.);
    EXPECT_DOUBLE_EQ(sum.z, 9.);
}

TEST(Vector3MinusVector3, ComponentWise) {
    const Vector3D diff = Vector3D{5., 7., 9.} - Vector3D{2., 3., 4.};
    EXPECT_DOUBLE_EQ(diff.x, 3.);
    EXPECT_DOUBLE_EQ(diff.y, 4.);
    EXPECT_DOUBLE_EQ(diff.z, 5.);
}

//------------------------------------------------------------------------------
// operator*(Vector3, T) / operator*(T, Vector3)
//------------------------------------------------------------------------------

TEST(Vector3TimesScalar, ScalesComponents) {
    const Vector3D v = Vector3D{3., 4., 5.} * 2.;
    EXPECT_DOUBLE_EQ(v.x, 6.);
    EXPECT_DOUBLE_EQ(v.y, 8.);
    EXPECT_DOUBLE_EQ(v.z, 10.);
}

TEST(ScalarTimesVector3, CommutativeWithVectorTimesScalar) {
    const Vector3D a = Vector3D{3., 4., 5.} * 2.5;
    const Vector3D b = 2.5 * Vector3D{3., 4., 5.};
    EXPECT_DOUBLE_EQ(a.x, b.x);
    EXPECT_DOUBLE_EQ(a.y, b.y);
    EXPECT_DOUBLE_EQ(a.z, b.z);
}

TEST(Vector3TimesScalar, ScaleByZero) {
    const Vector3D v = Vector3D{3., 4., 5.} * 0.;
    EXPECT_DOUBLE_EQ(v.x, 0.);
    EXPECT_DOUBLE_EQ(v.y, 0.);
    EXPECT_DOUBLE_EQ(v.z, 0.);
}

TEST(Vector3TimesScalar, ScaleByNegativeOne) {
    const Vector3D v = Vector3D{3., 4., 5.} * -1.;
    EXPECT_DOUBLE_EQ(v.x, -3.);
    EXPECT_DOUBLE_EQ(v.y, -4.);
    EXPECT_DOUBLE_EQ(v.z, -5.);
}

TEST(Vector3TimesScalar, ScaleByOneIsIdentity) {
    const Vector3D v{3., 4., 5.};
    const Vector3D u = v * 1.;
    EXPECT_DOUBLE_EQ(u.x, v.x);
    EXPECT_DOUBLE_EQ(u.y, v.y);
    EXPECT_DOUBLE_EQ(u.z, v.z);
}

//------------------------------------------------------------------------------
// Cross-cutting properties
//------------------------------------------------------------------------------

TEST(Properties3, CrossResultHasZeroDotWithOperands) {
    const Vector3D a{2., -1., 3.};
    const Vector3D b{0., 4., -2.};
    const Vector3D n = Cross(a, b);
    EXPECT_NEAR(Dot(n, a), 0., 1e-12);
    EXPECT_NEAR(Dot(n, b), 0., 1e-12);
}

TEST(Properties3, NormalizedLengthIsOneWhenNonZero) {
    const Vector3D v{1.5, -2.5, 3.5};
    EXPECT_NEAR(Length(Normalize(v)), 1., 1e-12);
}

TEST(Properties3, LengthSquaredIsDotWithSelf) {
    const Vector3D v{-7., 2., 1.};
    EXPECT_DOUBLE_EQ(LengthSquared(v), Dot(v, v));
}

}  // namespace geometry_kernel::test
