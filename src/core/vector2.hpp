#pragma once

#include <cmath>

#include "point2.hpp"
#include "tolerance.hpp"
#include "types.hpp"

namespace geometry_kernel::core {

//---------------------------------------------------------------------------
// Vector in the plane (2D)
//---------------------------------------------------------------------------

/// @brief 2D vector in the plane
template <ScalarType T>
struct Vector2 {
    T x{};
    T y{};
};

/// @brief 2D vector in the plane with double precision
using Vector2D = Vector2<double>;

/**
 * @brief Computes the dot product (scalar product) of two vectors.
 *
 * --------------------------------------------------------------------------------------
 *
 * Projection view -> drop a perpendicular from the tip of v1 onto the line of v2.
 *
 *               * v1
 *              /:
 *             / :
 *            /  :
 *           /θ  :
 *          *----+---------------> v2
 *        origin
 *          |<-->|  |v1| cos θ  ->  signed length of v1 measured along v2
 *
 *     Dot(v1, v2) = |v2| × (|v1| cos θ) = |v1| |v2| cos θ
 *
 * The sign alone answers "do these point the same general way" -> which is all
 * most callers need -> acute (positive), perpendicular (zero), obtuse (negative).
 *
 *        * v1              * v1           v1 *
 *       /                  |                  \
 *      *-----> v2          *-----> v2          *-----> v2
 *
 *      θ < 90°             θ = 90°             θ > 90°
 *      Dot > 0             Dot = 0             Dot < 0
 *      acute               perpendicular       obtuse
 *
 * --------------------------------------------------------------------------------------
 *
 * Note that callers rarely want the angle itself. The projection parameter used
 * throughout the queries layer divides by the target's own squared length
 *
 *     t = Dot(u, v) / Dot(v, v)   ->  u projected onto v, in units of v
 *
 * which is precisely how ClosestPointOnLine and ClosestPointOnSegment locate the
 * perpendicular foot without ever computing a cosine or a square root.
 *
 * --------------------------------------------------------------------------------------
 *
 * @param v1 First vector to compute the dot product with
 * @param v2 Second vector to compute the dot product with
 * @return Dot product: v1 · v2 = v1.x * v2.x + v1.y * v2.y
 *
 * @see https://en.wikipedia.org/wiki/Dot_product#Geometric_definition
 */
template <ScalarType T>
[[nodiscard]] inline T Dot(const Vector2<T>& v1, const Vector2<T>& v2) {
    return v1.x * v2.x + v1.y * v2.y;
}

/**
 * @brief Computes the "2D" cross product (vector product) of two vectors.
 *
 * --------------------------------------------------------------------------------------
 *
 * Not a vector, but the scalar z-component of the equivalent 3D cross product.
 * Magnitude gives the area of the parallelogram spanned by the two vectors;
 * sign gives which side of v1 that v2 falls on.
 *
 * Case 1: v1 × v2 > 0  (v2 is CCW from v1 / to its left)
 *
 *              *---------------* v1 + v2
 *             /               /
 *        v2  /       +       /        sweeping v1 towards v2 turns CCW
 *           /               /
 *          *------->-------*
 *       origin     v1
 *
 * Case 2: v1 × v2 < 0  (v2 is CW from v1 / to its right)
 *
 *       origin     v1
 *          *------->-------*
 *           \               \
 *        v2  \       -       \        sweeping v1 towards v2 turns CW
 *             \               \
 *              *---------------* v1 + v2
 *
 * --------------------------------------------------------------------------------------
 *
 * ** WHY CROSS PRODUCT EXPRESSION MEASURES THE AREA? **
 *
 * - Rotating v1 by 90° CCW gives the perpendicular vector
 *      perp(v1) = (-v1.y, v1.x), and Dot(perp(v1), v2) = v1.x * v2.y - v1.y * v2.x,
 *   which is the definition of the cross product.
 * - So the cross product is |v1| times the component of v2 perpendicular to v1 ->
 *   base × height -> area of the parallelogram.
 *
 * --------------------------------------------------------------------------------------
 *
 * @param v1 First vector to compute the cross product with
 * @param v2 Second vector to compute the cross product with
 * @return Cross product: v1 × v2 = v1.x * v2.y - v1.y * v2.x
 *
 * @see https://en.wikipedia.org/wiki/Cross_product#Two_dimensions
 */
template <ScalarType T>
[[nodiscard]] inline T Cross(const Vector2<T>& v1, const Vector2<T>& v2) {
    // z-component of the equivalent 3D cross product
    return v1.x * v2.y - v1.y * v2.x;
}

/**
 * @brief Computes the squared length (magnitude squared) of the vector.
 *
 * Prefer this over Length whenever the value is only compared against another
 * distance or used as a denominator -> it skips the square root entirely.
 * The queries layer is built on that preference: the projection parameter
 * divides by LengthSquared(direction), and PointToSegmentDistanceSquared exists
 * so callers can rank distances without ever taking a root.
 *
 * @param v Vector to compute the squared length of
 * @return Squared length: ||v||^2 = v.x^2 + v.y^2 = Dot(v, v)
 */
template <ScalarType T>
[[nodiscard]] inline T LengthSquared(const Vector2<T>& v) {
    return Dot(v, v);
}

/**
 * @brief Computes the length (magnitude) of the vector.
 *
 * @param v Vector to compute the length of
 * @return Length: ||v|| = sqrt(||v||^2) = sqrt(Dot(v, v))
 */
template <ScalarType T>
[[nodiscard]] inline T Length(const Vector2<T>& v) {
    return std::sqrt(LengthSquared(v));
}

/**
 * @brief Computes the normalized (unit) vector.
 *
 * @param v Vector to normalize
 * @return Normalized vector: v / ||v|| = v / Length(v), or {0, 0} if @p v is
 *         within tolerance of zero length
 */
template <ScalarType T>
[[nodiscard]] inline Vector2<T> Normalize(const Vector2<T>& v) {
    const T length{Length(v)};
    if (RobustSign(length) == 0) {
        return Vector2<T>{0, 0};
    }
    return Vector2<T>{v.x / length, v.y / length};
}

//---------------------------------------------------------------------------
// Point / Vector algebra —> respects the affine distinction
//---------------------------------------------------------------------------

/**
 * @brief Subtracts two points to get the vector between them.
 *
 * The one combination of two points the affine model permits: a difference of
 * positions is a displacement, not a position.
 *
 *          * a
 *         ^
 *        /
 *       /        a - b  ->  the displacement carrying b to a
 *      /
 *     * b
 *
 *   so that  b + (a - b) == a
 *
 * There is deliberately no operator+(Point2, Point2). Slide the origin by t and
 * every point becomes p + t, under which a - b is unchanged but a + b picks up
 * 2t instead of t -- so a sum of positions names no fixed location, while their
 * difference does.
 *
 * @param a Head point (the destination)
 * @param b Tail point (where the displacement starts)
 * @return Vector2 from b to a: a - b
 */
template <ScalarType T>
[[nodiscard]] inline Vector2<T> operator-(const Point2<T>& a, const Point2<T>& b) {
    return Vector2<T>{a.x - b.x, a.y - b.y};
}

/**
 * @brief Adds a vector to a point to translate the point.
 *
 *            * p + v
 *           ^
 *          /
 *     v   /      p slid along the displacement v
 *        /
 *       * p
 *
 *   Inverse of point subtraction:  (p + v) - p == v
 *
 * @param p Point to translate
 * @param v Vector to translate the point with
 * @return Translated point: p + v
 */
template <ScalarType T>
[[nodiscard]] inline Point2<T> operator+(const Point2<T>& p, const Vector2<T>& v) {
    return Point2<T>{p.x + v.x, p.y + v.y};
}

/**
 * @brief Adds a point to a vector to translate the point (commutative form of Point2 + Vector2).
 *
 * Provided so translation reads naturally in either order. The asymmetry that
 * matters is preserved either way: point + vector and vector + point both give
 * a point, vector + vector gives a vector, and point + point is not defined.
 *
 * @param v Vector to translate the point with
 * @param p Point to translate
 * @return Translated point: v + p = p + v
 */
template <ScalarType T>
[[nodiscard]] inline Point2<T> operator+(const Vector2<T>& v, const Point2<T>& p) {
    return p + v;
}

/**
 * @brief Adds two vectors to get the resulting vector.
 *
 * @param a First vector
 * @param b Second vector
 * @return Resulting vector: a + b
 */
template <ScalarType T>
[[nodiscard]] inline Vector2<T> operator+(const Vector2<T>& a, const Vector2<T>& b) {
    return Vector2<T>{a.x + b.x, a.y + b.y};
}

/**
 * @brief Subtracts two vectors to get the resulting vector.
 *
 * @param a First vector
 * @param b Second vector
 * @return Resulting vector: a - b
 */
template <ScalarType T>
[[nodiscard]] inline Vector2<T> operator-(const Vector2<T>& a, const Vector2<T>& b) {
    return Vector2<T>{a.x - b.x, a.y - b.y};
}

/**
 * @brief Multiplies a vector by a scalar to get the resulting vector.
 *
 * @param v Vector to multiply
 * @param s Scalar to multiply the vector by
 * @return Resulting vector: v * s
 */
template <ScalarType T>
[[nodiscard]] inline Vector2<T> operator*(const Vector2<T>& v, T s) {
    return Vector2<T>{v.x * s, v.y * s};
}

/**
 * @brief Multiplies a scalar by a vector to get the resulting vector.
 *
 * @param s Scalar to multiply the vector by
 * @param v Vector to multiply
 * @return Resulting vector: s * v
 */
template <ScalarType T>
[[nodiscard]] inline Vector2<T> operator*(T s, const Vector2<T>& v) {
    return v * s;
}

}  // namespace geometry_kernel::core
