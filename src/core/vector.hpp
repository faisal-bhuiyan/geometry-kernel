#pragma once

#include <cmath>

#include "point.hpp"
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
 * @brief Computes the dot product (scalar product) of two vectors
 *
 * @param v1 First vector to compute the dot product with.
 * @param v2 Second vector to compute the dot product with.
 * @return Dot product: v1 · v2 = v1.x * v2.x + v1.y * v2.y
 */
template <ScalarType T>
[[nodiscard]] inline T Dot(const Vector2<T>& v1, const Vector2<T>& v2) {
    return v1.x * v2.x + v1.y * v2.y;
}

/**
 * @brief Computes the "2D" cross product (vector product) of two vectors.
 *
 * Not a vector, but the scalar z-component of the equivalent 3D cross product.
 * Sign gives orientation -> magnitude gives the area of the parallelogram spanned by the two
 * vectors.
 *
 * @param v1 First vector to compute the cross product with.
 * @param v2 Second vector to compute the cross product with.
 * @return Cross product: v1 × v2 = v1.x * v2.y - v1.y * v2.x
 */
template <ScalarType T>
[[nodiscard]] inline T Cross(const Vector2<T>& v1, const Vector2<T>& v2) {
    // z-component of the equivalent 3D cross product
    return v1.x * v2.y - v1.y * v2.x;
}

/**
 * @brief Computes the squared length (magnitude squared) of the vector.
 *
 * @param v Vector to compute the squared length of.
 * @return Squared length: ||v||^2 = v.x^2 + v.y^2 = Dot(v, v)
 */
template <ScalarType T>
[[nodiscard]] inline T LengthSquared(const Vector2<T>& v) {
    return Dot(v, v);
}

/**
 * @brief Computes the length (magnitude) of the vector.
 *
 * @param v Vector to compute the length of.
 * @return Length: ||v|| = sqrt(||v||^2) = sqrt(Dot(v, v))
 */
template <ScalarType T>
[[nodiscard]] inline T Length(const Vector2<T>& v) {
    return std::sqrt(LengthSquared(v));
}

/**
 * @brief Computes the normalized (unit) vector.
 *
 * @param v Vector to normalize.
 * @return Normalized vector: v / ||v|| = v / Length(v)
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
 * @param a First point.
 * @param b Second point.
 *
 *        a
 *       /
 *      /  a - b  (vector from b -> a)
 *     /
 *    b
 *
 * @return Vector2: a - b
 */
template <ScalarType T>
[[nodiscard]] inline Vector2<T> operator-(const Point2<T>& a, const Point2<T>& b) {
    return Vector2<T>{a.x - b.x, a.y - b.y};
}

/**
 * @brief Adds a vector to a point to translate the point.
 *
 * @param p Point to translate.
 * @param v Vector to translate the point with.
 * @return Translated point: p + v
 */
template <ScalarType T>
[[nodiscard]] inline Point2<T> operator+(const Point2<T>& p, const Vector2<T>& v) {
    return Point2<T>{p.x + v.x, p.y + v.y};
}

/**
 * @brief Adds a point to a vector to translate the point (commutative form of Point2 + Vector2).
 *
 * @param v Vector to translate the point with.
 * @param p Point to translate.
 * @return Translated point: v + p = p + v
 */
template <ScalarType T>
[[nodiscard]] inline Point2<T> operator+(const Vector2<T>& v, const Point2<T>& p) {
    return p + v;
}

/**
 * @brief Adds two vectors to get the resulting vector.
 *
 * @param a First vector.
 * @param b Second vector.
 * @return Resulting vector: a + b
 */
template <ScalarType T>
[[nodiscard]] inline Vector2<T> operator+(const Vector2<T>& a, const Vector2<T>& b) {
    return Vector2<T>{a.x + b.x, a.y + b.y};
}

/**
 * @brief Subtracts two vectors to get the resulting vector.
 *
 * @param a First vector.
 * @param b Second vector.
 * @return Resulting vector: a - b
 */
template <ScalarType T>
[[nodiscard]] inline Vector2<T> operator-(const Vector2<T>& a, const Vector2<T>& b) {
    return Vector2<T>{a.x - b.x, a.y - b.y};
}

/**
 * @brief Multiplies a vector by a scalar to get the resulting vector.
 *
 * @param v Vector to multiply.
 * @param s Scalar to multiply the vector by.
 * @return Resulting vector: v * s
 */
template <ScalarType T>
[[nodiscard]] inline Vector2<T> operator*(const Vector2<T>& v, T s) {
    return Vector2<T>{v.x * s, v.y * s};
}

/**
 * @brief Multiplies a scalar by a vector to get the resulting vector.
 *
 * @param s Scalar to multiply the vector by.
 * @param v Vector to multiply.
 * @return Resulting vector: s * v
 */
template <ScalarType T>
[[nodiscard]] inline Vector2<T> operator*(T s, const Vector2<T>& v) {
    return v * s;
}

}  // namespace geometry_kernel::core
