#pragma once

#include <cmath>

#include "point.hpp"
#include "tolerance.hpp"
#include "types.hpp"

namespace geometry_kernel::core {

//---------------------------------------------------------------------------
// Vector in 3D space
//---------------------------------------------------------------------------

/// @brief 3D vector in space
template <ScalarType T>
struct Vector3 {
    T x{};
    T y{};
    T z{};
};

/// @brief 3D vector in space with double precision
using Vector3D = Vector3<double>;

/**
 * @brief Computes the dot product (scalar product) of two vectors
 *
 * @param v1 First vector to compute the dot product with
 * @param v2 Second vector to compute the dot product with
 * @return Dot product: v1 · v2 = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z
 */
template <ScalarType T>
[[nodiscard]] inline T Dot(const Vector3<T>& v1, const Vector3<T>& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

/**
 * @brief Computes the cross product (vector product) of two vectors.
 *
 * Unlike the 2D Cross (which returns the scalar z-component only), the 3D cross
 * product is itself a vector: perpendicular to both v1 and v2, with magnitude equal
 * to the area of the parallelogram they span, and direction given by the right-hand
 * rule.
 *
 * @param v1 First vector to compute the cross product with
 * @param v2 Second vector to compute the cross product with
 * @return Cross product: v1 × v2
 */
template <ScalarType T>
[[nodiscard]] inline Vector3<T> Cross(const Vector3<T>& v1, const Vector3<T>& v2) {
    return Vector3<T>{
        v1.y * v2.z - v1.z * v2.y,  // x-component
        v1.z * v2.x - v1.x * v2.z,  // y-component
        v1.x * v2.y - v1.y * v2.x   // z-component
    };
}

/**
 * @brief Computes the squared length (magnitude squared) of the vector.
 *
 * @param v Vector to compute the squared length of
 * @return Squared length: ||v||^2 = Dot(v, v)
 */
template <ScalarType T>
[[nodiscard]] inline T LengthSquared(const Vector3<T>& v) {
    return Dot(v, v);
}

/**
 * @brief Computes the length (magnitude) of the vector.
 *
 * @param v Vector to compute the length of
 * @return Length: ||v|| = sqrt(LengthSquared(v))
 */
template <ScalarType T>
[[nodiscard]] inline T Length(const Vector3<T>& v) {
    return std::sqrt(LengthSquared(v));
}

/**
 * @brief Computes the normalized (unit) vector.
 *
 * @param v Vector to normalize
 * @return Normalized vector: v / Length(v). Returns the zero vector if @p v is
 *         within tolerance of zero length (degenerate)
 */
template <ScalarType T>
[[nodiscard]] inline Vector3<T> Normalize(const Vector3<T>& v) {
    const T length{Length(v)};
    if (RobustSign(length) == 0) {
        return Vector3<T>{0, 0, 0};
    }
    return Vector3<T>{v.x / length, v.y / length, v.z / length};
}

//---------------------------------------------------------------------------
// Point / Vector algebra -> respects the affine distinction
//---------------------------------------------------------------------------

/**
 * @brief Subtracts two points to get the vector between them: a - b
 *
 * @param a First point
 * @param b Second point
 * @return Vector between the two points: a - b
 */
template <ScalarType T>
[[nodiscard]] inline Vector3<T> operator-(const Point3<T>& a, const Point3<T>& b) {
    return Vector3<T>{a.x - b.x, a.y - b.y, a.z - b.z};
}

/**
 * @brief Adds a vector to a point to translate the point: p + v
 *
 * @param p Point to translate
 * @param v Vector to translate the point with
 * @return Translated point: p + v
 */
template <ScalarType T>
[[nodiscard]] inline Point3<T> operator+(const Point3<T>& p, const Vector3<T>& v) {
    return Point3<T>{p.x + v.x, p.y + v.y, p.z + v.z};
}

/**
 * @brief Adds a point to a vector to translate the point (commutative form of Point3 + Vector3)
 *
 * @param v Vector to translate the point with
 * @param p Point to translate
 * @return Translated point: v + p = p + v
 */
template <ScalarType T>
[[nodiscard]] inline Point3<T> operator+(const Vector3<T>& v, const Point3<T>& p) {
    return p + v;
}

/**
 * @brief Adds two vectors to get the resulting vector: a + b
 *
 * @param a First vector
 * @param b Second vector
 * @return Resulting vector: a + b
 */
template <ScalarType T>
[[nodiscard]] inline Vector3<T> operator+(const Vector3<T>& a, const Vector3<T>& b) {
    return Vector3<T>{a.x + b.x, a.y + b.y, a.z + b.z};
}

/**
 * @brief Subtracts two vectors to get the resulting vector: a - b
 *
 * @param a First vector
 * @param b Second vector
 * @return Resulting vector: a - b
 */
template <ScalarType T>
[[nodiscard]] inline Vector3<T> operator-(const Vector3<T>& a, const Vector3<T>& b) {
    return Vector3<T>{a.x - b.x, a.y - b.y, a.z - b.z};
}

/**
 * @brief Multiplies a vector by a scalar: v * s
 *
 * @param v Vector to multiply
 * @param s Scalar to multiply the vector by
 * @return Resulting vector: v * s
 */
template <ScalarType T>
[[nodiscard]] inline Vector3<T> operator*(const Vector3<T>& v, T s) {
    return Vector3<T>{v.x * s, v.y * s, v.z * s};
}

/**
 * @brief Multiplies a scalar by a vector: s * v
 *
 * @param s Scalar to multiply the vector by
 * @param v Vector to multiply
 * @return Resulting vector: s * v
 */
template <ScalarType T>
[[nodiscard]] inline Vector3<T> operator*(T s, const Vector3<T>& v) {
    return v * s;
}

}  // namespace geometry_kernel::core
