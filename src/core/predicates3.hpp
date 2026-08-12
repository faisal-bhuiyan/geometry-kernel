#pragma once

#include "area.hpp"
#include "point3.hpp"
#include "tolerance.hpp"
#include "types.hpp"
#include "vector3.hpp"

namespace geometry_kernel::core {

//---------------------------------------------------------------------------
// Triangle predicates (3D)
//---------------------------------------------------------------------------

/**
 * @brief Signed half-space test for edge (vi, vj), relative to @p point, oriented by @p normal.
 *
 * The 3D generalization of SignedTriangleAreaTimes2's sign test. In 2D, the cross product of
 * an edge and a point-relative vector is itself a scalar with a usable sign. In 3D,
 * Cross(vj - vi, point - vi) is a full vector -> perpendicular to both the edge and
 * (point - vi) -> so it has no sign on its own. Dotting it with the face @p normal
 * projects out a scalar that answers "which side of this edge, as seen from the
 * normal direction, does point fall on" -- i.e. a half-space test for the plane that
 * contains the edge and is oriented by normal.
 *
 * This needs no explicit projection of @p point onto the triangle's plane first: any
 * component of (point - vi) parallel to normal is perpendicular to the in-plane edge,
 * so Cross(edge, that component) is itself perpendicular to normal and contributes
 * zero to the dot product. The result depends only on point's in-plane position.
 *
 * Viewed looking against @p normal (i.e. normal points out of the page, toward the
 * reader), the directed edge vi -> vj divides the plane into two half-spaces:
 *
 * Case 1: point left of directed edge vi -> vj -> positive
 *
 *              * point
 *
 *     vi ----------> vj
 *
 * EdgeSign3 > 0
 *
 * Case 2: point right of directed edge vi -> vj -> negative
 *
 *     vi ----------> vj
 *
 *              * point
 *
 * EdgeSign3 < 0
 *
 * Case 3: point collinear with the edge (in-plane) -> zero
 *
 *     vi ---- * ----> vj
 *           point
 *
 * EdgeSign3 == 0 (up to RobustSign / tolerance)
 *
 * Height off the triangle's plane does not change the sign: only the in-plane
 * side of the edge matters. That is why PointInTriangle (3D) can test a point
 * that is not coplanar with the face without projecting first.
 *
 * @param vi Edge start vertex
 * @param vj Edge end vertex
 * @param point Point to test
 * @param normal Face normal orienting the half-space test (need not be unit length --
 *        only its sign/direction matters, NOT its magnitude)
 * @return Signed scalar -> sign indicates which side of the edge (as seen from normal)
 *         point falls on. Magnitude is not independently meaningful (unlike
 *         SignedTriangleAreaTimes2, this is not twice a triangle area -> it is scaled by
 *         |normal| and is only useful through RobustSign()).
 */
template <ScalarType T>
[[nodiscard]] inline T EdgeSign3(
    const Point3<T>& vi, const Point3<T>& vj, const Point3<T>& point, const Vector3<T>& normal
) {
    // Scalar triple product [edge, point-vi, normal] -> signed volume of the
    // parallelepiped they span -> sign = which side of the edge the point is on
    return Dot(Cross(vj - vi, point - vi), normal);
}

}  // namespace geometry_kernel::core
