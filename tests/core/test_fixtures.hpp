#pragma once

#include <vector>

#include "core/point.hpp"
#include "core/types.hpp"

namespace geometry_kernel::test {

using namespace geometry_kernel::core;

/**
 * @brief Casts a `Point2D` (`Point<double>`) polygon to any `Point<T>` precision.
 *
 * Allows the `double`-precision fixture constants (kCcwSquare, kCcwLShape, etc.)
 * to be reused in typed tests parameterised over `float`, `double`, and `long double`
 * without duplicating fixture data.
 *
 * @tparam T Target scalar type satisfying `ScalarType`.
 * @param in Source polygon in `double` precision.
 * @return New polygon with every coordinate cast to `T`.
 */
template <ScalarType T>
[[nodiscard]] std::vector<Point<T>> CastPolygon(const std::vector<Point2D>& in) {
    std::vector<Point<T>> out;
    out.reserve(in.size());
    for (const auto& p : in) {
        out.push_back(Point<T>{static_cast<T>(p.x), static_cast<T>(p.y)});
    }
    return out;
}

/*
 * CCW Triangle
 *
 *        v3 (0,1)
 *          *
 *         / \
 *        /   \
 *       /     \
 *      *-------*
 *  v1 (0,0)   v2 (1,0)
 *
 *  Order: v1 -> v2 -> v3
 *  All vertices convex
 *  All vertices are ear tips
 */
static const std::vector<Point2D> kCcwTriangle = {{
    Point2D{0., 0.},  // v1
    Point2D{1., 0.},  // v2
    Point2D{0., 1.}   // v3
}};

/*
 * CW Triangle: same shape as CCW triangle, but reversed winding
 *
 *        v2 (0,1)
 *          *
 *         / \
 *        /   \
 *       /     \
 *      *-------*
 *  v1 (0,0)   v3 (1,0)
 *
 *  Order: v1 -> v2 -> v3
 *  All vertices reflex
 *  No ear tips
 */
static const std::vector<Point2D> kCwTriangle = {{
    Point2D{0., 0.},  // v1
    Point2D{0., 1.},  // v2
    Point2D{1., 0.}   // v3
}};

/*
 * CCW Square
 *
 *  v4 (0,1) *-------* v3 (1,1)
 *           |       |
 *           |       |
 *  v1 (0,0) *-------* v2 (1,0)
 *
 *  Order: v1 -> v2 -> v3 -> v4
 *  All vertices convex
 *  All vertices are ear tips
 */
static const std::vector<Point2D> kCcwSquare = {{
    Point2D{0., 0.},  // v1
    Point2D{1., 0.},  // v2
    Point2D{1., 1.},  // v3
    Point2D{0., 1.}   // v4
}};

/*
 * CW Square: same shape as CCW square, but reversed winding
 *
 *  v2 (0,1) *-------* v3 (1,1)
 *           |       |
 *           |       |
 *  v1 (0,0) *-------* v4 (1,0)
 *
 *  Order: v1 -> v2 -> v3 -> v4
 *  All vertices reflex
 *  No ear tips
 */
static const std::vector<Point2D> kCwSquare = {{
    Point2D{0., 0.},  // v1
    Point2D{0., 1.},  // v2
    Point2D{1., 1.},  // v3
    Point2D{1., 0.}   // v4
}};

/*
 * CCW L-Shape (concave polygon)
 *
 *  v6 (0,2) *-------* v5 (1,2)
 *           |       |
 *           |       |
 *           |       | v4 (1,1)  <-- reflex vertex
 *           |       *---------* v3 (2,1)
 *           |                 |
 *           |                 |
 *  v1 (0,0) *-----------------* v2 (2,0)
 *
 *  Order: v1 -> v2 -> v3 -> v4 -> v5 -> v6
 *  Reflex vertex: v4
 *  Five ear tips
 */
static const std::vector<Point2D> kCcwLShape = {{
    Point2D{0., 0.},  // v1
    Point2D{2., 0.},  // v2
    Point2D{2., 1.},  // v3
    Point2D{1., 1.},  // v4
    Point2D{1., 2.},  // v5
    Point2D{0., 2.}   // v6
}};

/*
 * CCW M-Shape
 *
 *   v5             v3
 * (0,3)           (4,3)
 *  *               *
 *  | \           / |
 *  |  \         /  |
 *  |   \   v4  /   |
 *  |     (2,1)     |
 *  |               |
 *  |               |
 *  |               |
 *  * ------------- *
 * v1 (0,0)      v2 (4,0)
 *
 *  Order: v1 -> v2 -> v3 -> v4 -> v5
 *  Reflex vertex: v4
 */
static const std::vector<Point2D> kCcwMShape = {
    Point2D{0., 0.},  // v1
    Point2D{4., 0.},  // v2
    Point2D{4., 3.},  // v3
    Point2D{2., 1.},  // v4
    Point2D{0., 3.}   // v5
};

// Contains a convex vertex that is not an ear - reflex vertex lies inside triangle of convex vertex
static const std::vector<Point2D> kReflexVertexInsideEarTriangle = kCcwMShape;

static const std::vector<Point2D> kCcwStarShape = {
    Point2D{3., 0.},  // v1
    Point2D{4., 1.},  // v2
    Point2D{5., 0.},  // v3
    Point2D{4., 2.},  // v4
    Point2D{5., 4.},  // v5
    Point2D{3., 3.},  // v6
    Point2D{1., 4.},  // v7
    Point2D{2., 2.},  // v8
    Point2D{1., 0.},  // v9
    Point2D{2., 1.}   // v10
};

/*
 * kCollinearShape Interior Angle
 *
 *           v5(0,1) *-------* v4 (2,1)
 *                           |
 *                           |
 *  v1 (0,0) *-------*-------* v3 (2,0)
 *              v2 (1,0)  <-- kCollinearShape (180 degrees)
 *
 *  Order: v1 -> v2 -> v3 -> v4 -> v5
 *
 *  NOTE: v2 treated as reflex vertex by PolygonMesh logic -> collinear vertices need
 *  to be filtered out upstream
 */
static const std::vector<Point2D> kCollinearShape = {{
    Point2D{0., 0.},  // v1
    Point2D{1., 0.},  // v2
    Point2D{2., 0.},  // v3
    Point2D{2., 1.},  // v4
    Point2D{0., 1.}   // v5
}};

// Narrow wedge: after clipping one ear, a neighbor may become reflex
static const std::vector<Point2D> kNarrowWedge = {
    Point2D{0., 0.},    // v1
    Point2D{10., 0.},   // v2
    Point2D{10., 0.5},  // v3
    Point2D{9., 0.05},  // v4
    Point2D{0., 1.}     // v5
};

/*
 * CCW simple concave polygon — based on simple_concave_poly.csv
 * Expected area = 0.65
 */
static const std::vector<Point2D> kSimpleConcavePolygon = {{
    Point2D{0., 0.},    // v1
    Point2D{0., 1.},    // v2
    Point2D{0.5, 0.5},  // v3
    Point2D{1., 1.},    // v4
    Point2D{0.8, 0.8},  // v5
    Point2D{1., 0.}     // v6
}};

}  // namespace geometry_kernel::test
