#pragma once

#include <algorithm>
#include <span>
#include <stdexcept>
#include <vector>

#include "core/area.hpp"
#include "core/point2.hpp"
#include "core/types.hpp"
#include "polygon_mesh.hpp"

namespace geometry_kernel::algorithms {

using namespace geometry_kernel::core;

//---------------------------------------------------------------------------
// Result types
//---------------------------------------------------------------------------

/// @brief Indices into a `TriangulationResult::vertices` array for one output triangle.
struct Triangle {
    std::size_t v0{0};  ///< first vertex index
    std::size_t v1{0};  ///< second vertex index
    std::size_t v2{0};  ///< third vertex index
};

/**
 * @brief Output of the ear-clipping triangulator.
 *
 * `vertices` is the CCW-normalised vertex list after removing any duplicate
 * closing vertex. Triangle indices reference positions in this list, not in
 * the original input.
 */
template <ScalarType T>
struct TriangulationResult {
    std::vector<Point2<T>> vertices;  ///< CCW-normalised input vertices
    std::vector<Triangle> triangles;  ///< Output triangles in clip order
    T total_area{};                   ///< Sum of all triangle areas
};

//---------------------------------------------------------------------------
// Ear Clipping Algorithm
//---------------------------------------------------------------------------

/**
 * @brief Triangulates a simple polygon using the ear-clipping algorithm.
 *
 * CW-wound input is silently reversed to CCW before triangulation.
 *
 * @param polygon_vertices Vertices of the polygon as an open ring (no closing vertex).
 * @return Triangulation result with CCW-normalised vertices and output triangles.
 * @throws std::invalid_argument If fewer than 3 vertices are provided.
 * @throws std::runtime_error    If no ear is found (indicates a non-simple or degenerate polygon).
 * @pre Vertices form a simple polygon (no self-intersections).
 * @pre No duplicate consecutive vertices or collinear triples.
 * @pre The last vertex does not duplicate the first (no explicit closing vertex).
 */
template <ScalarType T>
inline TriangulationResult<T> TriangulatePolygonWithEarClipping(
    std::span<const Point2<T>> polygon_vertices
) {
    // Copy the input vertices to a local vector
    std::vector<Point2<T>> vertices{polygon_vertices.begin(), polygon_vertices.end()};
    if (vertices.size() < 3U) {
        throw std::invalid_argument("TriangulatePolygonWithEarClipping: need at least 3 vertices");
    }

    // Reverse CW-wound input to CCW - O(n) winding check is acceptable since triangulation itself is
    // O(n^2)
    if (SignedPolygonArea(vertices) < T{}) {
        std::reverse(vertices.begin(), vertices.end());
    }

    TriangulationResult<T> result{};
    result.vertices = vertices;                // copy assignment
    PolygonMesh<T> mesh{std::move(vertices)};  // move constructor

    // triangulation -> iteratively clip ears until three vertices remain
    while (mesh.ActiveVertexCount() > 3U) {
        if (!mesh.HasEarVertex()) {
            throw std::runtime_error(
                "TriangulatePolygonWithEarClipping: no ear found - polygon may be non-simple"
            );
        }
        const std::size_t ear{mesh.EarHeadIndex()};
        const auto [v0, v1, v2] = mesh.EarVertexNeighbors(ear);
        result.triangles.emplace_back(Triangle{v0, v1, v2});
        result.total_area +=
            TriangleArea(result.vertices[v0], result.vertices[v1], result.vertices[v2]);
        mesh.ClipEarVertex(ear);
    }

    // Add final three active vertices to form the last triangle
    const std::vector<std::size_t> remaining{mesh.ActiveBoundaryVertices()};
    if (remaining.size() != 3U) {
        throw std::runtime_error(
            "TriangulatePolygonWithEarClipping: expected 3 vertices in final polygon"
        );
    }
    result.triangles.emplace_back(Triangle{remaining[0], remaining[1], remaining[2]});
    result.total_area += TriangleArea(
        result.vertices[remaining[0]], result.vertices[remaining[1]], result.vertices[remaining[2]]
    );

    return result;
}

/// @brief Convenience overload -> accepts std::vector<Point2<T>> directly.
template <ScalarType T>
inline TriangulationResult<T> TriangulatePolygonWithEarClipping(
    const std::vector<Point2<T>>& polygon_vertices
) {
    return TriangulatePolygonWithEarClipping(std::span<const Point2<T>>{polygon_vertices});
}

}  // namespace geometry_kernel::algorithms
