#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <set>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "algorithms/polygon_mesh.hpp"
#include "core/area.hpp"
#include "core/point2.hpp"
#include "core/test_fixtures.hpp"
#include "core/types.hpp"

namespace geometry_kernel::test {

using namespace geometry_kernel::core;
using namespace geometry_kernel::algorithms;

//------------------------------------------------------------------------------
// Invariant helpers
//------------------------------------------------------------------------------

/**
 * Verifies that the boundary ring contains exactly `expected_active` active vertices,
 * each appearing exactly once (no cycles or dropped nodes), and that the forward
 * and backward links between consecutive entries are consistent:
 *
 *   boundary[i].polygon_next        == boundary[i+1]
 *   boundary[i+1].polygon_previous  == boundary[i]   (for all i, wrapping around)
 *
 * This is the primary structural invariant for the cyclic boundary ring.
 * Called after construction and after every `ClipEarVertex` to ensure the ring
 * stays coherent.
 */
template <ScalarType T>
static void ExpectBoundaryRingValid(const PolygonMesh<T>& mesh, std::size_t expected_active) {
    const auto boundary = mesh.ActiveBoundaryVertices();
    ASSERT_EQ(boundary.size(), expected_active);
    // If there are no active vertices, the boundary ring is empty and valid
    if (expected_active == 0) {
        return;
    }

    // Track visited vertices to ensure each vertex appears exactly once
    std::set<std::size_t> visited;
    for (std::size_t index : boundary) {
        ASSERT_TRUE(mesh.Node(index).is_active);
        ASSERT_TRUE(visited.insert(index).second) << "duplicate index in boundary walk";
    }
    ASSERT_EQ(visited.size(), expected_active);

    // Verify that the forward and backward links are consistent for each vertex
    for (std::size_t i = 0; i < boundary.size(); ++i) {
        const std::size_t current = boundary[i];
        const std::size_t next = boundary[(i + 1) % boundary.size()];
        EXPECT_EQ(mesh.Node(current).polygon_next, next);
        EXPECT_EQ(mesh.Node(next).polygon_previous, current);
    }
}

/**
 * Asserts that no active vertex has a `kNone` boundary link.
 * Every vertex that is still part of the polygon ring must have both `polygon_previous`
 * and `polygon_next` pointing to valid neighbors. Called after construction and after
 * each `ClipEarVertex` to catch dangling links caused by incomplete boundary relinking.
 */
template <ScalarType T>
static void ExpectNoDeadLinksInActiveBoundary(const PolygonMesh<T>& mesh, std::size_t num_vertices) {
    for (std::size_t i = 0; i < num_vertices; ++i) {
        // Skip inactive vertices
        if (!mesh.Node(i).is_active) {
            continue;
        }
        // Assert that the vertex has valid previous and next neighbors
        EXPECT_NE(mesh.Node(i).polygon_previous, kNone);
        EXPECT_NE(mesh.Node(i).polygon_next, kNone);
    }
}

/**
 * Asserts the fundamental ear-tip invariant: an ear tip MUST be convex.
 * A reflex vertex can never be an ear because its interior angle exceeds
 * 180 deg, so the "triangle" it would form would fold back into the polygon.
 * Violation here indicates a bug in IsEarVertex or ReclassifyVertex.
 */
template <ScalarType T>
static void ExpectEarImpliesConvex(const PolygonMesh<T>& mesh, std::size_t num_vertices) {
    for (std::size_t i = 0; i < num_vertices; ++i) {
        if (mesh.Node(i).is_ear) {
            EXPECT_TRUE(mesh.Node(i).is_convex) << "vertex " << i;
        }
    }
}

/**
 * Asserts that ActiveVertexCount() matches the actual number of nodes flagged
 * `is_active == true`.
 * Guards against the active count becoming out of sync with the node flags, e.g.
 * if `ClipEarVertex` decrements the counter incorrectly or fails to set the flag.
 */
template <ScalarType T>
static void ExpectActiveVertexCountMatchesNodes(
    const PolygonMesh<T>& mesh, std::size_t num_vertices
) {
    std::size_t active_count{0};
    for (std::size_t i = 0; i < num_vertices; ++i) {
        if (mesh.Node(i).is_active) {
            ++active_count;
        }
    }
    EXPECT_EQ(mesh.ActiveVertexCount(), active_count);
}

/**
 * Asserts that no reflex vertex is also marked as an ear tip.
 * Convenience wrapper around `ExpectEarImpliesConvex` that checks
 * the contrapositive: `!is_convex -> !is_ear`.
 */
template <ScalarType T>
static void ExpectReflexVertexIsNeverEar(const PolygonMesh<T>& mesh, std::size_t num_vertices) {
    for (std::size_t i = 0; i < num_vertices; ++i) {
        if (!mesh.Node(i).is_convex) {
            EXPECT_FALSE(mesh.Node(i).is_ear) << "vertex " << i;
        }
    }
}

/**
 * Returns the number of vertices classified as reflex (is_convex == false)
 * across the first `num_vertices` nodes, including already-clipped (inactive) ones.
 * Use `CountReflexVertices(mesh, num_vertices)` on a freshly constructed mesh to check
 * the expected number of concavities.
 */
template <ScalarType T>
static std::size_t CountReflexVertices(const PolygonMesh<T>& mesh, std::size_t num_vertices) {
    std::size_t count{0};
    for (std::size_t i = 0; i < num_vertices; ++i) {
        if (!mesh.Node(i).is_convex) {
            ++count;
        }
    }
    return count;
}

/**
 * Returns the number of vertices currently flagged as ear tips
 * across the first `num_vertices` nodes.
 * Used to verify Meisters' theorem (>= 2 ears for `num_vertices >= 4`) and to
 * confirm that convex polygons have all vertices as ears.
 */
template <ScalarType T>
static std::size_t CountEars(const PolygonMesh<T>& mesh, std::size_t num_vertices) {
    std::size_t count{0};
    for (std::size_t i = 0; i < num_vertices; ++i) {
        if (mesh.Node(i).is_ear) {
            ++count;
        }
    }
    return count;
}

/**
 * Walks polygon_next exactly `num_steps` times starting from head and asserts
 * that the traversal returns to head.
 * Verifies that the cyclic boundary ring closes correctly after `num_steps` steps --
 * confirming no broken forward links and no extra or missing nodes.
 */
template <ScalarType T>
static void WalkForwardToHead(const PolygonMesh<T>& mesh, std::size_t head, std::size_t num_steps) {
    std::size_t current = head;
    for (std::size_t step = 0; step < num_steps; ++step) {
        current = mesh.Node(current).polygon_next;
    }
    EXPECT_EQ(current, head);
}

/**
 * Walks polygon_previous exactly `num_steps` times starting from head and asserts
 * that the traversal returns to head.
 * Mirrors `WalkForwardToHead` for the backward direction, ensuring the
 * reverse links are also consistent with a proper cyclic ring.
 */
template <ScalarType T>
static void WalkBackwardToHead(const PolygonMesh<T>& mesh, std::size_t head, std::size_t num_steps) {
    std::size_t current = head;
    for (std::size_t step = 0; step < num_steps; ++step) {
        current = mesh.Node(current).polygon_previous;
    }
    EXPECT_EQ(current, head);
}

/**
 * Scans the first `num_vertices` nodes and returns the indices of all current ear tips.
 * Uses the public `Node(i).is_ear` flag rather than traversing the private
 * ear list, so it works even when the ear head is not directly accessible.
 * Used by ExpectEarListIsCyclic to seed the traversal.
 */
template <ScalarType T>
static std::vector<std::size_t> CollectEarIndices(
    const PolygonMesh<T>& mesh, std::size_t num_vertices
) {
    std::vector<std::size_t> ears;
    for (std::size_t i = 0; i < num_vertices; ++i) {
        if (mesh.Node(i).is_ear) {
            ears.push_back(i);
        }
    }
    return ears;
}

/**
 * Asserts that the ear tip list forms a valid cyclic structure.
 *
 * Collects all ear indices via 'CollectEarIndices()', then starting from the
 * first ear, follows 'ear_next' links until the traversal returns to the
 * starting vertex. Checks that:
 *   - Every node visited has 'is_ear' == true.
 *   - No index appears twice (the cycle is simple).
 *   - The number of visited nodes equals the total ear count.
 *   - 'HasEarVertex()' agrees with whether any ears were found.
 *   - A single-element list is self-referential ('ear_previous' == 'ear_next' == self).
 */
template <ScalarType T>
static void ExpectEarListIsCyclic(const PolygonMesh<T>& mesh, std::size_t num_vertices) {
    const std::vector<std::size_t> ears = CollectEarIndices(mesh, num_vertices);
    if (ears.empty()) {
        EXPECT_FALSE(mesh.HasEarVertex());
        return;
    }

    EXPECT_TRUE(mesh.HasEarVertex());
    for (std::size_t start : ears) {
        std::set<std::size_t> visited;
        std::size_t current = start;
        do {
            ASSERT_TRUE(mesh.Node(current).is_ear);
            ASSERT_TRUE(visited.insert(current).second)
                << "ear list not simple cycle from " << start;
            current = mesh.Node(current).ear_next;
        } while (current != start);
        EXPECT_EQ(visited.size(), ears.size());
        break;  // one component
    }

    // Single-element list: self-referential
    if (ears.size() == 1U) {
        const std::size_t e = ears[0];
        EXPECT_EQ(mesh.Node(e).ear_previous, e);
        EXPECT_EQ(mesh.Node(e).ear_next, e);
    }
}

/**
 * Runs the ear-clipping loop to completion and returns all produced triangles.
 *
 * Repeatedly picks the current EarHeadIndex(), records its three boundary
 * neighbors as a triangle, then clips it. After each clip, calls
 * ExpectBoundaryRingValid, ExpectActiveVertexCountMatchesNodes, and
 * ExpectEarImpliesConvex to assert the mesh stays consistent. When only
 * three active vertices remain, collects them as the final triangle.
 *
 * Throws std::runtime_error if no ear is available while the active count
 * is still > 3, or if the final boundary does not contain exactly 3 vertices.
 * Tests wrap calls to this function in ASSERT_NO_THROW.
 */
template <ScalarType T>
static std::vector<std::array<std::size_t, 3>> TriangulateByEarClipping(
    PolygonMesh<T>& mesh, const std::vector<Point2<T>>& original_positions
) {
    std::vector<std::array<std::size_t, 3>> triangles;
    const std::size_t num_vertices = original_positions.size();
    while (mesh.ActiveVertexCount() > 3) {
        if (!mesh.HasEarVertex()) {
            throw std::runtime_error("TriangulateByEarClipping: no ear while count > 3");
        }
        const std::size_t ear = mesh.EarHeadIndex();
        const auto neigh = mesh.EarVertexNeighbors(ear);
        triangles.push_back({neigh[0], neigh[1], neigh[2]});
        mesh.ClipEarVertex(ear);
        ExpectBoundaryRingValid(mesh, mesh.ActiveVertexCount());
        ExpectActiveVertexCountMatchesNodes(mesh, num_vertices);
        ExpectEarImpliesConvex(mesh, num_vertices);
    }
    const auto boundary = mesh.ActiveBoundaryVertices();
    if (boundary.size() != 3U) {
        throw std::runtime_error("TriangulateByEarClipping: expected 3 boundary vertices");
    }
    triangles.push_back({boundary[0], boundary[1], boundary[2]});
    return triangles;
}

/**
 * Returns the total area of a triangulation by summing individual triangle areas.
 * Each entry in `triangles` holds three vertex indices into `vertices`. Used in
 * integration tests to verify that the sum of clipped triangle areas equals
 * the known area of the original polygon.
 */
template <ScalarType T>
static T SumTriangleAreas(
    const std::vector<std::array<std::size_t, 3>>& triangles, const std::vector<Point2<T>>& vertices
) {
    T sum{};
    for (const auto& t : triangles) {
        sum += TriangleArea(vertices[t[0]], vertices[t[1]], vertices[t[2]]);
    }
    return sum;
}

/**
 * Generates a CCW-wound regular convex polygon with `num_vertices` vertices
 * on a circle of the given `radius`, starting from the bottom (theta = -pi/2).
 * Produced vertices are guaranteed to be strictly convex and in CCW order,
 * making this polygon a valid input for stress and classification tests.
 */
template <ScalarType T>
static std::vector<Point2<T>> RegularConvexPolygon(std::size_t num_vertices, T radius) {
    std::vector<Point2<T>> out;
    out.reserve(num_vertices);
    for (std::size_t i = 0; i < num_vertices; ++i) {
        const T theta = -std::numbers::pi / 2. +
                        (2. * std::numbers::pi * static_cast<T>(i) / static_cast<T>(num_vertices));
        out.push_back({radius * std::cos(theta), radius * std::sin(theta)});
    }
    return out;
}

//------------------------------------------------------------------------------
// Construction tests
//------------------------------------------------------------------------------

TEST(PolygonMeshConstruction, TriangleActiveVertexCount) {
    PolygonMesh<double> mesh(kCcwTriangle);
    EXPECT_EQ(mesh.ActiveVertexCount(), 3U);
}

TEST(PolygonMeshConstruction, TriangleBoundaryRingCyclic) {
    PolygonMesh<double> mesh(kCcwTriangle);
    EXPECT_EQ(mesh.Node(0).polygon_next, 1U);
    EXPECT_EQ(mesh.Node(1).polygon_next, 2U);
    EXPECT_EQ(mesh.Node(2).polygon_next, 0U);
    EXPECT_EQ(mesh.Node(0).polygon_previous, 2U);
    EXPECT_EQ(mesh.Node(1).polygon_previous, 0U);
    EXPECT_EQ(mesh.Node(2).polygon_previous, 1U);
}

TEST(PolygonMeshConstruction, TriangleAllVerticesActive) {
    PolygonMesh<double> mesh(kCcwTriangle);
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_TRUE(mesh.Node(i).is_active);
    }
}

TEST(PolygonMeshConstruction, TriangleAllConvexNoReflex) {
    PolygonMesh<double> mesh(kCcwTriangle);
    EXPECT_EQ(CountReflexVertices(mesh, 3), 0U);
}

TEST(PolygonMeshConstruction, TriangleAllVerticesEarTips) {
    PolygonMesh<double> mesh(kCcwTriangle);
    EXPECT_EQ(CountEars(mesh, 3), 3U);
}

TEST(PolygonMeshConstruction, SquareAllConvexAllEars) {
    PolygonMesh<double> mesh(kCcwSquare);
    EXPECT_EQ(mesh.ActiveVertexCount(), 4U);
    EXPECT_EQ(CountReflexVertices(mesh, 4), 0U);
    EXPECT_EQ(CountEars(mesh, 4), 4U);
}

TEST(PolygonMeshConstruction, LShapeCorrectActiveCount) {
    PolygonMesh<double> mesh(kCcwLShape);
    // L-shape has 6 active vertices
    EXPECT_EQ(mesh.ActiveVertexCount(), 6U);
}

TEST(PolygonMeshConstruction, LShapeCorrectReflexCount) {
    PolygonMesh<double> mesh(kCcwLShape);
    // L-shape has exactly one reflex vertex
    EXPECT_EQ(CountReflexVertices(mesh, 6), 1U);
}

TEST(PolygonMeshConstruction, LShapeAtLeastTwoEarTips) {
    PolygonMesh<double> mesh(kCcwLShape);
    // L-shape has at least two ear tips
    EXPECT_GE(CountEars(mesh, 6), 2U);
}

TEST(PolygonMeshConstruction, ThrowsOnZeroVertices) {
    EXPECT_THROW(PolygonMesh<double>(std::vector<Point2<double>>{}), std::invalid_argument);
}

TEST(PolygonMeshConstruction, ThrowsOnOneVertex) {
    EXPECT_THROW(PolygonMesh<double>(std::vector<Point2<double>>{{0., 0.}}), std::invalid_argument);
}

TEST(PolygonMeshConstruction, ThrowsOnTwoVertices) {
    EXPECT_THROW(
        PolygonMesh<double>(std::vector<Point2<double>>{{0., 0.}, {1., 0.}}), std::invalid_argument
    );
}

//------------------------------------------------------------------------------
// Boundary ring tests
//------------------------------------------------------------------------------

TEST(PolygonMeshBoundaryRing, ForwardTraversalReturnsToHead) {
    // Test that the boundary ring is valid for a triangle
    PolygonMesh<double> triangle(kCcwTriangle);
    WalkForwardToHead(triangle, 0, 3);

    // Test that the boundary ring is valid for a square
    PolygonMesh<double> square(kCcwSquare);
    WalkForwardToHead(square, 0, 4);

    // Test that the boundary ring is valid for an L-shape
    PolygonMesh<double> lshape(kCcwLShape);
    WalkForwardToHead(lshape, 0, 6);
}

TEST(PolygonMeshBoundaryRing, BackwardTraversalReturnsToHead) {
    // Test that the boundary ring is valid for a triangle
    PolygonMesh<double> triangle(kCcwTriangle);
    WalkBackwardToHead(triangle, 0, 3);

    // Test that the boundary ring is valid for a square
    PolygonMesh<double> square(kCcwSquare);
    WalkBackwardToHead(square, 0, 4);

    // Test that the boundary ring is valid for an L-shape
    PolygonMesh<double> lshape(kCcwLShape);
    WalkBackwardToHead(lshape, 0, 6);
}

TEST(PolygonMeshBoundaryRing, NoInvalidLinksForActiveVertices) {
    PolygonMesh<double> lshape(kCcwLShape);
    ExpectNoDeadLinksInActiveBoundary(lshape, 6);
}

TEST(PolygonMeshBoundaryRing, BoundaryRingRemainsValidAfterClip) {
    PolygonMesh<double> square(kCcwSquare);
    square.ClipEarVertex(square.EarHeadIndex());
    ExpectBoundaryRingValid(square, 3);
    ExpectNoDeadLinksInActiveBoundary(square, 4);
}

TEST(PolygonMeshBoundaryRing, ClippedVertexHasKNoneLinks) {
    PolygonMesh<double> square(kCcwSquare);
    const std::size_t ear = square.EarHeadIndex();
    square.ClipEarVertex(ear);
    EXPECT_FALSE(square.Node(ear).is_active);
    EXPECT_EQ(square.Node(ear).polygon_previous, kNone);
    EXPECT_EQ(square.Node(ear).polygon_next, kNone);
}

TEST(PolygonMeshBoundaryRing, HeadAdvancesWhenHeadIsClipped) {
    PolygonMesh<double> square(kCcwSquare);
    ASSERT_EQ(square.Node(0).polygon_previous, 3U);

    // Clip vertex 0 if it is ear otherwise find an ear at index 0 by rotating -> square has all ears
    const std::size_t ear = 0;
    ASSERT_TRUE(square.Node(ear).is_ear);
    const std::size_t head_after_clip = square.Node(ear).polygon_next;
    square.ClipEarVertex(ear);

    // boundary walk should start from head_after_clip
    const auto b = square.ActiveBoundaryVertices();
    ASSERT_FALSE(b.empty());
    EXPECT_EQ(b[0], head_after_clip);
}

//------------------------------------------------------------------------------
// Reflex list tests
//------------------------------------------------------------------------------

TEST(PolygonMeshReflexList, ConvexPolygonHasNoReflexVertices) {
    PolygonMesh<double> square(kCcwSquare);
    EXPECT_EQ(CountReflexVertices(square, 4), 0U);
}

TEST(PolygonMeshReflexList, ConvexVerticesHaveNoReflexListLinks) {
    PolygonMesh<double> square(kCcwSquare);
    for (std::size_t i = 0; i < 4; ++i) {
        ASSERT_TRUE(square.Node(i).is_convex);
        EXPECT_EQ(square.Node(i).reflex_previous, kNone);
        EXPECT_EQ(square.Node(i).reflex_next, kNone);
    }
}

TEST(PolygonMeshReflexList, ReflexVertexAtKnownIndex) {
    PolygonMesh<double> lshape(kCcwLShape);
    // L-shape has exactly one reflex vertex at index 3
    EXPECT_FALSE(lshape.Node(3).is_convex);
}

TEST(PolygonMeshReflexList, CollinearShapeVertexIsClassifiedReflex) {
    PolygonMesh<double> col(kCollinearShape);
    EXPECT_FALSE(col.Node(1).is_convex);
}

TEST(PolygonMeshReflexList, CWPolygonAllVerticesReflex) {
    PolygonMesh<double> cw_sq(kCwSquare);
    EXPECT_EQ(CountReflexVertices(cw_sq, 4), 4U);
}

TEST(PolygonMeshReflexList, ReflexVertexRemovedFromListAfterBecomingConvex) {
    PolygonMesh<double> lshape(kCcwLShape);
    // Empirical: clip ears 5, 4, 0 then vertex 3 becomes convex
    lshape.ClipEarVertex(5);
    lshape.ClipEarVertex(4);
    lshape.ClipEarVertex(0);
    EXPECT_TRUE(lshape.Node(3).is_active);
    EXPECT_TRUE(lshape.Node(3).is_convex);
    ExpectReflexVertexIsNeverEar(lshape, 6);
}

//------------------------------------------------------------------------------
// Ear list tests
//------------------------------------------------------------------------------

TEST(PolygonMeshEarList, HasEarVertexTrueForValidPolygons) {
    EXPECT_TRUE(PolygonMesh<double>(kCcwTriangle).HasEarVertex());
    EXPECT_TRUE(PolygonMesh<double>(kCcwSquare).HasEarVertex());
    EXPECT_TRUE(PolygonMesh<double>(kCcwLShape).HasEarVertex());
}

TEST(PolygonMeshEarList, HasEarVertexFalseForCWPolygon) {
    EXPECT_FALSE(PolygonMesh<double>(kCwSquare).HasEarVertex());
}

TEST(PolygonMeshEarList, EarListCyclicAfterConstruction) {
    PolygonMesh<double> square(kCcwSquare);
    ExpectEarListIsCyclic(square, 4);
}

TEST(PolygonMeshEarList, EarListCyclicAfterOneClip) {
    PolygonMesh<double> square(kCcwSquare);
    square.ClipEarVertex(square.EarHeadIndex());
    ExpectEarListIsCyclic(square, 4);
}

TEST(PolygonMeshEarList, EarListEarImpliesConvex) {
    PolygonMesh<double> lshape(kCcwLShape);
    ExpectEarImpliesConvex(lshape, 6);
}

TEST(PolygonMeshEarList, ActiveEarVerticesHaveValidEarLinks) {
    PolygonMesh<double> square(kCcwSquare);
    for (std::size_t i = 0; i < 4; ++i) {
        if (square.Node(i).is_ear) {
            EXPECT_NE(square.Node(i).ear_previous, kNone);
            EXPECT_NE(square.Node(i).ear_next, kNone);
        }
    }
}

//------------------------------------------------------------------------------
// Ear detection tests
//------------------------------------------------------------------------------

TEST(PolygonMeshEarDetection, ConvexPolygonAllVerticesAreEars) {
    PolygonMesh<double> square(kCcwSquare);
    // Square has all vertices as ears
    EXPECT_EQ(CountEars(square, 4), 4U);
}

TEST(PolygonMeshEarDetection, ReflexVertexIsNeverAnEar) {
    PolygonMesh<double> lshape(kCcwLShape);
    ExpectReflexVertexIsNeverEar(lshape, 6);
    PolygonMesh<double> mshape(kCcwMShape);
    ExpectReflexVertexIsNeverEar(mshape, 5);
}

TEST(PolygonMeshEarDetection, ConcavePolygonHasAtLeastTwoEarTips) {
    PolygonMesh<double> lshape(kCcwLShape);
    EXPECT_GE(CountEars(lshape, 6), 2U);
    PolygonMesh<double> mshape(kCcwMShape);
    EXPECT_GE(CountEars(mshape, 5), 2U);
}

TEST(PolygonMeshEarDetection, ReflexVertexInsideTriangleBlocksEar) {
    PolygonMesh<double> reflex_inside_ear(kReflexVertexInsideEarTriangle);
    EXPECT_TRUE(reflex_inside_ear.Node(0).is_convex);
    EXPECT_FALSE(reflex_inside_ear.Node(0).is_ear);
}

//------------------------------------------------------------------------------
// ClipEarVertex success tests
//------------------------------------------------------------------------------

TEST(PolygonMeshClipEar, ClipDecrementsActiveVertexCount) {
    PolygonMesh<double> square(kCcwSquare);
    const std::size_t before = square.ActiveVertexCount();
    square.ClipEarVertex(square.EarHeadIndex());
    // Active vertex count should decrement by 1
    EXPECT_EQ(square.ActiveVertexCount(), before - 1U);
}

TEST(PolygonMeshClipEar, ClipMarksVertexInactive) {
    PolygonMesh<double> square(kCcwSquare);
    const std::size_t ear = square.EarHeadIndex();
    square.ClipEarVertex(ear);
    // Clipped vertex should be marked as inactive
    EXPECT_FALSE(square.Node(ear).is_active);
}

TEST(PolygonMeshClipEar, ClipRelinksPreviousAndNext) {
    PolygonMesh<double> square(kCcwSquare);
    const std::size_t ear = square.EarHeadIndex();
    const std::size_t prev = square.Node(ear).polygon_previous;
    const std::size_t next = square.Node(ear).polygon_next;
    square.ClipEarVertex(ear);
    // Previous and next vertices should be relinked
    EXPECT_EQ(square.Node(prev).polygon_next, next);
    EXPECT_EQ(square.Node(next).polygon_previous, prev);
}

TEST(PolygonMeshClipEar, ClipRemovesFromEarList) {
    PolygonMesh<double> square(kCcwSquare);
    const std::size_t ear = square.EarHeadIndex();
    square.ClipEarVertex(ear);
    // Clipped vertex should be removed from the ear list
    EXPECT_FALSE(square.Node(ear).is_ear);
}

TEST(PolygonMeshClipEar, ClipHeadEarAdvancesHead) {
    PolygonMesh<double> square(kCcwSquare);
    const std::size_t first = square.EarHeadIndex();
    square.ClipEarVertex(first);
    if (square.HasEarVertex()) {
        // Head should advance to the next ear
        EXPECT_NE(square.EarHeadIndex(), first);
    }
}

//------------------------------------------------------------------------------
// ClipEarVertex failure tests
//------------------------------------------------------------------------------

TEST(PolygonMeshClipEar, ThrowsOnNonEarVertex) {
    PolygonMesh<double> lshape(kCcwLShape);
    ASSERT_FALSE(lshape.Node(3).is_ear);
    // Should throw if we try to clip a non-ear vertex
    EXPECT_THROW(lshape.ClipEarVertex(3), std::invalid_argument);
}

TEST(PolygonMeshClipEar, ThrowsOnAlreadyClippedVertex) {
    PolygonMesh<double> square(kCcwSquare);
    const std::size_t ear = square.EarHeadIndex();
    square.ClipEarVertex(ear);
    // Should throw if we try to clip an already clipped vertex
    EXPECT_THROW(square.ClipEarVertex(ear), std::invalid_argument);
}

//------------------------------------------------------------------------------
// ActiveBoundaryVertices tests
//------------------------------------------------------------------------------

TEST(PolygonMeshBoundary, ReturnsAllVerticesAfterConstruction) {
    PolygonMesh<double> m(kCcwLShape);
    const auto b = m.ActiveBoundaryVertices();
    EXPECT_EQ(b.size(), 6U);
    std::set<std::size_t> s(b.begin(), b.end());
    EXPECT_EQ(s.size(), 6U);
}

TEST(PolygonMeshBoundary, OrderMatchesBoundaryRing) {
    PolygonMesh<double> m(kCcwSquare);
    const auto b = m.ActiveBoundaryVertices();
    for (std::size_t i = 0; i < b.size(); ++i) {
        EXPECT_EQ(m.Node(b[i]).polygon_next, b[(i + 1) % b.size()]);
    }
}

TEST(PolygonMeshBoundary, ShrinksByOneAfterClip) {
    PolygonMesh<double> m(kCcwSquare);
    m.ClipEarVertex(m.EarHeadIndex());
    EXPECT_EQ(m.ActiveBoundaryVertices().size(), 3U);
}

//------------------------------------------------------------------------------
// Invariants after each clip tests
//------------------------------------------------------------------------------

TEST(PolygonMeshInvariants, LShapeClipsPreserveInvariants) {
    PolygonMesh<double> m(kCcwLShape);
    const std::size_t n = 6;
    int steps = 0;
    while (m.ActiveVertexCount() > 3 && m.HasEarVertex()) {
        ExpectBoundaryRingValid(m, m.ActiveVertexCount());
        ExpectNoDeadLinksInActiveBoundary(m, n);
        ExpectEarImpliesConvex(m, n);
        ExpectActiveVertexCountMatchesNodes(m, n);
        m.ClipEarVertex(m.EarHeadIndex());
        ++steps;
        ASSERT_LT(steps, 20);
    }
    ExpectBoundaryRingValid(m, m.ActiveVertexCount());
    ExpectActiveVertexCountMatchesNodes(m, n);
}

//------------------------------------------------------------------------------
// Integration triangulation tests
//------------------------------------------------------------------------------

TEST(PolygonMeshIntegration, SquareProducesTwoTriangles) {
    std::vector<Point2D> vertices = kCcwSquare;
    PolygonMesh<double> m(vertices);
    std::vector<std::array<std::size_t, 3>> triangles;
    ASSERT_NO_THROW(triangles = TriangulateByEarClipping(m, vertices));
    EXPECT_EQ(triangles.size(), 2U);
    const double poly_area = PolygonArea(vertices);
    EXPECT_NEAR(SumTriangleAreas(triangles, vertices), poly_area, 1e-9);
}

TEST(PolygonMeshIntegration, LShapeProducesFourTriangles) {
    std::vector<Point2D> vertices = kCcwLShape;
    PolygonMesh<double> m(vertices);
    std::vector<std::array<std::size_t, 3>> triangles;
    ASSERT_NO_THROW(triangles = TriangulateByEarClipping(m, vertices));
    EXPECT_EQ(triangles.size(), 4U);
    const double poly_area = PolygonArea(vertices);
    EXPECT_NEAR(SumTriangleAreas(triangles, vertices), poly_area, 1e-9);
}

TEST(PolygonMeshIntegration, ArrowProducesThreeTriangles) {
    std::vector<Point2D> vertices = kCcwMShape;
    PolygonMesh<double> m(vertices);
    std::vector<std::array<std::size_t, 3>> triangles;
    ASSERT_NO_THROW(triangles = TriangulateByEarClipping(m, vertices));
    EXPECT_EQ(triangles.size(), 3U);
    const double poly_area = PolygonArea(vertices);
    EXPECT_NEAR(SumTriangleAreas(triangles, vertices), poly_area, 1e-9);
}

TEST(PolygonMeshIntegration, TriangleProducesOneTriangle) {
    PolygonMesh<double> m(kCcwTriangle);
    const auto b = m.ActiveBoundaryVertices();
    ASSERT_EQ(b.size(), 3U);
    const double a =
        TriangleArea(m.Node(b[0]).position, m.Node(b[1]).position, m.Node(b[2]).position);
    const double poly_area =
        PolygonArea(std::span<const Point2D>(kCcwTriangle.data(), kCcwTriangle.size()));
    EXPECT_NEAR(a, poly_area, 1e-9);
}

//------------------------------------------------------------------------------
// Stress/degenerate corner cases tests
//------------------------------------------------------------------------------

TEST(PolygonMeshStress, LargeConvexPolygon) {
    constexpr std::size_t kN = 100;
    const auto vertices = RegularConvexPolygon(kN, 10.);
    PolygonMesh<double> m(vertices);
    EXPECT_EQ(m.ActiveVertexCount(), kN);
    EXPECT_EQ(CountReflexVertices(m, kN), 0U);
    EXPECT_TRUE(m.HasEarVertex());
    std::vector<std::array<std::size_t, 3>> triangles;
    ASSERT_NO_THROW(triangles = TriangulateByEarClipping(m, vertices));
    EXPECT_EQ(triangles.size(), kN - 2U);
    const double poly_area = PolygonArea(vertices);
    EXPECT_NEAR(SumTriangleAreas(triangles, vertices), poly_area, 1e-6);
}

TEST(PolygonMeshStress, NearlyCollinearShapeVertices) {
    std::vector<Point2D> vertices = {{0., 0.}, {1., 0.}, {2., 1e-12}, {1., 1.}, {0., 1.}};
    PolygonMesh<double> m(vertices);
    EXPECT_NO_THROW(CountEars(m, vertices.size()));
    ExpectEarImpliesConvex(m, vertices.size());
}

TEST(PolygonMeshStress, ThinEarPolygon) {
    std::vector<Point2D> vertices = {{0., 0.}, {100., 0.}, {100., 1e-6}, {50., 5e-7}, {0., 1.}};
    PolygonMesh<double> m(vertices);
    std::vector<std::array<std::size_t, 3>> triangles;
    ASSERT_NO_THROW(triangles = TriangulateByEarClipping(m, vertices));
    const double poly_area = PolygonArea(vertices);
    EXPECT_NEAR(SumTriangleAreas(triangles, vertices), poly_area, 1e-3);
}

TEST(PolygonMeshStress, DuplicateAdjacentVertices) {
    std::vector<Point2D> vertices = {{0., 0.}, {0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}};
    PolygonMesh<double> m(vertices);
    EXPECT_EQ(m.ActiveVertexCount(), 5U);
    // Coincident vertex typically classified non-convex
    EXPECT_FALSE(m.Node(1).is_convex);
}

TEST(PolygonMeshStress, CollinearShapeInteriorVertexReflexOthersStillHaveEars) {
    PolygonMesh<double> m(kCollinearShape);
    // Index 1 is a 180 deg interior angle (kCollinearShape boundary), classified reflex.
    EXPECT_FALSE(m.Node(1).is_convex);
    EXPECT_FALSE(m.Node(1).is_ear);
    // Other corners remain convex and can still be ear tips.
    EXPECT_TRUE(m.HasEarVertex());
    ExpectReflexVertexIsNeverEar(m, kCollinearShape.size());
}

TEST(PolygonMeshAccessors, EarVertexNeighborsOrdered) {
    PolygonMesh<double> m(kCcwSquare);
    const std::size_t i = 2;
    const auto n = m.EarVertexNeighbors(i);
    EXPECT_EQ(n[0], m.Node(i).polygon_previous);
    EXPECT_EQ(n[1], i);
    EXPECT_EQ(n[2], m.Node(i).polygon_next);
}

}  // namespace geometry_kernel::test
