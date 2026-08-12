#include <array>
#include <cmath>
#include <span>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "algorithms/ear_clipping.hpp"
#include "core/area.hpp"
#include "core/point2.hpp"
#include "core/test_fixtures.hpp"
#include "core/types.hpp"

namespace geometry_kernel::test {

using namespace geometry_kernel::core;
using namespace geometry_kernel::algorithms;

//------------------------------------------------------------------------------
// Validation helpers
//------------------------------------------------------------------------------

/**
 * @brief Validates the structural and geometric invariants of a triangulation result.
 *
 * Checks:
 * - `result.triangles.size() == n - 2`
 * - All triangle indices `v0`, `v1`, `v2` are in `[0, result.vertices.size())`
 * - Sum of per-triangle areas equals `result.total_area` within `kTolerance`
 * - `SignedPolygonArea(result.vertices) > 0` (output is CCW, provided input is CCW)
 *
 * @param result The triangulation result to validate.
 * @param num_vertices The original vertex count of the input polygon.
 */
template <ScalarType T>
void ExpectValidTriangulation(const TriangulationResult<T>& result, std::size_t num_vertices) {
    // Check that the number of triangles is n - 2
    ASSERT_EQ(result.triangles.size(), num_vertices - 2);

    // Check that all triangle indices are in range
    const std::size_t num_output_vertices = result.vertices.size();
    ASSERT_GT(num_output_vertices, 0u);
    T sum_area{};
    for (const auto& triangle : result.triangles) {
        ASSERT_LT(triangle.v0, num_output_vertices);
        ASSERT_LT(triangle.v1, num_output_vertices);
        ASSERT_LT(triangle.v2, num_output_vertices);
        sum_area += TriangleArea<T>(
            result.vertices[triangle.v0], result.vertices[triangle.v1], result.vertices[triangle.v2]
        );
    }
    // Check that the sum of the triangle areas equals the total area
    EXPECT_NEAR(sum_area, result.total_area, kTolerance<T>);

    // Check that the output vertices are CCW, provided input is CCW
    EXPECT_GT(SignedPolygonArea<T>(result.vertices), T{});
}

//------------------------------------------------------------------------------
// Triangulation tests
//------------------------------------------------------------------------------

TEST(TriangulationBasic, Triangle) {
    const auto result = TriangulatePolygonWithEarClipping<double>(kCcwTriangle);
    EXPECT_EQ(result.triangles.size(), 1u);
    EXPECT_NEAR(result.total_area, 0.5, kTolerance<double>);
    ExpectValidTriangulation<double>(result, kCcwTriangle.size());
}

TEST(TriangulationBasic, Square) {
    const auto result = TriangulatePolygonWithEarClipping<double>(kCcwSquare);
    EXPECT_EQ(result.triangles.size(), 2u);
    EXPECT_NEAR(result.total_area, 1., kTolerance<double>);
    ExpectValidTriangulation<double>(result, kCcwSquare.size());
}

TEST(TriangulationBasic, ConcavePolygon) {
    const auto result = TriangulatePolygonWithEarClipping<double>(kCcwLShape);
    EXPECT_EQ(result.triangles.size(), 4u);
    const double expected =
        SignedPolygonArea(std::vector<Point2D>(kCcwLShape.begin(), kCcwLShape.end()));
    EXPECT_NEAR(result.total_area, expected, kTolerance<double>);
    ExpectValidTriangulation<double>(result, kCcwLShape.size());
}

//------------------------------------------------------------------------------
// Triangulation invariants tests
//------------------------------------------------------------------------------

TEST(TriangulationInvariants, TriangleCountIsNMinus2) {
    for (std::span<const Point2D> span_data :
         {std::span<const Point2D>(kCcwSquare), std::span<const Point2D>(kCcwLShape)}) {
        const auto result = TriangulatePolygonWithEarClipping<double>(span_data);
        EXPECT_EQ(result.triangles.size(), span_data.size() - 2);
    }
}

TEST(TriangulationInvariants, AllIndicesInRange) {
    const auto result = TriangulatePolygonWithEarClipping<double>(kCcwLShape);
    const std::size_t nv = result.vertices.size();
    for (const Triangle& triangle : result.triangles) {
        EXPECT_LT(triangle.v0, nv);
        EXPECT_LT(triangle.v1, nv);
        EXPECT_LT(triangle.v2, nv);
    }
}

TEST(TriangulationInvariants, AreaSumMatchesTotalArea) {
    const auto result = TriangulatePolygonWithEarClipping<double>(kCcwSquare);
    double sum = 0.;
    for (const Triangle& triangle : result.triangles) {
        sum += TriangleArea(
            result.vertices[triangle.v0], result.vertices[triangle.v1], result.vertices[triangle.v2]
        );
    }
    EXPECT_NEAR(sum, result.total_area, kTolerance<double>);
}

//------------------------------------------------------------------------------
// Triangulation error path tests
//------------------------------------------------------------------------------

TEST(TriangulationErrors, EmptyInput) {
    const std::vector<Point2D> empty{};
    EXPECT_THROW(TriangulatePolygonWithEarClipping<double>(empty), std::invalid_argument);
}

TEST(TriangulationErrors, OneVertex) {
    const std::vector<Point2D> one{Point2D{0., 0.}};
    EXPECT_THROW(TriangulatePolygonWithEarClipping<double>(one), std::invalid_argument);
}

TEST(TriangulationErrors, TwoVertices) {
    const std::vector<Point2D> two{Point2D{0., 0.}, Point2D{1., 0.}};
    EXPECT_THROW(TriangulatePolygonWithEarClipping<double>(two), std::invalid_argument);
}

//------------------------------------------------------------------------------
// Triangulation integration tests
//------------------------------------------------------------------------------

// Based on the fixture's own "simple_concave_poly.csv"-derived geometry (see test_fixtures.hpp)
TEST(TriangulationIntegration, SimpleConcavePolygon) {
    const auto result = TriangulatePolygonWithEarClipping<double>(kSimpleConcavePolygon);
    EXPECT_NEAR(result.total_area, 0.65, 1e-12);
    ExpectValidTriangulation<double>(result, kSimpleConcavePolygon.size());
}

}  // namespace geometry_kernel::test
