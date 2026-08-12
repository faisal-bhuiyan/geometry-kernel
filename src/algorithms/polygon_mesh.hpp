#pragma once

#include <array>
#include <limits>
#include <stdexcept>
#include <vector>

#include "core/point2.hpp"
#include "core/predicates2.hpp"
#include "core/types.hpp"
#include "queries/containment2.hpp"

namespace geometry_kernel::algorithms {

using namespace geometry_kernel::core;
using geometry_kernel::queries::PointInTriangle;

/**
 * @brief Sentinel index meaning "no linked-list neighbor".
 *
 * Used in VertexNode link fields (polygon_previous, polygon_next, etc.) to indicate that a vertex
 * has no predecessor or successor in a given list -- either because it was never a member, or
 * because it has been unlinked. Chosen as the maximum value of std::size_t so it can never alias a
 * valid array index. In practice the polygon size is bounded by available memory long before this
 * sentinel value could be reached as a legitimate index.
 */
inline constexpr std::size_t kNone{std::numeric_limits<std::size_t>::max()};

/**
 * @brief Implements the polygon data structure to support the ear-clipping triangulation algorithm.
 *
 * A naive ear-clipping implementation is O(n^3): for each of the n ears clipped,
 * it scans all n vertices to find an ear, and tests all n reflex vertices for each
 * candidate. This class reduces that to O(n^2) by maintaining three doubly-linked
 * lists simultaneously over a single array of VertexNode objects:
 *
 *   - polygon: cyclic list of active boundary vertices (the shrinking polygon)
 *   - reflex:  linear list of reflex vertices (scanned during ear validity tests)
 *   - ear:     cyclic list of current ear tips (picked from directly each iteration)
 *
 * Each clip removes a vertex in O(1) and reclassifies only its two neighbors,
 * so the total cost per clip is O(n) for the reflex scan, giving O(n^2) overall
 * time complexity.
 *
 * References:
 *   - Eberly, David. "Triangulation by Ear Clipping." 2002.
 *     https://www.geometrictools.com/Documentation/TriangulationByEarClipping.pdf
 */
template <ScalarType T>
class PolygonMesh {
public:
    //-----------------------------------------------------------------------------------------------------------------
    // VertexNode
    //-----------------------------------------------------------------------------------------------------------------

    /**
     * @brief Storage record for a single polygon vertex and its linked-list membership.
     *
     * Each node stores a 2D position and participates in up to three doubly-linked lists. All lists
     * share the same underlying VertexNode array, with links represented as array indices rather
     * than pointers:
     *
     *  - **Polygon boundary** (cyclic): ordered traversal of the active polygon edge loop.
     *      Cyclic because the polygon has no natural start or end vertex.
     *
     *  - **Reflex list** (linear): vertices whose interior angle exceeds 180 deg, iterated
     *      during the point-in-triangle test to check ear validity.
     *      Linear because it is only scanned, never picked from by position.
     *
     *  - **Ear tip list** (cyclic): convex vertices whose triangle contains no reflex vertex.
     *      Cyclic so that the ear-clipping loop can pick and remove any ear in O(1)
     *      without a special case for the head.
     *
     * @note Index links use `kNone` (sentinel) to mean "not a member of this list".
     */
    struct VertexNode {
        Point2<T> position{};  ///< 2D coordinates of this vertex

        // Polygon boundary vertices -> cyclic doubly-linked list
        // Both links are always set to valid indices while the vertex is active.
        // Set to kNone when the vertex is clipped from the polygon.
        std::size_t polygon_previous{kNone};  ///< Previous vertex in the boundary vertices list
        std::size_t polygon_next{kNone};      ///< Next vertex in the boundary vertices list

        // Reflex vertices list -> linear doubly-linked list
        // Only populated when this vertex is reflex (interior angle > 180 deg).
        // Both links are kNone if this vertex is convex or not yet classified.
        std::size_t reflex_previous{kNone};  ///< Previous vertex in the reflex vertices list
        std::size_t reflex_next{kNone};      ///< Next vertex in the reflex vertices list

        // Ear tip vertices list -> cyclic doubly-linked list
        // Only populated when this vertex is a current ear tip.
        // Both links are kNone if this vertex is not an ear tip.
        std::size_t ear_previous{kNone};  ///< Previous vertex in the ear tips list
        std::size_t ear_next{kNone};      ///< Next vertex in the ear tips list

        bool is_convex{false};  ///< True when interior angle < 180 deg (left turn at this vertex)
        bool is_ear{false};     ///< True when this vertex is a current ear tip
        bool is_active{true};   ///< True if the vertex is active in the polygon boundary
                                ///< ring, false after the vertex has been clipped
    };

    /**
     * @brief Constructs the polygon mesh and initialises all three linked lists.
     *
     * Runs in two passes over the vertex array:
     *  1. Classifies every vertex as convex or reflex and seeds the reflex list.
     *  2. Tests each convex vertex for ear-tip status and seeds the ear list.
     *
     * @param vertices CCW-wound simple polygon + no duplicate closing vertex.
     *                 behavior is undefined for self-intersecting polygons.
     * @throws std::invalid_argument if fewer than 3 vertices are provided.
     */
    explicit PolygonMesh(std::vector<Point2<T>> vertices) {
        if (vertices.size() < 3U) {
            throw std::invalid_argument("PolygonMesh: need at least 3 vertices");
        }

        // Allocate nodes and wire up the cyclic polygon boundary
        const std::size_t num_vertices{vertices.size()};
        this->nodes_.resize(num_vertices);
        this->active_vertex_count_ = num_vertices;
        this->polygon_head_ = 0;

        for (std::size_t index{0}; index < num_vertices; ++index) {
            this->nodes_[index].position = vertices[index];
            // indexing is modulo num_vertices to wrap around the polygon boundary
            this->nodes_[index].polygon_previous = (index + num_vertices - 1) % num_vertices;
            this->nodes_[index].polygon_next = (index + 1) % num_vertices;
        }

        // Step 1: classify each vertex as convex or reflex and build the reflex list
        // Reflex vertices list must be complete before ear tip vertices list is seeded (Step 2)
        for (std::size_t index{0}; index < num_vertices; ++index) {
            this->nodes_[index].is_convex = this->IsConvexVertex(index);
            // Vertex is not convex -> insert into reflex list
            if (!this->nodes_[index].is_convex) {
                this->InsertReflexVertex(index);
            }
        }

        // Step 2: test each convex vertex for ear-tip status and seed the ear list
        for (std::size_t index{0}; index < num_vertices; ++index) {
            if (this->nodes_[index].is_convex && this->IsEarVertex(index)) {
                this->InsertEarVertex(index);
            }
        }
    }

    /// Returns the number of active vertices in the polygon
    [[nodiscard]] std::size_t ActiveVertexCount() const { return this->active_vertex_count_; }

    /// Returns True if there exists any ear tip in the polygon
    [[nodiscard]] bool HasEarVertex() const { return this->ear_head_ != kNone; }

    /// Returns the index of the first ear tip in the polygon
    [[nodiscard]] std::size_t EarHeadIndex() const { return this->ear_head_; }

    /// Returns the vertex node at index
    [[nodiscard]] const VertexNode& Node(std::size_t index) const { return this->nodes_[index]; }

    /// Returns the previous, current, and next indices of the ear tip on the polygon boundary
    [[nodiscard]] std::array<std::size_t, 3> EarVertexNeighbors(std::size_t index) const {
        return {
            this->nodes_[index].polygon_previous,  // previous vertex on polygon boundary
            index,                                 // current vertex
            this->nodes_[index].polygon_next       // next vertex on polygon boundary
        };
    }

    /**
     * @brief Clips the ear tip at @p index, removes it from the polygon boundary,
     * reduces ActiveVertexCount() by one, and updates ear/reflex membership
     * of its two former neighbors.
     * @pre Node(index).is_ear && Node(index).is_active
     * @throws std::invalid_argument if the precondition is violated.
     */
    void ClipEarVertex(std::size_t index) {
        if (!this->nodes_[index].is_active || !this->nodes_[index].is_ear) {
            throw std::invalid_argument("ClipEarVertex: index is not an active ear tip");
        }

        const std::size_t previous{this->nodes_[index].polygon_previous};
        const std::size_t next{this->nodes_[index].polygon_next};

        // Unlink vertex from ear and reflex vertices lists
        this->UnlinkEarVertex(index);
        this->UnlinkReflexVertex(index);

        // Update the polygon boundary vertices list
        this->nodes_[previous].polygon_next = next;
        this->nodes_[next].polygon_previous = previous;
        if (this->polygon_head_ == index) {
            this->polygon_head_ = next;
        }

        // Mark the vertex as clipped, nullify its boundary links, and decrement the active count
        this->nodes_[index].is_active = false;
        this->nodes_[index].polygon_previous = kNone;
        this->nodes_[index].polygon_next = kNone;
        --this->active_vertex_count_;

        // Reclassify the previous and next vertices
        this->ReclassifyVertex(previous);
        this->ReclassifyVertex(next);
    }

    /**
     * @brief Returns the ordered indices of all active vertices on the polygon boundary.
     *
     * Traverses the cyclic polygon boundary starting from the polygon head and collects
     * each active vertex index in order. Used to collect the final triangle when only
     * three vertices remain.
     *
     * @return Vector of active vertex indices in boundary-ring order, or empty if
     *         the polygon has no active vertices.
     */
    [[nodiscard]] std::vector<std::size_t> ActiveBoundaryVertices() const {
        std::vector<std::size_t> output;

        // Guard clause -> return empty if there are no active vertices
        if (this->polygon_head_ == kNone || this->active_vertex_count_ == 0) {
            return output;
        }

        // Traverse the cyclic polygon boundary, collecting each vertex index
        std::size_t current = this->polygon_head_;
        do {
            output.emplace_back(current);
            current = this->nodes_[current].polygon_next;
        } while (current != this->polygon_head_);

        return output;
    }

private:
    //---------------------------------------------------------------------------
    // Data members
    //---------------------------------------------------------------------------

    std::vector<VertexNode> nodes_;       ///< Array of vertex nodes in the polygon
                                          ///< (used for all three lists)
    std::size_t active_vertex_count_{0};  ///< Number of active vertices in the polygon
    std::size_t polygon_head_{kNone};     ///< Index of first vertex in polygon boundary (list 1)
    std::size_t reflex_head_{kNone};      ///< Index of first vertex in reflex list (list 2)
    std::size_t ear_head_{kNone};         ///< Index of first vertex in ear tip list (list 3)

    //---------------------------------------------------------------------------
    // Reflex vertices list helpers
    //---------------------------------------------------------------------------

    /**
     * @brief Inserts a vertex at the front of the reflex vertices list.
     *
     * The reflex list is a linear doubly-linked list. Prepending keeps insertion
     * O(1) and iteration order does not matter for the ear-tip test.
     */
    void InsertReflexVertex(std::size_t index) {
        // Point the new node forward to the current head
        this->nodes_[index].reflex_next = this->reflex_head_;

        // Back-link the old head to the new node if the list was non-empty
        if (this->reflex_head_ != kNone) {
            this->nodes_[this->reflex_head_].reflex_previous = index;
        }

        // The new node has no predecessor and becomes the new head
        this->nodes_[index].reflex_previous = kNone;
        this->reflex_head_ = index;
    }

    /**
     * @brief Removes a vertex from the reflex vertices list.
     *
     * @note Idempotent: safe to call on a vertex that is not currently a list member
     * (e.g. a convex vertex). This allows unconditional calls in ReclassifyVertex
     * and ClipEarVertex without a prior membership check at the call site.
     */
    void UnlinkReflexVertex(std::size_t index) {
        const std::size_t previous{this->nodes_[index].reflex_previous};
        const std::size_t next{this->nodes_[index].reflex_next};

        // Guard: not a member if no links and not the head -> return
        if ((previous == kNone) && (next == kNone) && (this->reflex_head_ != index)) {
            return;
        }

        // Patch the predecessor or advance the head
        if (previous != kNone) {
            this->nodes_[previous].reflex_next = next;
        } else {
            this->reflex_head_ = next;
        }

        // Patch the successor
        if (next != kNone) {
            this->nodes_[next].reflex_previous = previous;
        }

        // Clear the vertex's own links to signal it is no longer a list member
        this->nodes_[index].reflex_previous = kNone;
        this->nodes_[index].reflex_next = kNone;
    }

    //---------------------------------------------------------------------------
    // Ear tip vertices list helpers
    //---------------------------------------------------------------------------

    /**
     * @brief Inserts a vertex at the front of the ear tip vertices list.
     *
     * The ear list is a cyclic doubly-linked list. Prepending keeps insertion O(1).
     * A single-element list is self-referential (previous and next both point to itself).
     */
    void InsertEarVertex(std::size_t index) {
        this->nodes_[index].is_ear = true;

        // First ear tip -> create a single-element cyclic list pointing to itself
        if (this->ear_head_ == kNone) {
            this->ear_head_ = index;
            this->nodes_[index].ear_previous = index;
            this->nodes_[index].ear_next = index;
            return;
        }

        // Prepend to the front: wire the new node between the tail and the current head,
        // then advance the head to the new node
        const std::size_t tail = this->nodes_[this->ear_head_].ear_previous;
        this->nodes_[index].ear_next = this->ear_head_;
        this->nodes_[index].ear_previous = tail;
        this->nodes_[tail].ear_next = index;
        this->nodes_[this->ear_head_].ear_previous = index;
        this->ear_head_ = index;
    }

    /**
     * @brief Removes a vertex from the ear tip vertices list.
     *
     * @note Idempotent: safe to call on a vertex that is not currently a list member.
     * This allows unconditional calls in ReclassifyVertex and ClipEarVertex
     * without a prior membership check at the call site.
     */
    void UnlinkEarVertex(std::size_t index) {
        // Guard: vertex is not a list member if is_ear is false
        if (!this->nodes_[index].is_ear) {
            return;
        }

        const std::size_t previous{this->nodes_[index].ear_previous};
        const std::size_t next{this->nodes_[index].ear_next};

        // Patch the predecessor and successor, or empty the list if this is the last ear tip
        if (next == index) {
            // Cyclic list: last ear tip points to itself -> empty the list
            this->ear_head_ = kNone;
        } else {
            this->nodes_[previous].ear_next = next;
            this->nodes_[next].ear_previous = previous;
            // Advance the head if this vertex was at the front
            if (this->ear_head_ == index) {
                this->ear_head_ = next;
            }
        }

        // Clear the vertex's own links and flag to signal it is no longer a list member
        this->nodes_[index].ear_previous = kNone;
        this->nodes_[index].ear_next = kNone;
        this->nodes_[index].is_ear = false;
    }

    //---------------------------------------------------------------------------
    // Vertex classification helpers
    //---------------------------------------------------------------------------

    /**
     * @brief Returns true if the interior angle at @p index is less than 180 deg.
     *
     * Convexity is determined by checking that the turn from the previous vertex
     * through the current vertex to the next vertex is a strict left turn (CCW).
     */
    [[nodiscard]] bool IsConvexVertex(std::size_t index) const {
        const auto& previous = this->nodes_[this->nodes_[index].polygon_previous].position;
        const auto& current = this->nodes_[index].position;
        const auto& next = this->nodes_[this->nodes_[index].polygon_next].position;
        return IsLeftTurn(previous, current, next);
    }

    /**
     * @brief Returns true if the vertex at @p index is an ear tip.
     *
     * A convex vertex is an ear tip if no reflex vertex (other than the triangle's
     * own corners) lies strictly inside the triangle formed by its two neighbors.
     *
     * @note Assumes the caller has already verified the vertex is convex.
     */
    [[nodiscard]] bool IsEarVertex(std::size_t index) const {
        const auto& previous = this->nodes_[this->nodes_[index].polygon_previous].position;
        const auto& current = this->nodes_[index].position;
        const auto& next = this->nodes_[this->nodes_[index].polygon_next].position;

        // Reject if any reflex vertex other than the ear triangle's own corners lies inside
        //
        // TODO: Add a check based on AABB of the triangle to rule out reflex vertices that are
        // definitely outside the triangle -> should be more efficient for large polygons
        for (std::size_t reflex_index{this->reflex_head_}; reflex_index != kNone;
             reflex_index = this->nodes_[reflex_index].reflex_next) {
            if (reflex_index == this->nodes_[index].polygon_previous || reflex_index == index ||
                reflex_index == this->nodes_[index].polygon_next) {
                continue;
            }
            if (PointInTriangle(this->nodes_[reflex_index].position, previous, current, next)) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Recomputes the convex/reflex classification of a vertex and updates
     * its membership in the reflex and ear tip lists accordingly.
     *
     * Called on the two neighbors of a just-clipped ear tip, since removing a
     * vertex from the boundary can change the interior angle at adjacent vertices.
     * Both Unlink calls are unconditional and idempotent -- see UnlinkReflexVertex
     * and UnlinkEarVertex for the membership guards.
     */
    void ReclassifyVertex(std::size_t index) {
        // Guard: skip vertices that have already been clipped
        if (!this->nodes_[index].is_active) {
            return;
        }

        // Remove from whichever lists the vertex currently belongs to
        this->UnlinkEarVertex(index);
        this->UnlinkReflexVertex(index);

        // Re-test convexity and re-insert into appropriate lists
        this->nodes_[index].is_convex = this->IsConvexVertex(index);
        if (!this->nodes_[index].is_convex) {
            this->InsertReflexVertex(index);
            return;
        }
        if (this->IsEarVertex(index)) {
            this->InsertEarVertex(index);
        }
    }
};

}  // namespace geometry_kernel::algorithms
