# geometry-kernel

A header-only C++20 library of computational geometry primitives. Some of the features
include: orientation tests, intersections, containment, closest-point queries, and
polygon triangulation.

---

## Design principles

- **Header-only.** Every file can be copied individually into another project
  with minimal transitive includes — no build step required to consume a
  single primitive.
- **2D and 3D are separate, not defaulted.** `Point2`/`Vector2` and
  `Point3`/`Vector3` are distinct types in parallel files (`point2.hpp` /
  `point3.hpp`, and so on up the stack). A consumer who only needs planar
  geometry never pulls in 3D machinery and vice versa.
- **One tolerance policy.** All floating-point sign comparisons go through
  `RobustSign()` / `kTolerance<T>` in `core/tolerance.hpp` — no function
  anywhere else compares floats with `==` or a hand-rolled epsilon.
- **Concept-constrained templates.** Every primitive is templated on
  `ScalarType` (`float`, `double`, `long double`), constrained via a C++20
  concept rather than SFINAE.
- **Free functions, not methods.** `Dot(a, b)`, `Cross(a, b)`,
  `TriangleOrientation(a, b, c)` — consistent with how the geometry itself is
  expressed mathematically, and matches the `std::` style of operating on
  value types.
- **Layered folders, one-directional dependencies.** `core/` has no
  dependencies beyond itself. `queries/` builds on `core/`. `algorithms/`
  builds on both. Nothing lower in the stack ever depends on something above
  it.

---

## Build

**Prerequisites:** CMake 3.20 or higher, a C++20 compiler (GCC 11 or higher, Clang 13 or higher, MSVC 19.29 or higher).

```bash
# Configure and build (tests enabled by default, fetches GoogleTest)
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build -j

# Build without tests
cmake -S . -B build -DBUILD_TESTING=OFF
cmake --build build -j
```

`geometry_kernel` itself is a header-only `INTERFACE` target — there is
nothing to link against. Point your include path at `src/` and include only
the headers you need.

## Running the tests

```bash
ctest --test-dir build --output-on-failure
```

**233 tests across 46 suites**, currently all passing, zero warnings under
`-Wall -Wextra -Wpedantic` (`/W4` on MSVC).

---

## What's implemented

| Folder | Contents |
|---|---|
| `core/` | `types.hpp` (`ScalarType`, `Orientation`), `tolerance.hpp` (`RobustSign`), `point2.hpp`/`point3.hpp`, `vector2.hpp`/`vector3.hpp` (`Dot`, `Cross`, affine operators), `predicates2.hpp`/`predicates3.hpp` (`TriangleOrientation`, `IsLeftTurn`, `EdgeSign3`), `area.hpp` (triangle & shoelace polygon area) |
| `queries/` | `closest_point2.hpp`/`closest_point3.hpp` (point–line, point–segment, point–triangle), `containment2.hpp`/`containment3.hpp` (`PointInTriangle`, `PointInPolygon`, `InPlaneTriangleContainment`), `intersection.hpp` (`SegmentsIntersect`, including collinear/touching cases) |
| `algorithms/` | `polygon_mesh.hpp` + `ear_clipping.hpp` — O(n²) ear-clipping triangulation of simple polygons (see below) |

Every function above has direct unit test coverage, mirrored 1:1 under
`tests/`.

---

## Algorithms & data structures

### Ear-clipping triangulation

A vertex of a simple polygon is an *ear tip* if it is convex (interior angle
< 180°) and no reflex vertex lies strictly inside the triangle formed by it
and its two polygon neighbors. By the two-ears theorem (Meisters, 1975),
every simple polygon with more than three vertices has at least two ear
tips. The algorithm iteratively clips an ear tip — adding its triangle to
the result and removing it from the polygon — until three vertices remain,
which form the final triangle.

**Naive complexity: O(n³)** — n clips, each requiring O(n) to find the next
ear and O(n) to check reflex vertices against the candidate triangle.

### `PolygonMesh` — the O(n²) optimization

Clipping one ear only changes the classification of its two neighbors.
Rather than re-scanning every vertex after each clip, `PolygonMesh`
maintains three doubly-linked lists threaded through a single flat array of
`VertexNode` objects:

| List | Contents | Purpose |
|---|---|---|
| boundary | all active polygon vertices | the shrinking ring |
| reflex | reflex vertices only | scanned during ear validity tests |
| ear | current ear tips | O(1) access to the next ear to clip |

Each clip removes the ear tip from all lists it belongs to in O(1), then
reclassifies only its two neighbors — re-testing convex/reflex status and
ear validity against the current reflex list, O(r) where r is the number of
reflex vertices. Total cost: **O(n) clips × O(n) per clip = O(n²)**, with
**O(n) space** — one `VertexNode` per input vertex, the three lists are
links within the same array.

**Reference:** Eberly, D. "Triangulation by Ear Clipping." 2002.
https://www.geometrictools.com/Documentation/TriangulationByEarClipping.pdf

### Shoelace formula

`SignedPolygonArea` (`core/area.hpp`) computes signed polygon area and
winding order in O(n) time, O(1) space. Positive is CCW, negative is CW —
`TriangulatePolygonWithEarClipping` uses this to silently normalize CW
input before triangulating.

---

## Project layout

```
src/
  core/
    types.hpp              ScalarType concept, Orientation enum
    tolerance.hpp          kTolerance, RobustSign
    point2.hpp             Point2<T>, Point2D
    point3.hpp             Point3<T>, Point3D
    vector2.hpp            Vector2<T>, Vector2D, Dot, Cross, affine transformations
    vector3.hpp            3D counterparts of vector2.hpp
    predicates2.hpp        2D predicates: TriangleOrientation, IsLeftTurn
    predicates3.hpp        3D predicates: EdgeSign3
    area.hpp               TriangleArea, SignedPolygonArea, PolygonArea
  queries/
    closest_point2.hpp     ClosestPointOnSegment2<T>, ClosestPointOnLine2<T>, ClosestPointOnPlane2<T>
    closest_point3.hpp     3D counterparts of closest_point2.hpp
    containment2.hpp       2D containment predicates: PointInTriangle<T>, PointInPolygon<T>, InPlaneTriangleContainment<T>
    containment3.hpp       3D counterparts of containment2.hpp
    intersection.hpp       2D intersection predicates: SegmentsIntersect, PointOnSegment
  algorithms/
    polygon_mesh.hpp       2D polygon mesh data structure: PolygonMesh<T> — boundary/reflex/ear linked lists
    ear_clipping.hpp       2D ear-clipping triangulation: TriangulatePolygonWithEarClipping<T>, TriangulationResult<T>
tests/
  core/                    per-file unit tests + test_fixtures.hpp (shared polygon fixtures)
  queries/                 per-file unit tests
  algorithms/              test_polygon_mesh.cpp (structural invariants), test_ear_clipping.cpp (end-to-end)
```
