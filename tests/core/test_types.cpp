#include <complex>

#include <gtest/gtest.h>

#include "core/types.hpp"
#include "test_fixtures.hpp"

namespace geometry_kernel::test {

using namespace geometry_kernel::core;

//------------------------------------------------------------------------------
// Concept correctness (compile-time assertions)
//------------------------------------------------------------------------------
//
// These static_asserts compile (or fail at compile-time) — no runtime TEST body needed

// Admitted types
static_assert(ScalarType<float>);
static_assert(ScalarType<double>);
static_assert(ScalarType<long double>);

// Rejected: integer family
static_assert(!ScalarType<int>);
static_assert(!ScalarType<unsigned int>);
static_assert(!ScalarType<long>);
static_assert(!ScalarType<bool>);
static_assert(!ScalarType<char>);

// Note: std::floating_point (and therefore ScalarType) admits cv-qualified floating-point types
// because std::is_floating_point<const double> is true. This is standard-conforming behavior.
// The concept does not strip cv qualifiers before checking.
static_assert(ScalarType<const double>);
static_assert(ScalarType<volatile float>);

// Rejected: non-arithmetic types
static_assert(!ScalarType<std::complex<double>>);
static_assert(!ScalarType<void>);

struct UserStruct {
    double x{};
};
static_assert(!ScalarType<UserStruct>);

}  // namespace geometry_kernel::test
