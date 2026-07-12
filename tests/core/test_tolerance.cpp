#include <gtest/gtest.h>

#include "core/tolerance.hpp"

namespace geometry_kernel::test {

using namespace geometry_kernel::core;

//------------------------------------------------------------------------------
// robust_sign tests
//------------------------------------------------------------------------------

TEST(RobustSign, ClearlyPositive) {
    EXPECT_EQ(RobustSign(1.), 1);
    EXPECT_EQ(RobustSign(1e-3), 1);
    EXPECT_EQ(RobustSign(1e6), 1);
}

TEST(RobustSign, ClearlyNegative) {
    EXPECT_EQ(RobustSign(-1.), -1);
    EXPECT_EQ(RobustSign(-1e-3), -1);
    EXPECT_EQ(RobustSign(-1e6), -1);
}

TEST(RobustSign, ExactlyZero) {
    EXPECT_EQ(RobustSign(0.), 0);
}

TEST(RobustSign, ExactlyAtDefaultToleranceBoundary) {
    // Strict inequalities -> exactly +/-kTolerance rounds to zero
    EXPECT_EQ(RobustSign(kTolerance<double>), 0);
    EXPECT_EQ(RobustSign(-kTolerance<double>), 0);
}

TEST(RobustSign, JustOutsideDefaultToleranceBoundary) {
    const double eps = std::numeric_limits<double>::epsilon();
    EXPECT_EQ(RobustSign(kTolerance<double> + eps), 1);
    EXPECT_EQ(RobustSign(-kTolerance<double> - eps), -1);
}

TEST(RobustSign, JustInsideDefaultToleranceBoundary) {
    const double eps = std::numeric_limits<double>::epsilon();
    EXPECT_EQ(RobustSign(kTolerance<double> - eps), 0);
    EXPECT_EQ(RobustSign(-kTolerance<double> + eps), 0);
}

TEST(RobustSign, CustomToleranceWiderBand) {
    EXPECT_EQ(RobustSign(0.5, 1.), 0);  // inside wider band
    EXPECT_EQ(RobustSign(1.5, 1.), 1);  // outside wider band
    EXPECT_EQ(RobustSign(-1.5, 1.), -1);
}

TEST(RobustSign, ZeroToleranceOnlyExactZeroIsZero) {
    EXPECT_EQ(RobustSign(0., 0.), 0);
    EXPECT_EQ(RobustSign(1e-15, 0.), 1);
    EXPECT_EQ(RobustSign(-1e-15, 0.), -1);
}

}  // namespace geometry_kernel::test
