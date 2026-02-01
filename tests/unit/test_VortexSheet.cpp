#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "VortexSheet.hpp"
#include "Exceptions.hpp"

class VortexSheetTest : public ::testing::Test {
protected:
    void SetUp() override {
    }
};

TEST_F(VortexSheetTest, DefaultConstructorWorks) {
    VortexSheet sheet;
    // Default constructor doesn't initialize, so just check it doesn't crash
    SUCCEED();
}

// Note: More comprehensive VortexSheet tests would require XML configuration
// and body coordinate files. These are better handled as integration tests.

TEST_F(VortexSheetTest, GetForcesReturnsZeroForDefaultSheet) {
    VortexSheet sheet;
    auto [fx, fz] = sheet.get_forces();
    // Default constructed sheet has zero forces
    EXPECT_DOUBLE_EQ(fx, 0.0);
    EXPECT_DOUBLE_EQ(fz, 0.0);
}
