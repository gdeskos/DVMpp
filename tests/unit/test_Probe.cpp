#include <gtest/gtest.h>
#include "Probe.hpp"

class ProbeTest : public ::testing::Test {
protected:
    void SetUp() override {
    }
};

TEST_F(ProbeTest, DefaultConstructorCreatesEmptyProbe) {
    Probe probe;
    EXPECT_EQ(probe.size(), 0u);
}

TEST_F(ProbeTest, SizeReturnsCorrectCount) {
    Probe probe;
    // Default probe has size 0
    EXPECT_EQ(probe.size(), 0u);
}

// Note: More comprehensive Probe tests would require XML configuration.
// These are better handled as integration tests.
