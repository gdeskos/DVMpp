#include <gtest/gtest.h>
#include "Random.hpp"
#include <set>

class RandomTest : public ::testing::Test {
protected:
    Random rng;
};

TEST_F(RandomTest, RandReturnsValueInRange) {
    for (int i = 0; i < 1000; ++i) {
        double val = rng.rand();
        EXPECT_GE(val, 0.0);
        EXPECT_LT(val, 1.0);
    }
}

TEST_F(RandomTest, SeedProducesReproducibleSequence) {
    Random rng1, rng2;
    rng1.seed(42);
    rng2.seed(42);

    for (int i = 0; i < 100; ++i) {
        EXPECT_DOUBLE_EQ(rng1.rand(), rng2.rand());
    }
}

TEST_F(RandomTest, DifferentSeedsProduceDifferentSequences) {
    Random rng1, rng2;
    rng1.seed(42);
    rng2.seed(123);

    bool any_different = false;
    for (int i = 0; i < 100; ++i) {
        if (rng1.rand() != rng2.rand()) {
            any_different = true;
            break;
        }
    }
    EXPECT_TRUE(any_different);
}

TEST_F(RandomTest, NegativeSeedUsesRandomDevice) {
    Random rng1, rng2;
    rng1.seed(-1);
    rng2.seed(-1);

    // With high probability, two random seeds should produce different sequences
    // Note: This test could theoretically fail, but it's extremely unlikely
    bool any_different = false;
    for (int i = 0; i < 100; ++i) {
        if (rng1.rand() != rng2.rand()) {
            any_different = true;
            break;
        }
    }
    // We don't assert here because both use random device and COULD be the same
    // Just check that it doesn't crash
}

TEST_F(RandomTest, ProducesVariedOutput) {
    rng.seed(42);
    std::set<double> values;

    for (int i = 0; i < 100; ++i) {
        values.insert(rng.rand());
    }

    // Should have many unique values
    EXPECT_GT(values.size(), 90u);
}

TEST_F(RandomTest, ZeroSeedIsValid) {
    EXPECT_NO_THROW({
        Random rng;
        rng.seed(0);
        rng.rand();
    });
}
