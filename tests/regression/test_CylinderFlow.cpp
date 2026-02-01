#include <gtest/gtest.h>
#include "VortexBlobs.hpp"
#include "VelocityKernel.hpp"
#include "Random.hpp"
#include "BaseTypes.hpp"
#include <cmath>
#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;

class CylinderFlowRegressionTest : public ::testing::Test {
protected:
    void SetUp() override {
    }
};

// Note: Full cylinder flow regression tests would require:
// 1. Setting up VortexSheet from body coordinates
// 2. Running multiple timesteps
// 3. Comparing drag coefficients to reference values
//
// For now, we include simpler regression tests that verify
// the core algorithms produce consistent results.

TEST_F(CylinderFlowRegressionTest, DiffusionRandomWalkStatistics) {
    // Test that diffusion random walk produces expected statistical behavior
    // After many steps, particles should spread according to sqrt(4*nu*dt)

    const int N = 1000;
    VortexBlobs blobs(N);
    Random rng;
    rng.seed(42);  // Fixed seed for reproducibility

    double nu = 0.001;
    double dt = 0.01;

    // Initialize all particles at origin
    for (int i = 0; i < N; ++i) {
        blobs.m_ID[i] = i;
        blobs.m_x(i) = 0.0;
        blobs.m_z(i) = 0.0;
        blobs.m_circ(i) = 1.0;
        blobs.m_sigma(i) = 0.1;
    }

    // Apply diffusion
    blobs.diffusion_random_walk(rng, nu, dt);

    // Calculate mean position (should be near zero)
    double mean_x = arma::mean(blobs.m_x);
    double mean_z = arma::mean(blobs.m_z);

    EXPECT_NEAR(mean_x, 0.0, 0.1);  // Allow some statistical variation
    EXPECT_NEAR(mean_z, 0.0, 0.1);

    // Calculate standard deviation
    double std_x = arma::stddev(blobs.m_x);
    double std_z = arma::stddev(blobs.m_z);

    // Expected std should be proportional to sqrt(4*nu*dt)
    // Note: The actual distribution is not exactly Gaussian due to the
    // log transformation, but the spread should be of this order
    double expected_spread = std::sqrt(4.0 * nu * dt);

    // Check that spread is reasonable (within factor of 2)
    EXPECT_GT(std_x, expected_spread * 0.5);
    EXPECT_LT(std_x, expected_spread * 2.0);
    EXPECT_GT(std_z, expected_spread * 0.5);
    EXPECT_LT(std_z, expected_spread * 2.0);
}

TEST_F(CylinderFlowRegressionTest, CirculationConservation) {
    // Total circulation should be conserved during Biot-Savart computation

    const int N = 50;
    VortexBlobs blobs(N);
    Random rng;
    rng.seed(123);

    double total_initial_circ = 0.0;
    for (int i = 0; i < N; ++i) {
        blobs.m_ID[i] = i;
        blobs.m_x(i) = rng.rand() * 10.0 - 5.0;
        blobs.m_z(i) = rng.rand() * 10.0 - 5.0;
        blobs.m_circ(i) = rng.rand() * 2.0 - 1.0;
        blobs.m_sigma(i) = 0.1;
        total_initial_circ += blobs.m_circ(i);
    }

    // Compute Biot-Savart (doesn't change circulation)
    blobs.biotsavart();

    // Total circulation should be unchanged
    EXPECT_DOUBLE_EQ(blobs.totalcirc(), total_initial_circ);
}

TEST_F(CylinderFlowRegressionTest, VortexAppendPreservesData) {
    // Test that appending vortices preserves original data

    VortexBlobs blobs1(3);
    for (int i = 0; i < 3; ++i) {
        blobs1.m_ID[i] = i;
        blobs1.m_x(i) = i * 1.0;
        blobs1.m_z(i) = i * 2.0;
        blobs1.m_circ(i) = i * 0.1;
        blobs1.m_sigma(i) = 0.1;
    }

    VortexBlobs blobs2(2);
    for (int i = 0; i < 2; ++i) {
        blobs2.m_ID[i] = i;
        blobs2.m_x(i) = 10.0 + i;
        blobs2.m_z(i) = 20.0 + i;
        blobs2.m_circ(i) = 1.0 + i * 0.1;
        blobs2.m_sigma(i) = 0.2;
    }

    double total_circ_before = blobs1.totalcirc() + blobs2.totalcirc();

    blobs1.append_vortices(blobs2);

    EXPECT_EQ(blobs1.size(), 5u);
    EXPECT_DOUBLE_EQ(blobs1.totalcirc(), total_circ_before);

    // Check original data is preserved
    EXPECT_DOUBLE_EQ(blobs1.m_x(0), 0.0);
    EXPECT_DOUBLE_EQ(blobs1.m_x(1), 1.0);
    EXPECT_DOUBLE_EQ(blobs1.m_x(2), 2.0);

    // Check appended data
    EXPECT_DOUBLE_EQ(blobs1.m_x(3), 10.0);
    EXPECT_DOUBLE_EQ(blobs1.m_x(4), 11.0);
}

TEST_F(CylinderFlowRegressionTest, ReproducibilityWithSeed) {
    // Test that computations are reproducible with fixed seed

    auto run_simulation = [](int seed) {
        VortexBlobs blobs(10);
        Random rng;
        rng.seed(seed);

        for (int i = 0; i < 10; ++i) {
            blobs.m_ID[i] = i;
            blobs.m_x(i) = rng.rand();
            blobs.m_z(i) = rng.rand();
            blobs.m_circ(i) = rng.rand() - 0.5;
            blobs.m_sigma(i) = 0.1;
        }

        blobs.biotsavart();
        blobs.diffusion_random_walk(rng, 0.001, 0.01);

        return blobs.m_x(0);  // Return first x position as checksum
    };

    double result1 = run_simulation(42);
    double result2 = run_simulation(42);
    double result3 = run_simulation(123);

    // Same seed should give same result
    EXPECT_DOUBLE_EQ(result1, result2);

    // Different seed should give different result
    EXPECT_NE(result1, result3);
}
