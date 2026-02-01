#include <gtest/gtest.h>
#include "VortexBlobs.hpp"
#include "VelocityKernel.hpp"
#include "BaseTypes.hpp"
#include <cmath>
#include <vector>

class BiotSavartIntegrationTest : public ::testing::Test {
protected:
    void SetUp() override {
    }
};

// Test: Vortex pair translating perpendicular to connecting line
TEST_F(BiotSavartIntegrationTest, VortexPairTranslation) {
    // Two vortices at distance d apart with equal but opposite circulation
    // should translate perpendicular to the line joining them

    VortexBlobs blobs(2);
    double d = 2.0;
    double Gamma = 1.0;
    double sigma = 0.01;

    blobs.m_ID[0] = 0;
    blobs.m_ID[1] = 1;
    blobs.m_x(0) = -d/2;
    blobs.m_x(1) = d/2;
    blobs.m_z(0) = 0.0;
    blobs.m_z(1) = 0.0;
    blobs.m_circ(0) = Gamma;
    blobs.m_circ(1) = -Gamma;
    blobs.m_sigma(0) = sigma;
    blobs.m_sigma(1) = sigma;

    blobs.biotsavart();

    // Analytical solution: w = Gamma / (2*pi*d)
    double expected_w = Gamma / (2.0 * math::pi * d);

    // Both vortices should move in positive z
    EXPECT_GT(blobs.m_w(0), 0.0);
    EXPECT_GT(blobs.m_w(1), 0.0);

    // Velocity should match analytical solution
    EXPECT_NEAR(blobs.m_w(0), expected_w, 0.01 * expected_w);
    EXPECT_NEAR(blobs.m_w(1), expected_w, 0.01 * expected_w);

    // No x-velocity (symmetric configuration)
    EXPECT_NEAR(blobs.m_u(0), 0.0, 1e-10);
    EXPECT_NEAR(blobs.m_u(1), 0.0, 1e-10);
}

// Test: Co-rotating vortex pair
TEST_F(BiotSavartIntegrationTest, CoRotatingVortexPair) {
    // Two vortices with same sign circulation rotate around their centroid

    VortexBlobs blobs(2);
    double d = 2.0;
    double Gamma = 1.0;
    double sigma = 0.01;

    blobs.m_ID[0] = 0;
    blobs.m_ID[1] = 1;
    blobs.m_x(0) = -d/2;
    blobs.m_x(1) = d/2;
    blobs.m_z(0) = 0.0;
    blobs.m_z(1) = 0.0;
    blobs.m_circ(0) = Gamma;
    blobs.m_circ(1) = Gamma;  // Same sign
    blobs.m_sigma(0) = sigma;
    blobs.m_sigma(1) = sigma;

    blobs.biotsavart();

    // Analytical: Each vortex rotates about center with velocity Gamma/(2*pi*d)
    double expected_w = Gamma / (2.0 * math::pi * d);

    // They should move in opposite z-directions
    // Left vortex (at -d/2) moves down (-z), right vortex (at +d/2) moves up (+z)
    EXPECT_NEAR(blobs.m_w(0), -expected_w, 0.01 * expected_w);
    EXPECT_NEAR(blobs.m_w(1), expected_w, 0.01 * expected_w);

    // No x-velocity
    EXPECT_NEAR(blobs.m_u(0), 0.0, 1e-10);
    EXPECT_NEAR(blobs.m_u(1), 0.0, 1e-10);
}

// Test: Regularization prevents singularity
TEST_F(BiotSavartIntegrationTest, RegularizationPreventsInfinity) {
    // Two vortices very close together
    VortexBlobs blobs(2);
    double d = 0.001;  // Very small distance
    double Gamma = 1.0;
    double sigma = 0.1;  // Regularization radius much larger than d

    blobs.m_ID[0] = 0;
    blobs.m_ID[1] = 1;
    blobs.m_x(0) = -d/2;
    blobs.m_x(1) = d/2;
    blobs.m_z(0) = 0.0;
    blobs.m_z(1) = 0.0;
    blobs.m_circ(0) = Gamma;
    blobs.m_circ(1) = -Gamma;
    blobs.m_sigma(0) = sigma;
    blobs.m_sigma(1) = sigma;

    blobs.biotsavart();

    // Velocity should be finite (regularized), not infinity
    EXPECT_TRUE(std::isfinite(blobs.m_u(0)));
    EXPECT_TRUE(std::isfinite(blobs.m_w(0)));
    EXPECT_TRUE(std::isfinite(blobs.m_u(1)));
    EXPECT_TRUE(std::isfinite(blobs.m_w(1)));
}

// Test: Four vortex leapfrogging configuration
TEST_F(BiotSavartIntegrationTest, FourVortexLeapfrog) {
    // Classic leapfrogging configuration: two pairs of opposite vortices
    VortexBlobs blobs(4);
    double d = 1.0;
    double Gamma = 1.0;
    double sigma = 0.05;

    // First pair
    blobs.m_ID[0] = 0;
    blobs.m_x(0) = -d;
    blobs.m_z(0) = 0.5;
    blobs.m_circ(0) = Gamma;
    blobs.m_sigma(0) = sigma;

    blobs.m_ID[1] = 1;
    blobs.m_x(1) = -d;
    blobs.m_z(1) = -0.5;
    blobs.m_circ(1) = -Gamma;
    blobs.m_sigma(1) = sigma;

    // Second pair
    blobs.m_ID[2] = 2;
    blobs.m_x(2) = d;
    blobs.m_z(2) = 0.5;
    blobs.m_circ(2) = Gamma;
    blobs.m_sigma(2) = sigma;

    blobs.m_ID[3] = 3;
    blobs.m_x(3) = d;
    blobs.m_z(3) = -0.5;
    blobs.m_circ(3) = -Gamma;
    blobs.m_sigma(3) = sigma;

    blobs.biotsavart();

    // All velocities should be finite
    for (unsigned i = 0; i < 4; ++i) {
        EXPECT_TRUE(std::isfinite(blobs.m_u(i)));
        EXPECT_TRUE(std::isfinite(blobs.m_w(i)));
    }

    // Due to symmetry, same circulation vortices should have same x-velocity
    EXPECT_NEAR(blobs.m_u(0), blobs.m_u(2), 1e-10);
    EXPECT_NEAR(blobs.m_u(1), blobs.m_u(3), 1e-10);
}

// Test: Circular vortex ring approximation
TEST_F(BiotSavartIntegrationTest, CircularRingApproximation) {
    // N vortices arranged in a circle should produce roughly circular motion
    const int N = 20;
    VortexBlobs blobs(N);
    double R = 1.0;  // Ring radius
    double total_Gamma = 1.0;
    double sigma = 0.05;

    for (int i = 0; i < N; ++i) {
        double theta = 2.0 * math::pi * i / N;
        blobs.m_ID[i] = i;
        blobs.m_x(i) = R * std::cos(theta);
        blobs.m_z(i) = R * std::sin(theta);
        blobs.m_circ(i) = total_Gamma / N;  // Equal distribution
        blobs.m_sigma(i) = sigma;
    }

    blobs.biotsavart();

    // Total circulation should be conserved
    EXPECT_NEAR(blobs.totalcirc(), total_Gamma, 1e-10);

    // All velocities should be finite
    for (int i = 0; i < N; ++i) {
        EXPECT_TRUE(std::isfinite(blobs.m_u(i)));
        EXPECT_TRUE(std::isfinite(blobs.m_w(i)));
    }
}

// Test: Velocity kernel consistency
TEST_F(BiotSavartIntegrationTest, KernelConsistencyWithVortexBlobs) {
    // Test that VortexBlobs.biotsavart() and DirectKernel give same results

    VortexBlobs blobs(3);
    blobs.m_ID[0] = 0;
    blobs.m_ID[1] = 1;
    blobs.m_ID[2] = 2;
    blobs.m_x(0) = 0.0;
    blobs.m_x(1) = 1.0;
    blobs.m_x(2) = 0.5;
    blobs.m_z(0) = 0.0;
    blobs.m_z(1) = 0.0;
    blobs.m_z(2) = 1.0;
    blobs.m_circ(0) = 1.0;
    blobs.m_circ(1) = -0.5;
    blobs.m_circ(2) = 0.3;
    blobs.m_sigma(0) = 0.1;
    blobs.m_sigma(1) = 0.1;
    blobs.m_sigma(2) = 0.1;

    // Compute using VortexBlobs
    blobs.biotsavart();
    Vector u_blobs = blobs.m_u;
    Vector w_blobs = blobs.m_w;

    // Compute using DirectKernel directly
    auto kernel = dvm::VelocityKernel::create("direct");
    std::vector<dvm::VortexParticle> particles(3);
    for (int i = 0; i < 3; ++i) {
        particles[i] = {blobs.m_x(i), blobs.m_z(i), blobs.m_circ(i), blobs.m_sigma(i)};
    }

    std::vector<dvm::VelocityResult> velocities;
    kernel->computeSelfInduced(particles, velocities);

    // Results should match
    for (int i = 0; i < 3; ++i) {
        EXPECT_DOUBLE_EQ(u_blobs(i), velocities[i].u);
        EXPECT_DOUBLE_EQ(w_blobs(i), velocities[i].w);
    }
}
