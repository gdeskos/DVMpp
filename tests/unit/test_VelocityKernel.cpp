#include <gtest/gtest.h>
#include "VelocityKernel.hpp"
#include "Exceptions.hpp"
#include <cmath>

class VelocityKernelTest : public ::testing::Test {
protected:
    void SetUp() override {
    }
};

TEST_F(VelocityKernelTest, CreateDirectKernel) {
    auto kernel = dvm::VelocityKernel::create("direct");
    ASSERT_NE(kernel, nullptr);
    EXPECT_EQ(kernel->name(), "direct");
}

TEST_F(VelocityKernelTest, CreateInvalidKernelThrows) {
    EXPECT_THROW({
        dvm::VelocityKernel::create("invalid");
    }, dvm::AlgorithmException);
}

TEST_F(VelocityKernelTest, DirectKernelComputesSelfInducedVelocity) {
    auto kernel = dvm::VelocityKernel::create("direct");

    // Create two vortices
    std::vector<dvm::VortexParticle> particles = {
        {0.0, 0.0, 1.0, 0.1},   // x, z, circ, sigma
        {1.0, 0.0, -1.0, 0.1}
    };

    std::vector<dvm::VelocityResult> velocities;
    kernel->computeSelfInduced(particles, velocities);

    ASSERT_EQ(velocities.size(), 2u);

    // Both vortices should have non-zero z-velocity (moving perpendicular to line connecting them)
    // They should move in the same direction since they have opposite circulations
    EXPECT_NE(velocities[0].w, 0.0);
    EXPECT_NE(velocities[1].w, 0.0);

    // Sign should be the same (both moving in same z-direction)
    EXPECT_GT(velocities[0].w * velocities[1].w, 0.0);
}

TEST_F(VelocityKernelTest, DirectKernelEmptyParticlesReturnsEmpty) {
    auto kernel = dvm::VelocityKernel::create("direct");

    std::vector<dvm::VortexParticle> particles;
    std::vector<dvm::VelocityResult> velocities;

    kernel->computeSelfInduced(particles, velocities);

    EXPECT_EQ(velocities.size(), 0u);
}

TEST_F(VelocityKernelTest, DirectKernelSingleParticleHasZeroSelfVelocity) {
    auto kernel = dvm::VelocityKernel::create("direct");

    std::vector<dvm::VortexParticle> particles = {
        {0.0, 0.0, 1.0, 0.1}
    };

    std::vector<dvm::VelocityResult> velocities;
    kernel->computeSelfInduced(particles, velocities);

    ASSERT_EQ(velocities.size(), 1u);
    EXPECT_DOUBLE_EQ(velocities[0].u, 0.0);
    EXPECT_DOUBLE_EQ(velocities[0].w, 0.0);
}

TEST_F(VelocityKernelTest, DirectKernelComputesInducedAtTargets) {
    auto kernel = dvm::VelocityKernel::create("direct");

    // Single vortex at origin
    std::vector<dvm::VortexParticle> sources = {
        {0.0, 0.0, 1.0, 0.1}
    };

    // Target point at (1, 0)
    std::vector<double> target_x = {1.0};
    std::vector<double> target_z = {0.0};

    std::vector<dvm::VelocityResult> velocities;
    kernel->computeInducedAt(sources, target_x, target_z, velocities);

    ASSERT_EQ(velocities.size(), 1u);

    // Vortex at origin with positive circulation should induce positive w at (1, 0)
    // u should be zero (target is on the x-axis)
    EXPECT_NEAR(velocities[0].u, 0.0, 1e-10);

    // w should be positive and approximately Gamma / (2*pi*r)
    double expected_w = 1.0 / (2.0 * math::pi * 1.0);
    EXPECT_NEAR(velocities[0].w, expected_w, 0.01);
}

TEST_F(VelocityKernelTest, BiotSavartVortexPairAnalytical) {
    // Two vortices of equal and opposite strength at distance d
    // should move perpendicular to the line joining them
    auto kernel = dvm::VelocityKernel::create("direct");

    double d = 2.0;
    double Gamma = 1.0;
    double sigma = 0.01; // Small sigma for nearly singular kernel

    std::vector<dvm::VortexParticle> particles = {
        {-d/2, 0.0, Gamma, sigma},
        {d/2, 0.0, -Gamma, sigma}
    };

    std::vector<dvm::VelocityResult> velocities;
    kernel->computeSelfInduced(particles, velocities);

    // Analytical: w = Gamma / (2*pi*d) for both particles
    double expected_w = Gamma / (2.0 * math::pi * d);

    EXPECT_NEAR(velocities[0].w, expected_w, 0.01 * expected_w);
    EXPECT_NEAR(velocities[1].w, expected_w, 0.01 * expected_w);

    // u should be nearly zero
    EXPECT_NEAR(velocities[0].u, 0.0, 1e-6);
    EXPECT_NEAR(velocities[1].u, 0.0, 1e-6);
}

#ifdef USE_FMM
TEST_F(VelocityKernelTest, CreateFMMKernel) {
    auto kernel = dvm::VelocityKernel::create("fmm");
    ASSERT_NE(kernel, nullptr);
    EXPECT_EQ(kernel->name(), "fmm");
}

TEST_F(VelocityKernelTest, FMMKernelMatchesDirectForSmallN) {
    auto direct_kernel = dvm::VelocityKernel::create("direct");
    auto fmm_kernel = dvm::VelocityKernel::create("fmm");

    // Create N vortices in a random configuration
    const int N = 100;
    std::vector<dvm::VortexParticle> particles(N);

    std::mt19937 rng(42);
    std::uniform_real_distribution<> pos_dist(-10.0, 10.0);
    std::uniform_real_distribution<> circ_dist(-1.0, 1.0);

    for (int i = 0; i < N; ++i) {
        particles[i] = {
            pos_dist(rng),
            pos_dist(rng),
            circ_dist(rng),
            0.1
        };
    }

    std::vector<dvm::VelocityResult> direct_vel, fmm_vel;
    direct_kernel->computeSelfInduced(particles, direct_vel);
    fmm_kernel->computeSelfInduced(particles, fmm_vel);

    // FMM should match direct within tolerance
    for (int i = 0; i < N; ++i) {
        EXPECT_NEAR(fmm_vel[i].u, direct_vel[i].u, 1e-6);
        EXPECT_NEAR(fmm_vel[i].w, direct_vel[i].w, 1e-6);
    }
}
#endif
