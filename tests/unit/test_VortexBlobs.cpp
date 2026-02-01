#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "VortexBlobs.hpp"
#include "Exceptions.hpp"
#include <cmath>

class VortexBlobsTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Common setup for all tests
    }
};

// Constructor tests
TEST_F(VortexBlobsTest, DefaultConstructorCreatesEmptyBlobs) {
    VortexBlobs blobs;
    EXPECT_EQ(blobs.size(), 0u);
}

TEST_F(VortexBlobsTest, SizedConstructorCreatesCorrectSize) {
    VortexBlobs blobs(10);
    EXPECT_EQ(blobs.size(), 10u);
}

TEST_F(VortexBlobsTest, SizedConstructorInitializesVectors) {
    VortexBlobs blobs(5);
    EXPECT_EQ(blobs.m_x.size(), 5u);
    EXPECT_EQ(blobs.m_z.size(), 5u);
    EXPECT_EQ(blobs.m_circ.size(), 5u);
    EXPECT_EQ(blobs.m_sigma.size(), 5u);
    EXPECT_EQ(blobs.m_u.size(), 5u);
    EXPECT_EQ(blobs.m_w.size(), 5u);
    EXPECT_EQ(blobs.m_ID.size(), 5u);
}

// Append tests
TEST_F(VortexBlobsTest, AppendVorticesIncreasesSize) {
    VortexBlobs blobs1(3);
    VortexBlobs blobs2(2);

    // Initialize blobs2
    for (unsigned i = 0; i < 2; ++i) {
        blobs2.m_ID[i] = i;
        blobs2.m_x(i) = i * 1.0;
        blobs2.m_z(i) = i * 2.0;
        blobs2.m_circ(i) = 0.1;
        blobs2.m_sigma(i) = 0.05;
    }

    blobs1.append_vortices(blobs2);
    EXPECT_EQ(blobs1.size(), 5u);
}

// Biot-Savart tests
TEST_F(VortexBlobsTest, BiotSavartTwoVorticesOppositeCirculation) {
    VortexBlobs blobs(2);

    // Two vortices at distance d apart, opposite circulation
    double d = 1.0;
    double circ = 1.0;
    double sigma = 0.1;

    blobs.m_ID[0] = 0;
    blobs.m_ID[1] = 1;
    blobs.m_x(0) = -d/2;
    blobs.m_x(1) = d/2;
    blobs.m_z(0) = 0.0;
    blobs.m_z(1) = 0.0;
    blobs.m_circ(0) = circ;
    blobs.m_circ(1) = -circ;
    blobs.m_sigma(0) = sigma;
    blobs.m_sigma(1) = sigma;

    blobs.biotsavart();

    // Both vortices should move in the same z-direction
    // The velocity magnitude should be approximately Gamma/(2*pi*d)
    double expected_w = circ / (2.0 * math::pi * d);

    // Allow for some numerical tolerance due to regularization
    EXPECT_NEAR(blobs.m_w(0), expected_w, 0.1 * expected_w);
    EXPECT_NEAR(blobs.m_w(1), expected_w, 0.1 * expected_w);
}

TEST_F(VortexBlobsTest, BiotSavartEmptyBlobsDoesNotCrash) {
    VortexBlobs blobs;
    EXPECT_NO_THROW(blobs.biotsavart());
}

// Total circulation tests
TEST_F(VortexBlobsTest, TotalCirculationSumsCorrectly) {
    VortexBlobs blobs(3);
    blobs.m_circ(0) = 1.0;
    blobs.m_circ(1) = 2.0;
    blobs.m_circ(2) = -0.5;

    EXPECT_DOUBLE_EQ(blobs.totalcirc(), 2.5);
}

TEST_F(VortexBlobsTest, TotalCirculationEmptyBlobsIsZero) {
    VortexBlobs blobs;
    EXPECT_DOUBLE_EQ(blobs.totalcirc(), 0.0);
}

// Copy and move tests
TEST_F(VortexBlobsTest, CopyConstructorCreatesIndependentCopy) {
    VortexBlobs original(3);
    original.m_x(0) = 1.0;
    original.m_x(1) = 2.0;
    original.m_x(2) = 3.0;

    VortexBlobs copy(original);

    EXPECT_EQ(copy.size(), original.size());
    EXPECT_DOUBLE_EQ(copy.m_x(0), 1.0);

    // Modify original, copy should be unchanged
    original.m_x(0) = 100.0;
    EXPECT_DOUBLE_EQ(copy.m_x(0), 1.0);
}

TEST_F(VortexBlobsTest, MoveConstructorTransfersOwnership) {
    VortexBlobs original(3);
    original.m_x(0) = 1.0;

    VortexBlobs moved(std::move(original));

    EXPECT_EQ(moved.size(), 3u);
    EXPECT_DOUBLE_EQ(moved.m_x(0), 1.0);
}

// Kernel type tests
TEST_F(VortexBlobsTest, SetKernelTypeDirectWorks) {
    VortexBlobs blobs;
    EXPECT_NO_THROW(blobs.setKernelType("direct"));
}

TEST_F(VortexBlobsTest, SetKernelTypeInvalidThrows) {
    VortexBlobs blobs;
    EXPECT_THROW(blobs.setKernelType("invalid"), dvm::AlgorithmException);
}
