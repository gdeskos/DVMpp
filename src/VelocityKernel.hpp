#ifndef VELOCITY_KERNEL_HPP
#define VELOCITY_KERNEL_HPP

/** \file VelocityKernel.hpp
 * \brief Interface for velocity computation kernels (Direct, FMM) */

#include "BaseTypes.hpp"
#include <memory>
#include <string>
#include <vector>

namespace dvm {

/// Result of velocity computation for a single particle
struct VelocityResult {
    double u = 0.0;  ///< x-velocity
    double w = 0.0;  ///< z-velocity
};

/// Particle data for velocity computation
struct VortexParticle {
    double x;      ///< x-coordinate
    double z;      ///< z-coordinate
    double circ;   ///< circulation
    double sigma;  ///< cutoff radius
};

/// Abstract base class for velocity computation kernels
class VelocityKernel {
public:
    virtual ~VelocityKernel() = default;

    /// Compute self-induced velocities for all particles
    /** This is the Biot-Savart velocity computation
     * \param particles Vector of vortex particles
     * \param velocities Output vector of computed velocities (must be pre-sized)
     */
    virtual void computeSelfInduced(
        const std::vector<VortexParticle>& particles,
        std::vector<VelocityResult>& velocities) = 0;

    /// Compute velocities induced by sources on targets
    /** \param sources Vortex particles that induce velocity
     * \param target_x x-coordinates of target points
     * \param target_z z-coordinates of target points
     * \param velocities Output vector of computed velocities
     */
    virtual void computeInducedAt(
        const std::vector<VortexParticle>& sources,
        const std::vector<double>& target_x,
        const std::vector<double>& target_z,
        std::vector<VelocityResult>& velocities) = 0;

    /// Factory method to create a kernel by name
    /** \param type "direct" for O(N^2) method, "fmm" for Fast Multipole Method
     * \return Unique pointer to the created kernel
     */
    [[nodiscard]] static std::unique_ptr<VelocityKernel> create(const std::string& type);

    /// Get kernel name
    [[nodiscard]] virtual std::string name() const = 0;
};

/// Direct O(N^2) Biot-Savart computation
class DirectKernel : public VelocityKernel {
public:
    DirectKernel() = default;

    void computeSelfInduced(
        const std::vector<VortexParticle>& particles,
        std::vector<VelocityResult>& velocities) override;

    void computeInducedAt(
        const std::vector<VortexParticle>& sources,
        const std::vector<double>& target_x,
        const std::vector<double>& target_z,
        std::vector<VelocityResult>& velocities) override;

    [[nodiscard]] std::string name() const override { return "direct"; }
};

#ifdef USE_FMM
/// Fast Multipole Method O(N) computation using exafmm-t
class FMMKernel : public VelocityKernel {
public:
    /// Constructor with FMM parameters
    /** \param expansion_order Multipole expansion order (default 10)
     * \param theta Opening angle criterion (default 0.5)
     */
    explicit FMMKernel(int expansion_order = 10, double theta = 0.5);

    void computeSelfInduced(
        const std::vector<VortexParticle>& particles,
        std::vector<VelocityResult>& velocities) override;

    void computeInducedAt(
        const std::vector<VortexParticle>& sources,
        const std::vector<double>& target_x,
        const std::vector<double>& target_z,
        std::vector<VelocityResult>& velocities) override;

    [[nodiscard]] std::string name() const override { return "fmm"; }

private:
    int m_expansion_order;
    double m_theta;
};
#endif // USE_FMM

} // namespace dvm

#endif // VELOCITY_KERNEL_HPP
