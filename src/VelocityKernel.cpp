#include "VelocityKernel.hpp"
#include "Exceptions.hpp"
#include <cmath>
#include <omp.h>

namespace dvm {

std::unique_ptr<VelocityKernel> VelocityKernel::create(const std::string& type)
{
    if (type == "direct") {
        return std::make_unique<DirectKernel>();
    }
#ifdef USE_FMM
    else if (type == "fmm") {
        return std::make_unique<FMMKernel>();
    }
#endif
    else {
        throw AlgorithmException("Unknown velocity kernel type: " + type +
                                ". Available: direct"
#ifdef USE_FMM
                                ", fmm"
#endif
                                );
    }
}

void DirectKernel::computeSelfInduced(
    const std::vector<VortexParticle>& particles,
    std::vector<VelocityResult>& velocities)
{
    const auto N = particles.size();
    velocities.resize(N);

    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        double u_i = 0.0;
        double w_i = 0.0;

        for (size_t j = 0; j < N; ++j) {
            if (i != j) {
                const double dx_ij = particles[i].x - particles[j].x;
                const double dz_ij = particles[i].z - particles[j].z;
                const double dr_ij2 = dx_ij * dx_ij + dz_ij * dz_ij;

                const double sigma_j_sq = particles[j].sigma * particles[j].sigma;
                const double threshold = 10.0 * sigma_j_sq;
                const double rsigmasqr = 1.0 / sigma_j_sq;

                double dK_ij;
                if (dr_ij2 < threshold) {
                    dK_ij = (1.0 - std::exp(-dr_ij2 * rsigmasqr)) / dr_ij2;
                } else {
                    dK_ij = 1.0 / dr_ij2;
                }

                u_i -= dK_ij * dz_ij * particles[j].circ;
                w_i += dK_ij * dx_ij * particles[j].circ;
            }
        }

        velocities[i].u = u_i * math::rpi2;
        velocities[i].w = w_i * math::rpi2;
    }
}

void DirectKernel::computeInducedAt(
    const std::vector<VortexParticle>& sources,
    const std::vector<double>& target_x,
    const std::vector<double>& target_z,
    std::vector<VelocityResult>& velocities)
{
    const auto Nt = target_x.size();
    const auto Ns = sources.size();
    velocities.resize(Nt);

    #pragma omp parallel for
    for (size_t i = 0; i < Nt; ++i) {
        double u_i = 0.0;
        double w_i = 0.0;

        for (size_t j = 0; j < Ns; ++j) {
            const double dx_ij = target_x[i] - sources[j].x;
            const double dz_ij = target_z[i] - sources[j].z;
            const double dr_ij2 = dx_ij * dx_ij + dz_ij * dz_ij;

            if (dr_ij2 > 1e-15) {  // Avoid self-interaction issues
                const double sigma_j_sq = sources[j].sigma * sources[j].sigma;
                const double threshold = 10.0 * sigma_j_sq;
                const double rsigmasqr = 1.0 / sigma_j_sq;

                double dK_ij;
                if (dr_ij2 < threshold) {
                    dK_ij = (1.0 - std::exp(-dr_ij2 * rsigmasqr)) / dr_ij2;
                } else {
                    dK_ij = 1.0 / dr_ij2;
                }

                u_i -= dK_ij * dz_ij * sources[j].circ;
                w_i += dK_ij * dx_ij * sources[j].circ;
            }
        }

        velocities[i].u = u_i * math::rpi2;
        velocities[i].w = w_i * math::rpi2;
    }
}

#ifdef USE_FMM

FMMKernel::FMMKernel(int expansion_order, double theta)
    : m_expansion_order(expansion_order), m_theta(theta)
{
}

void FMMKernel::computeSelfInduced(
    const std::vector<VortexParticle>& particles,
    std::vector<VelocityResult>& velocities)
{
    // FMM implementation using exafmm-t would go here
    // For now, fall back to direct computation with a warning
    // This will be implemented when exafmm-t integration is complete

    // TODO: Implement actual FMM computation using exafmm-t
    // For the vortex blob kernel, we need to use a regularized Biot-Savart
    // which requires some adaptation of the standard 2D FMM

    // Fall back to direct method for now
    DirectKernel direct;
    direct.computeSelfInduced(particles, velocities);
}

void FMMKernel::computeInducedAt(
    const std::vector<VortexParticle>& sources,
    const std::vector<double>& target_x,
    const std::vector<double>& target_z,
    std::vector<VelocityResult>& velocities)
{
    // FMM implementation for target evaluation
    // Fall back to direct method for now
    DirectKernel direct;
    direct.computeInducedAt(sources, target_x, target_z, velocities);
}

#endif // USE_FMM

} // namespace dvm
