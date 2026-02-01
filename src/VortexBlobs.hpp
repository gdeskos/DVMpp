#ifndef VORTEXBLOBS_HPP
#define VORTEXBLOBS_HPP

#include "BaseTypes.hpp"
#include "Exceptions.hpp"
#include "Random.hpp"
#include "XmlHandler.hpp"

#include <memory>
#include <vector>

// Forward declaration
namespace dvm {
class VelocityKernel;
}

/** \class VortexBlobs
 * \brief Describes the vortex blobs
 * \file VortexBlobs.hpp */

class VortexBlobs
{
    public:
    std::vector<unsigned> m_ID; ///< Vortex blob ID
    Vector m_x;                 ///< x-coordinate
    Vector m_z;                 ///< z-coordinate
    Vector m_circ;              ///< circulation
    Vector m_sigma;             ///< vortex cut-off
    Vector m_u;                 ///< local x-velocity
    Vector m_w;                 ///< local z-velocity
    Vector m_uvs;
    Vector m_wvs;

    private:
    /// Output filename for the blob data
    std::string m_blobsfile;

    /// Output filename for the number of blobs at each timestep
    std::string m_numfile;

    /// Velocity kernel for Biot-Savart computation
    std::unique_ptr<dvm::VelocityKernel> m_kernel;

    public:
    VortexBlobs();

    /// Create vortex blobs instance that will write to file
    VortexBlobs(const XmlHandler &xml, const std::string &stamp);

    /// Constructor
    /** \creates N number of vortices and sets them to zero*/
    explicit VortexBlobs(unsigned N);

    /// Destructor
    ~VortexBlobs();

    /// Move constructor
    VortexBlobs(VortexBlobs&& other) noexcept;

    /// Move assignment
    VortexBlobs& operator=(VortexBlobs&& other) noexcept;

    /// Copy constructor
    VortexBlobs(const VortexBlobs& other);

    /// Copy assignment
    VortexBlobs& operator=(const VortexBlobs& other);

    /// Appends vortexblobs
    void append_vortices(VortexBlobs &NewVortBlobs);

    /// Biot-Savart relationship
    void biotsavart();

    /// Diffusion through random walk
    void diffusion_random_walk(Random &_rand, double nu, double dt);

    /// Find the total circulation
    [[nodiscard]] double totalcirc() const;

    /// Current number of vortex blobs
    [[nodiscard]] unsigned size() const;

    /// Print the location of each vortex blob
    void print_location() const;

    /// Print the velocity of each vortex blob
    void print_velocity() const;

    /// Print the circulation of each vortex blob
    void print_circulation() const;

    /// Write the blob info and number of vortices to file
    /** \param time Simulation time [s] */
    void write_step(double time, unsigned step);

    /// Set velocity kernel type ("direct" or "fmm")
    void setKernelType(const std::string& type);
};

#endif // VORTEXBLOBS_HPP
