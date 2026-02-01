#include "VortexBlobs.hpp"
#include "VelocityKernel.hpp"

#include <filesystem>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

//************************* Constructors *************************************//
VortexBlobs::VortexBlobs()
    : m_kernel(dvm::VelocityKernel::create("direct"))
{
}

VortexBlobs::VortexBlobs(const XmlHandler &xml, const std::string &stamp)
    : m_kernel(dvm::VelocityKernel::create("direct"))
{
    // Initialise the output files with the headers
    auto outdir = xml.getStringAttribute("io", "output_dir", true);
    auto dt = xml.getValueAttribute("time", "dt");
    auto steps = xml.getValueAttribute("time", "steps");

    fs::path outpath(outdir);
    m_blobsfile = (outpath / (stamp + "_vortex.dat")).string();
    auto bf = std::ofstream(m_blobsfile);
    if (!bf) {
        throw dvm::FileIOException::cannotOpen(m_blobsfile);
    }
    bf << size() << " # Number of nodes\n"
       << dt << " # Time step\n"
       << steps << " # Steps\n"
       << "Time [s]\tx-position [m]\tz-position [m]\tcirculation\n";

    m_numfile = (outpath / (stamp + "_vortex_num.dat")).string();
    auto nf = std::ofstream(m_numfile);
    if (!nf) {
        throw dvm::FileIOException::cannotOpen(m_numfile);
    }
    nf << "Time [s]\tNumber of vortices\n";
}

//****************************** Public Methods ******************************//
VortexBlobs::VortexBlobs(unsigned N)
    : m_kernel(dvm::VelocityKernel::create("direct"))
{
    m_ID.resize(N);
    m_x.set_size(N);
    m_z.set_size(N);
    m_circ.set_size(N);
    m_sigma.set_size(N);
    m_u.set_size(N);
    m_w.set_size(N);
    m_uvs.set_size(N);
    m_wvs.set_size(N);
}

VortexBlobs::~VortexBlobs() = default;

VortexBlobs::VortexBlobs(VortexBlobs&& other) noexcept = default;

VortexBlobs& VortexBlobs::operator=(VortexBlobs&& other) noexcept = default;

VortexBlobs::VortexBlobs(const VortexBlobs& other)
    : m_ID(other.m_ID)
    , m_x(other.m_x)
    , m_z(other.m_z)
    , m_circ(other.m_circ)
    , m_sigma(other.m_sigma)
    , m_u(other.m_u)
    , m_w(other.m_w)
    , m_uvs(other.m_uvs)
    , m_wvs(other.m_wvs)
    , m_blobsfile(other.m_blobsfile)
    , m_numfile(other.m_numfile)
    , m_kernel(dvm::VelocityKernel::create(other.m_kernel ? other.m_kernel->name() : "direct"))
{
}

VortexBlobs& VortexBlobs::operator=(const VortexBlobs& other)
{
    if (this != &other) {
        m_ID = other.m_ID;
        m_x = other.m_x;
        m_z = other.m_z;
        m_circ = other.m_circ;
        m_sigma = other.m_sigma;
        m_u = other.m_u;
        m_w = other.m_w;
        m_uvs = other.m_uvs;
        m_wvs = other.m_wvs;
        m_blobsfile = other.m_blobsfile;
        m_numfile = other.m_numfile;
        m_kernel = dvm::VelocityKernel::create(other.m_kernel ? other.m_kernel->name() : "direct");
    }
    return *this;
}

void VortexBlobs::setKernelType(const std::string& type)
{
    m_kernel = dvm::VelocityKernel::create(type);
}

void VortexBlobs::copyOutputConfig(const VortexBlobs& other)
{
    m_blobsfile = other.m_blobsfile;
    m_numfile = other.m_numfile;
}

unsigned VortexBlobs::merge_blobs(double threshold)
{
    if (size() < 2) {
        return 0;
    }

    // Track which blobs have been merged (marked for deletion)
    std::vector<bool> merged(size(), false);
    unsigned merge_count = 0;

    // Iterate over all pairs
    for (unsigned i = 0; i < size(); ++i) {
        if (merged[i]) continue;

        for (unsigned j = i + 1; j < size(); ++j) {
            if (merged[j]) continue;

            const double gamma_i = m_circ(i);
            const double gamma_j = m_circ(j);
            const double gamma_sum = gamma_i + gamma_j;

            // Skip if circulations would cancel (avoid division by zero)
            if (std::abs(gamma_sum) < 1e-15) continue;

            // Compute distance between vortices
            const double dx = m_x(i) - m_x(j);
            const double dz = m_z(i) - m_z(j);
            const double dist = std::sqrt(dx * dx + dz * dz);

            // Merge criterion: |Γi*Γj/(Γi+Γj)| * |xi - xj| < threshold
            const double criterion = std::abs(gamma_i * gamma_j / gamma_sum) * dist;

            if (criterion < threshold) {
                // Merge j into i, conserving moments
                // Zeroth moment: Γ_new = Γi + Γj
                // First moment: x_new = (Γi*xi + Γj*xj) / (Γi + Γj)

                const double new_x = (gamma_i * m_x(i) + gamma_j * m_x(j)) / gamma_sum;
                const double new_z = (gamma_i * m_z(i) + gamma_j * m_z(j)) / gamma_sum;
                const double new_circ = gamma_sum;

                // Weighted average for sigma
                const double new_sigma = (std::abs(gamma_i) * m_sigma(i) +
                                          std::abs(gamma_j) * m_sigma(j)) /
                                         (std::abs(gamma_i) + std::abs(gamma_j));

                // Update vortex i with merged values
                m_x(i) = new_x;
                m_z(i) = new_z;
                m_circ(i) = new_circ;
                m_sigma(i) = new_sigma;
                m_u(i) = 0.0;  // Velocities will be recomputed
                m_w(i) = 0.0;
                m_uvs(i) = 0.0;
                m_wvs(i) = 0.0;

                // Mark j for removal
                merged[j] = true;
                ++merge_count;
            }
        }
    }

    if (merge_count == 0) {
        return 0;
    }

    // Create compacted array without merged vortices
    std::vector<unsigned> keep_indices;
    keep_indices.reserve(size() - merge_count);

    for (unsigned i = 0; i < size(); ++i) {
        if (!merged[i]) {
            keep_indices.push_back(i);
        }
    }

    VortexBlobs new_vortex(keep_indices.size());

    for (unsigned i = 0; i < keep_indices.size(); ++i) {
        unsigned src = keep_indices[i];
        new_vortex.m_ID[i] = i;
        new_vortex.m_x(i) = m_x(src);
        new_vortex.m_z(i) = m_z(src);
        new_vortex.m_circ(i) = m_circ(src);
        new_vortex.m_sigma(i) = m_sigma(src);
        new_vortex.m_u(i) = m_u(src);
        new_vortex.m_w(i) = m_w(src);
        new_vortex.m_uvs(i) = m_uvs(src);
        new_vortex.m_wvs(i) = m_wvs(src);
    }

    // Preserve output config
    new_vortex.copyOutputConfig(*this);

    // Replace with compacted array
    *this = std::move(new_vortex);

    return merge_count;
}

void VortexBlobs::append_vortices(VortexBlobs &NewVortBlobs)
{
    // Get the number of Vortices
    auto Nold = size();
    auto Nnew = NewVortBlobs.size();
    for (unsigned i = 0; i < Nnew; i++) {
        m_ID.push_back(Nold + NewVortBlobs.m_ID[i]);
    }
    m_x.insert_rows(m_x.n_elem, NewVortBlobs.m_x);
    m_z.insert_rows(m_z.n_elem, NewVortBlobs.m_z);
    m_circ.insert_rows(m_circ.n_elem, NewVortBlobs.m_circ);
    m_sigma.insert_rows(m_sigma.n_elem, NewVortBlobs.m_sigma);
    m_u.insert_rows(m_u.n_elem, NewVortBlobs.m_u);
    m_w.insert_rows(m_w.n_elem, NewVortBlobs.m_w);
    m_uvs.insert_rows(m_uvs.n_elem, NewVortBlobs.m_uvs);
    m_wvs.insert_rows(m_wvs.n_elem, NewVortBlobs.m_wvs);
}

void VortexBlobs::biotsavart()
{
    const auto N = size();
    if (N == 0) return;

    // Convert to VelocityKernel format
    std::vector<dvm::VortexParticle> particles(N);
    for (unsigned i = 0; i < N; ++i) {
        particles[i] = {m_x(i), m_z(i), m_circ(i), m_sigma(i)};
    }

    // Compute velocities using kernel
    std::vector<dvm::VelocityResult> velocities;
    m_kernel->computeSelfInduced(particles, velocities);

    // Copy back results
    for (unsigned i = 0; i < N; ++i) {
        m_u(i) = velocities[i].u;
        m_w(i) = velocities[i].w;
    }
}

void VortexBlobs::diffusion_random_walk(Random &_rand, double nu, double dt)
{
    for (unsigned i = 0; i < size(); i++) {
        // Generate two random numbers in the range 0...1
        const double R1 = _rand.rand();
        const double R2 = _rand.rand();

        // Calculate r and theta for the random walk
        const double rrw = std::sqrt(4.0 * nu * dt * std::log(1.0 / R1));
        const double thetarw = 2.0 * math::pi * R2;

        m_x(i) += rrw * std::cos(thetarw);
        m_z(i) += rrw * std::sin(thetarw);
    }
}

double VortexBlobs::totalcirc() const
{
    return arma::sum(m_circ);
}

unsigned VortexBlobs::size() const
{
    const auto ID_sz = m_ID.size();
    if ((ID_sz != m_x.size()) || (ID_sz != m_z.size())
        || (ID_sz != m_circ.size())
        || (ID_sz != m_sigma.size())
        || (ID_sz != m_u.size())
        || (ID_sz != m_w.size())
        || (ID_sz != m_uvs.size())
        || (ID_sz != m_wvs.size())) {
        throw dvm::SizeMismatchException("VortexBlobs");
    }
    return static_cast<unsigned>(ID_sz);
}

void VortexBlobs::print_location() const
{
    for (unsigned i = 0; i < size(); i++) {
        std::cout << " The location of Vortex Blob " << m_ID[i]
                  << " is (x,z) = (" << m_x(i) << "," << m_z(i) << ")"
                  << std::endl;
    }
}

void VortexBlobs::print_velocity() const
{
    for (unsigned i = 0; i < size(); i++) {
        std::cout << "(u,w) = (" << m_u(i) << "," << m_w(i) << ")" << std::endl;
    }
}

void VortexBlobs::print_circulation() const
{
    for (unsigned i = 0; i < size(); i++) {
        std::cout << " circ = " << m_circ(i) << std::endl;
    }
}

void VortexBlobs::write_step(double time, unsigned step)
{
    if (m_blobsfile.empty() || m_numfile.empty()) {
        throw dvm::FileIOException("VortexBlobs::write_step: must set filenames before writing");
    }

    // Make sure we open up to end of the file
    auto bf = std::ofstream(m_blobsfile, std::ios_base::app);
    if (!bf) {
        throw dvm::FileIOException::cannotWrite(m_blobsfile);
    }

    for (unsigned i = 0; i < size(); ++i) {
        bf << time << "\t" << m_x(i) << "\t" << m_z(i) << "\t" << m_circ(i)
           << "\n";
    }

    auto nf = std::ofstream(m_numfile, std::ios_base::app);
    if (!nf) {
        throw dvm::FileIOException::cannotWrite(m_numfile);
    }
    nf << step << "\t" << size() << "\n";
}
