#include "Probe.hpp"

#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;

Probe::Probe()
{
}

Probe::Probe(const XmlHandler &xml, const std::string &stamp)
{
    // Initialise the data members
    m_x = xml.getList("probe", "x");
    m_z = xml.getList("probe", "z");

    m_u.copy_size(m_x);
    m_w.copy_size(m_x);
    m_uvs.copy_size(m_x);
    m_wvs.copy_size(m_x);

    // Initialise the output file
    auto outdir = xml.getStringAttribute("io", "output_dir", true);
    fs::path outpath(outdir);
    m_outfile = (outpath / (stamp + "_probe.dat")).string();

    std::ofstream of(m_outfile);
    if (!of) {
        throw dvm::FileIOException::cannotOpen(m_outfile);
    }
    of << xml.getValueAttribute("time", "dt") << " # Time step\n"
       << xml.getValueAttribute("time", "steps") << " # Steps\n"
       << "Time [s]\tu [m/s]\tw [m/s]\n";
}

unsigned Probe::size() const
{
    return static_cast<unsigned>(m_x.n_elem);
}

void Probe::write_step(double time)
{
    std::ofstream of(m_outfile, std::ios_base::app);
    if (!of) {
        throw dvm::FileIOException::cannotWrite(m_outfile);
    }

    for (unsigned i = 0; i < size(); ++i) {
        of << time << "\t" << m_u(i) << "\t" << m_w(i) << "\n";
    }
}
