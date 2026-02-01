#ifndef PROBE_HPP
#define PROBE_HPP

#include "BaseTypes.hpp"
#include "Exceptions.hpp"
#include "XmlHandler.hpp"

/** \class Probe
 * \brief Describes a velocity probe.
 * \file Probe.hpp */

class Probe
{
    public:
    Vector m_x;   ///< x-location of the probe points
    Vector m_z;   ///< z-location of the probe points
    Vector m_u;   ///< x-direction velocity of the points
    Vector m_w;   ///< z-direction velocity of the points
    Vector m_uvs; ///< x-direction velocity induced by vortex sheet
    Vector m_wvs; ///< z-direction velocity induced by vortex sheet

    private:
    std::string m_outfile; ///< The name of the output file

    public:
    /// Default constructor
    Probe();

    /// Create a probe
    Probe(const XmlHandler &xml, const std::string &stamp);

    /// Return the size of the probe
    [[nodiscard]] unsigned size() const;

    /// Write the step to file
    void write_step(double time);
};

#endif // PROBE_HPP
