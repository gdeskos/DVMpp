#include "DVMBase.hpp"
#include <armadillo>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

DVMBase::DVMBase(XmlHandler &xml, const std::string &timestamp)
    : m_vortex(xml, timestamp), m_vortsheet(xml, timestamp),
      m_probe(xml, timestamp)
{
    m_step = 0;

    // Some helpers so we can do less typing
    auto getVal = [&xml](const char *l, const char *p) {
        return xml.getValueAttribute(l, p);
    };

    auto getStr = [&xml](const char *l, const char *p) {
        return xml.getStringAttribute(l, p);
    };

    m_nu = getVal("constants", "nu");

    m_dt = getVal("time", "dt");
    m_steps = getVal("time", "steps");

    m_Ux = getVal("flow", "ux");
    m_Uz = getVal("flow", "uz");
    m_Ur = std::sqrt(m_Ux * m_Ux + m_Uz * m_Uz);

    auto seed = xml.getIntAttribute("constants", "seed");
    m_rand.seed(seed);

    auto scheme = getStr("time", "scheme");
    if (scheme == "euler") {
        m_scheme = Scheme::Euler;
    } else if (scheme == "RK3") {
        m_scheme = Scheme::RK3;
    } // Invalid cases dealt with by the xml handler

    auto surfacecross = getStr("algorithms", "surface_crossing");
    if (surfacecross == "DELETE") {
        m_surfcross = SurfaceCross::Delete;
    } else if (surfacecross == "ABSORB") {
        m_surfcross = SurfaceCross::Absorb;
    } else if (surfacecross == "REFLECT") {
        m_surfcross = SurfaceCross::Reflect;
    } // Invalid cases dealt with by the xml handler
}

void DVMBase::solve()
{
    // Timeloop
    for (unsigned j = 1; j <= get_steps(); j++) {

        increment_step();
        compute_step();

        // Output
        write_outputs();

        // Screen output
        std::cout << "Simulation time          = " << get_time() << "\tStep "
                  << j << "/" << get_steps() << std::endl;
        std::cout << "Number of vortex blobs   = " << m_vortex.size()
                  << std::endl;
    }
}

void DVMBase::compute_step()
{
    //************** Advection substep*****************//
    convect();

    // remesh();
    // For remeshing both triangulated meshes (with RBF) and structured meshes
    // with Morgenthal M6 interpolators will be used.

    //************** Diffusion substep ****************//
    // The diffusion problem in an infinite domain
    m_vortex.diffusion_random_walk(m_rand, m_nu, m_dt);

    // A diffusion problem with only a flux of vorticity in the boundaries
    // dgamma/dn=a
    VortexBlobs NewVortices = m_vortsheet.release_nascent_vortices_rw(m_rand);
    m_vortex.append_vortices(NewVortices);

    // If a large time step is used some vortices may cross the boundary due to
    // random walk!
    // Care is taken of these vortices by 1) Deleting them, 2) Absorbing them,
    // 3) Reflecting them back to the flow
    m_vortsheet.reflect(m_vortex);

    //********************* Computing Loads/ Moving body *********************//
    // We compute the loads at the end of the time step
    m_vortsheet.compute_loads(m_Ur);

    //********************* Merge/Delete Vortices ****************************//
}

double DVMBase::get_time() const
{
    return m_time;
}

unsigned DVMBase::get_steps() const
{
    return m_steps;
}

void DVMBase::increment_step()
{
    m_step++;
    m_time = m_step * m_dt;
}

void DVMBase::write_outputs()
{
    m_vortex.write_step(m_time, m_step);

    m_vortsheet.write_step(m_time);

    m_probe.write_step(m_time);
}

void DVMBase::convect()
{
    switch (m_scheme) {
    case Scheme::Euler:
        // Find the free-space velocity
        m_vortex.biotsavart();
        // Find the boundary-imposed velocities-these two should be combined
        // together into one
        m_vortsheet.solvevortexsheet(m_vortex);
        m_vortsheet.vortexsheetbc(m_vortex);

        // Added them together
        m_vortex.m_x += (m_vortex.m_u + m_vortex.m_uvs + m_Ux) * m_dt;
        m_vortex.m_z += (m_vortex.m_w + m_vortex.m_wvs + m_Uz) * m_dt;
        break;
    case Scheme::RK3:
        // coefficients for Low Storage-Runge Kutta 3rd
        const double a[3] = {0, -17. / 32, -32. / 27};
        const double b[3] = {1. / 4, 8. / 9, 3. / 4};

        // execute 3 stages of Low Storage Runge-Kutta 3rd
        Vector q1 = arma::zeros(m_vortex.size());
        Vector q2 = arma::zeros(m_vortex.size());

        for (int i = 0; i < 3; ++i) {

            // Compute the right-hand side of the system
            // Find the free-space velocity
            m_vortex.biotsavart();
            // Find the boundary-imposed velocities-these two may be combined
            // together into one
            m_vortsheet.solvevortexsheet(m_vortex);
            m_vortsheet.vortexsheetbc(m_vortex);

            q1 = a[i] * q1 + (m_vortex.m_u + m_vortex.m_uvs + m_Ux) * m_dt;
            q2 = a[i] * q2 + (m_vortex.m_w + m_vortex.m_wvs + m_Uz) * m_dt;

            m_vortex.m_x += b[i] * q1;
            m_vortex.m_z += b[i] * q2;
        }
        break;
    }
}
