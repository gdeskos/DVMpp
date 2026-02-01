#ifndef DVMBASE_HPP
#define DVMBASE_HPP

#include "BackgroundMesh.hpp"
#include "BaseTypes.hpp"
#include "Exceptions.hpp"
#include "Probe.hpp"
#include "Random.hpp"
#include "VortexBlobs.hpp"
#include "VortexSheet.hpp"
#include "XmlHandler.hpp"
#include "pugiconfig.hpp"
#include "pugixml.hpp"
#include <algorithm>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

/** \class DVMBase
 * \brief Contains all the elements of the DVM
 * \file DVMBase.hpp */

/// Time-stepping scheme enumeration
enum class Scheme {
    Euler, ///< Euler explicit
    RK3,   ///< Low-storage Runge-Kutta 3rd order
};

/// Surface crossing algorithm enumeration
enum class SurfaceCross {
    Delete,  ///< Deletes the vortices entering the body
    /// Absorbs the vortices crossing the boundary.
    /** Vortices detained for one time step and their strength is added to
       the next set of nascent vortices switching also the release to the
       Q-distribution */
    Absorb,
    Reflect  ///< Reflects the vortices entering the body back to their image position
};

class DVMBase
{
    private:
    VortexBlobs m_vortex;    ///< The point vortices
    VortexSheet m_vortsheet; ///< The vortex sheet
    Probe m_probe;           ///< Velocity probe
    Random m_rand;           ///< Random number generation
    XmlHandler m_xml;        ///< Xml handler
    BackgroundMesh m_mesh;   ///< Background mesh for remeshing

    double m_nu; ///< Kinematic viscosity

    double m_Ux; ///< Free stream velocity in x
    double m_Uz; ///< Free stream velocity in y
    double m_Ur; ///< Magnitude of the free stream velocity

    double m_dt;      ///< Timestep
    double m_time;    ///< Current simulation time
    unsigned m_step;  ///< Current timestep
    unsigned m_steps; ///< Total timesteps

    /// Computes a particular timestep
    void compute_step();

    /** @name Timeloop control */
    ///@{

    /// Increment the time step
    void increment_step();

    /// Get number of total timesteps
    [[nodiscard]] unsigned get_steps() const;

    /// Simulation time of the current time step
    [[nodiscard]] double get_time() const;
    ///@}

    Scheme m_scheme; ///< Time-stepping scheme

    SurfaceCross m_surfcross; ///< Surface crossing algorithm

    double m_merge_threshold;   ///< Threshold for vortex merging (0 = disabled)
    unsigned m_merge_frequency; ///< Merge every N steps (0 = disabled)

    /// Convect point vortices
    void convect();

    /// Perform remeshing if enabled and at correct frequency
    void remesh();

    /// Write the output file
    /** Happens at each time step */
    void write_outputs();

    /// Get the velocities at a given point
    void probe_velocities();

    public:
    DVMBase(XmlHandler &xml, const std::string &timestamp);

    /// Solve the particular problem - contains the timeloop
    void solve();
};

#endif // DVMBASE_HPP
