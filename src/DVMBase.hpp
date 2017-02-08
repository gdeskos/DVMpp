#ifndef DVMBASE
#define DVMBASE

#include "BaseTypes.hpp"
#include "VortexBlobs.hpp"
#include "VortexSheet.hpp"
#include "Probe.hpp"
#include "XmlHandler.hpp"
#include "Random.hpp"
#include "pugiconfig.hpp"
#include "pugixml.hpp"
#include <algorithm>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <string>

/// DMV class contaning all the elements of the DVM

class DVMBase
{
	private:
	VortexBlobs m_vortex;    ///< The point vortices
	VortexSheet m_vortsheet; ///< The vortex sheet
	Probe m_probe;           ///< Velocity probe
	Random m_rand; ///< Random number generation
	XmlHandler m_xml; ///< Xml handler

	double m_nu; ///< Kinematic viscosity

	double m_Ux; ///< Free stream velocity in x
	double m_Uz; ///< Free stream velocity in ym_n, m_np;
	double m_Ur; ///< Magnitute of the free stream velocity in the magnitude of the flow

	double m_dt;            ///< Timestep
	double m_time;          ///< Current simulation time
	unsigned m_step;        ///< Current timestep
	unsigned m_steps;       ///< Total timesteps

	/// Computes a particular timestep
	void compute_step();

	/** @name Timeloop control */
	///@{

	/// Increment the time step
	void increment_step();

	/// Get number of total timesteps
	unsigned get_steps();

	/// Simulation time of the current time step
	double get_time();
	///@}

	/// Definition time-stepping scheme
	enum Scheme {
		Euler, ///< Euler explicit
		RK3,   ///< Low-storage Runge-Kutta 3nd order
		RK4    ///< Low-storage Runge-Kutta 4th order
	};
	Scheme m_scheme; ///< Time-stepping scheme
    
    /// Definition surface_crossing algorithm
	enum SurfaceCross {
		DELETE, ///< Deletes the vortices entering the body
		/// Absorbs the vortices the vortices crossing the boundary.
		/** Vorticies detained for one time step and their strngth is added to
	       the next set of nascent vortices switching also the release to the
	       Q-distribution*/
		ABSORB,
		REFLECT ///< Reflects the vortices entering the body back to their image
		        /// position
	};
	SurfaceCross m_surfcross; ///< Surface crossing algorithm 

	/// Convect point vortices
	void convect();

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

#endif
