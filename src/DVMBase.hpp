#ifndef DVMBASE
#define DVMBASE

#include "BaseTypes.hpp"
#include "Body.hpp"
#include "VortexBlobs.hpp"
#include "VortexSheet.hpp"
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
	Body m_body;             ///< The solid body
	VortexBlobs m_probe;     ///< The probe???

	Matrix m_infM; ///< Influence matrix

	Random m_rand; ///< Random number generation

	/** @name Output files */
	///@{
	std::ofstream dev_dvm;   ///< Vortex position
	std::ofstream dev_Num;   ///< Number of vortices
	std::ofstream dev_gamma; ///< Vortex sheet strength
	std::ofstream dev_loads; ///< Forces
	std::ofstream dev_probe; ///< Velocity
	///@}

	double m_dt; ///< Timestep
	double m_nu; ///< Kinematic viscosity
	double m_Ux; ///< Free stream velocity in x
	double m_Uz; ///< Free stream velocity in ym_n, m_np;
	double m_n; ///< Number of body points

	double m_kernel_threshold; ///< Kernel threshold
	double m_sigma_cutoff;     ///< Cut-off distance

	double m_time;       ///< Current simulation time
	unsigned m_step;     ///< Current timestep
	unsigned m_steps;    ///< Total timesteps
	double m_maxGamma;   ///< Maximum value of the vorticity


	/** @name Reading and writing */

	///@{
	XmlHandler m_xml;          ///< xml handler
	std::string m_in_dir;      ///< input directory
	std::string m_out_dir;     ///< output directory
	std::string m_domain_file; ///< body geometry file
	///@}

	std::string m_timestamp;

	// some constants
	double m_pi;   ///< pi
	double m_rpi2; ///< 1 / 2pi

	// Loads related variables
	double m_rho; ///< Density of the fluid

	private:
	public:
	DVMBase(XmlHandler &xml);

	/// Definition time-stepping scheme
	enum Scheme {
		Euler, ///< Euler explicit
		RK2,   ///< Runge-Kutta 2nd order
		RK4    ///< Runge-Kutta 4th order
	};
	Scheme m_scheme; ///< Time-stepping scheme


	/** @name General Initialisation */
	///@{

	/// Initialise the object
	void init(XmlHandler &xml, std::string timestamp);

	/// Read the input geometry file
	void read_input_coord();

	/// Initialise the output files
	void init_outputs();
	///@}

	/// Compute the influence matrix
	/** Computed using the coefficients after Morgenthal */
	void compute_influence_matrix();

	/** @name Vortex sheet related methods */
	///@{

	/// Initialise the vortex sheet
	void form_vortex_sheet();

	/// Solve for the vortex sheet
	void solvevortexsheet(VortexBlobs &blobs);
	///@}

	/// Convect point vortices
	void convect();

	/** @name Methods related to diffusion */
	///@{

	/// Move the point vortices using the random walk
	void diffrw();

	/// Generate vorticity on the boundary and release into the flow
	void diffuse_vs_rw();
	///@}

	/// Reflect particles generated released inside the body outside
	void reflect();

	/// Write the output file
	/** Happens at each time step */
	void write_outputs();

	/** @name Timeloop control */
	///@{

	/// Increment the time step
	void increment_step();

	/// Get the size of the vortex sheet
	unsigned get_vs_size();

	/// Get number of total timesteps
	unsigned get_steps();

	/// Simulation time of the current time step
	double get_time();
	///@}

	/// Total number of point vortices at a given time
	unsigned get_size();

	/// Get the velocities at a given point
	void probe_velocities();

	/// Solve the particular problem - contains the timeloop
	void solve();

	private:
	/// Mirror a particle from one side of the boundary to the other
	std::vector<double> mirror(double x_init,
	                           double z_init,
	                           double x_0,
	                           double z_0,
	                           double x_1,
	                           double z_1);

	/// Determine if a particle lies inside the solid body
	int inside_body(double x, double z);

	/// Computes a particular timestep
	void compute_step();

	/// This is all the useless stuff below //////////////////////

	/*
	// read paramters
	void init_dipole(double zdistance,
	                 double xdistance,
	                 unsigned nvb1,
	                 unsigned nvb2,
	                 double radius1,
	                 double radius2);

	// Vortex sheet stuff
	void imagebc();

	// Diffusion stuff
	void diffpse(double nu, double dt);

	// Housekeeping stuff
	void delete_vort();
	void absorb();

	// Diagnostics
	void compute_strouhal();

	std::vector<double> uvs, wvs, uI, wI;
	// std::vector<double> xi_rw, eta_rw; // Why were these here???
	// std::vector<double> m_Gamma_abs;
	// double m_GammaDel,
	// double m_St;
	*/
};

#endif
