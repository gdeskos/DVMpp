#ifndef VORTEXSHEET_H
#define VORTEXSHEET_H

#include "BaseTypes.hpp"
#include "Probe.hpp"
#include "Random.hpp"
#include "VortexBlobs.hpp"
#include "XmlHandler.hpp"
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>

/** \class VortexSheet
 * \brief Describes the vortex sheet
 * \file VortexSheet.hpp */

class VortexSheet
{
	private:
	Vector m_gamma; ///< surface vorticity
	Vector m_x;     ///< x-coordinate point
	Vector m_z;     ///< z-coordinate point
	Vector m_xc;    ///< x-collocation point
	Vector m_zc;    ///< z-collocation point
	Vector m_theta; ///< angle
	Vector m_ds;    ///< panel
	Vector m_enx;   ///< x-coor normal unit vector
	Vector m_enz;   ///< z-coor normal unit vector
	Vector m_etx;   ///< x-coor tangential unit vector
	Vector m_etz;   ///< z-coor tangential unit vector

	double m_A; ///< Reference area for drag coefficient

	Matrix m_infM; ///< Influence matrix

	double m_rho; ///< Density of fluid
	double m_nu;  ///< Kinematic viscosity of the fluid
	double m_dt;  ///< Timestep

	/// Maximum number of vorticies allow in each panel
	unsigned m_maxNumPanelVort;
	/// Cutoff distance coefficient q, according to Perlman 1985 0.5<q<1
	double m_cutoff_exp;

	double m_fx; ///< Force in x-direction
	double m_fz; ///< Force in z-direction

	double m_Ux; ///< Freestream x-velocity
	double m_Uz; ///< Freestream z-velocity

	std::ofstream m_gammafile; ///< Stream for gamma output file
	std::ofstream m_forcefile; ///< Stream for force output file

	/// Read the input coordinate file
	/** \param domain domain file to be read */
	void read_input_coord(std::string domain);

	/// Form the vortex sheet from the body coordinates
	void form_vortex_sheet();

	/// Compute the influence matrix
	/** Computed using the coefficients after Morgenthal */
	void compute_influence_matrix();

	/// Find Inside the vortex sheet vortices
	bool inside_body(double xcoor, double zcoor);

	/// Mirror the particle from a panel
	Vector mirror(double x_init,
	              double z_init,
	              double x_0,
	              double z_0,
	              double x_1,
	              double z_1);

	public:
	VortexSheet();

	VortexSheet(const XmlHandler &xml, const std::string &timestamp);

	/// Impose the vortex sheet boundary condition
	/** Updates the blobs in place.
	 * \param blobs vortex blobs to update */
	void vortexsheetbc(VortexBlobs &blobs);

	/// Solve for the vortex sheet
	void solvevortexsheet(VortexBlobs &blobs);

	/// Release nascent vortices from the vortexsheet using random walk
	VortexBlobs release_nascent_vortices_rw(Random &_rand);

	/// Efficient Surface Algorithm for Random Walk (Smith and Stansby 1989, JCP
	/// paper)

	/// Reflects the vortices the the image location from the panels
	void reflect(VortexBlobs &vortex);

	/// Compute the forces on the body
	void compute_loads(double Urel);

	/// Number of vortex sheets
	unsigned size();

	/// Print location of each collocation point
	void print_collocation();

	/// Print the unit vector for each collocation pont
	void print_unit_vectors();

	/// Print the surface vorticity for each collocation point
	void print_gamma();

	/// Print the relevant quantities to file
	/** \param time Time step (s) */
	void write_step(double time);

	/// Return the forces in the x and z directions
	std::tuple<double, double> get_forces();

	/// Calculate the velocity at the probe points
	/** Updates the probe in place */
	void probe_velocities(const VortexBlobs &blobs, Probe &probe);
};

#endif
