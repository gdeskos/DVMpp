#ifndef VORTEXSHEET_H
#define VORTEXSHEET_H

#include "BaseTypes.hpp"
#include "XmlHandler.hpp"
#include "VortexBlobs.hpp"
#include <iostream>
#include <tuple>

/// VortexSheet describes the vortex sheet

class VortexSheet
{
	public:
	Vector gamma;      ///< surface vorticity
	Vector x;          ///< x-coordinate point
	Vector z;          ///< z-coordinate point
	Vector xc;         ///< x-collocation point
	Vector zc;         ///< z-collocation point
	Vector theta;      ///< angle
	Vector ds;         ///< panel
	Vector enx;        ///< x-coor normal unit vector
	Vector enz;        ///< z-coor normal unit vector
	Vector etx;        ///< x-coor tangential unit vector
	Vector etz;        ///< z-coor tangential unit vector

	private:
	double m_rho; ///< Density of fluid
	double m_dt;  ///< Timestep
	double m_fx; ///< Force in x-direction
	double m_fz; ///< Force in z-direction

	double m_pi;   ///< pi
	double m_rpi2; ///< 1 / (2pi)

	public:
	VortexSheet();

	VortexSheet(const XmlHandler &xml);

	/// Resize all of the data members
	/** \param size New size */
	void resize(unsigned size);

	/// Impose the vortex sheet boundary condition
	/** Updates the blobs in place.
	 * \param blobs vortex blobs to update */
	void vortexsheetbc(VortexBlobs &blobs);

	/// Compute the forces on the body
	void compute_loads();

	/// Number of vortex sheets
	unsigned size();

	/// Print location of each collocation point
	void print_collocation();

	/// Print the unit vector for each collocation pont
	void print_unit_vectors();

	/// Print the surface vorticity for each collocation point
	void print_gamma();

	/// Return the forces in the x and z directions
	std::tuple<double, double> get_forces();
};

#endif
