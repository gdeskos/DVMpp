#ifndef VORTEXSHEET_H
#define VORTEXSHEET_H

#include "XmlHandler.hpp"
#include <iostream>
#include <vector>
#include <tuple>

/// VortexSheet describes the vortex sheet

class VortexSheet
{
	public:
	std::vector<double> gamma;      ///< surface vorticity
	std::vector<double> x;          ///< x-coordinate point
	std::vector<double> z;          ///< z-coordinate point
	std::vector<double> xc;         ///< x-collocation point
	std::vector<double> zc;         ///< z-collocation point
	std::vector<double> theta;      ///< angle
	std::vector<double> ds;         ///< panel
	std::vector<double> enx;        ///< x-coor normal unit vector
	std::vector<double> enz;        ///< z-coor normal unit vector
	std::vector<double> etx;        ///< x-coor tangential unit vector
	std::vector<double> etz;        ///< z-coor tangential unit vector

	private:
	double m_rho; ///< Density of fluid
	double m_dt;  ///< Timestep
	double m_fx; ///< Force in x-direction
	double m_fz; ///< Force in z-direction

	public:
	VortexSheet();

	VortexSheet(const XmlHandler &xml);

	/// Resize all of the data members
	/** \param size New size */
	void resize(unsigned size);

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
