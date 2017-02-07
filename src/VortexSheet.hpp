#ifndef VORTEXSHEET_H
#define VORTEXSHEET_H

#include "BaseTypes.hpp"
#include "XmlHandler.hpp"
#include "VortexBlobs.hpp"
#include "Random.hpp"
#include <string>
#include <sstream>
#include <string.h>
#include <fstream>
#include <iostream>
#include <tuple>
#include <cassert>

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
    double m_nu;  ///< Kinematic viscosity of the fluid
	double m_dt;  ///< Timestep
    double m_maxGamma; ///< maximum circulation of the nascent vortex
    unsigned m_maxNumPanelVort; /// Maximum number of vortices allowed in each panel
    double   m_cutoff_exp;      /// Cutoff distance coefficient q, according to Perlman 1985 0.5<q<1 
	double m_fx; ///< Force in x-direction
	double m_fz; ///< Force in z-direction
	

	public:
	VortexSheet();

	VortexSheet(const XmlHandler &xml);

	/// Resize all of the data members
	/** \param size New size */
	void resize(unsigned size);

    /// Release nascent vortices from the vortexsheet using random walk
    VortexBlobs release_nascent_vortices_rw(Random& _rand);
   
    /// Efficient Surface Algorithm for Random Walk (Smith and Stansby 1989, JCP paper)

    /// Reflects the vortices the the image location from the panels
    void reflect(VortexBlobs &vortex);

    /// Find Inside the vortex sheet vortices
    int inside_body(const double &xcoor, const double &zcoor);

    /// Mirror the particle from a panel 
    Vector mirror(const double &x_init,
                  const double &z_init,
                  const double &x_0,
                  const double &z_0,
                  const double &x_1,
                  const double &z_1);

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

	/// Return the forces in the x and z directions
	std::tuple<double, double> get_forces();
};

#endif
