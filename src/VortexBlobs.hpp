#ifndef VORTEXBLOBS_H
#define VORTEXBLOBS_H

#include "BaseTypes.hpp"
#include "XmlHandler.hpp"

#include <vector>

/// VortexBlobs describes the vortex blobs

class VortexBlobs
{
	public:
	std::vector<unsigned> m_ID;  ///< Vortex blob ID
	Vector m_x;     ///< x-coordinate
	Vector m_z;     ///< z-coordinate
	Vector m_circ;  ///< circulation
	Vector m_sigma; ///< vortex cut-off
	Vector m_u;     ///< local x-velocity
	Vector m_w;     ///< local z-velocity
	Vector m_uvs;
	Vector m_wvs;
	Vector m_omega; ///< Local vorticity

	double m_pi;   ///< pi
	double m_rpi2; ///< 1 / (2pi)

	double m_kernel_threshold; ///< Threshold for the kernel

	public:
	VortexBlobs();

	/// Constructor
	/** \param xml input xml file*/
	VortexBlobs(const XmlHandler &xml);

	/// Biot-Savart relationship
	void biotsavart();

	/// Find the total circulation
	double totalcirc();

	/// Resize all of the data members
	/** \param size New size */
	void resize(unsigned size);

	/// Current number of vortex blobs
	unsigned size();

	/// Print the location of each vortex blob
	void print_location();

	/// Print the velocity of each vortex blob
	void print_velocity();

	/// Print the circulation of each vortex blob
	void print_circulation();
};

#endif
