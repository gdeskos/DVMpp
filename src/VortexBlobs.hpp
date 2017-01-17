#ifndef VORTEXBLOBS_H
#define VORTEXBLOBS_H

#include "BaseTypes.hpp"

#include <vector>

/// VortexBlobs describes the vortex blobs

class VortexBlobs
{
	public:
	std::vector<unsigned> ID;  ///< Vortex blob ID
	Vector x;     ///< x-coordinate
	Vector z;     ///< z-coordinate
	Vector circ;  ///< circulation
	Vector sigma; ///< vortex cut-off
	Vector u;     ///< local x-velocity
	Vector w;     ///< local z-velocity
	Vector uvs;
	Vector wvs;
	Vector omega; ///< Local vorticity

	public:
	VortexBlobs();

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
