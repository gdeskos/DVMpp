#ifndef VORTEXBLOBS_H
#define VORTEXBLOBS_H

#include <vector>

/// VortexBlobs describes the vortex blobs

class VortexBlobs
{
	public:
	std::vector<unsigned> ID;  ///< Vortex blob ID
	std::vector<double> x;     ///< x-coordinate
	std::vector<double> z;     ///< z-coordinate
	std::vector<double> circ;  ///< circulation
	std::vector<double> sigma; ///< vortex cut-off
	std::vector<double> u;     ///< local x-velocity
	std::vector<double> w;     ///< local z-velocity
	std::vector<double> uvs;
	std::vector<double> wvs;
	std::vector<double> omega; ///< Local vorticity

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
