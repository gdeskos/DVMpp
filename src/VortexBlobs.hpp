#ifndef VORTEXBLOBS_H
#define VORTEXBLOBS_H

#include "BaseTypes.hpp"
#include "Random.hpp"
#include "XmlHandler.hpp"

#include <vector>

/** \class VortexBlobs
 * \brief Describes the vortex blobs
 * \file VortexBlobs.hpp */

class VortexBlobs
{
	public:
	std::vector<unsigned> m_ID; ///< Vortex blob ID
	Vector m_x;                 ///< x-coordinate
	Vector m_z;                 ///< z-coordinate
	Vector m_circ;              ///< circulation
	Vector m_sigma;             ///< vortex cut-off
	Vector m_u;                 ///< local x-velocity
	Vector m_w;                 ///< local z-velocity
	Vector m_uvs;
	Vector m_wvs;

	private:
	/// Output filename for the blob data;
	/** Needs to hold a string as a stream cannot be copy-constructed. This
	 * means we cannot return a VortexBlobs instance from a function. */
	std::string m_blobsfile;

	/// Output filename for the number of blobs at each timestep
	/** Needs to hold a string as a stream cannot be copy-constructed. This
	 * means we cannot return a VortexBlobs instance from a function. */
	std::string m_numfile;

	public:
	VortexBlobs();

	/// Create vortex blobs instance that will write to file
	VortexBlobs(const XmlHandler &xml, const std::string &stamp);

	/// Constructor
	/** \creates N number of vortices and sets them to zero*/
	VortexBlobs(const unsigned &N);

	/// Appends vortexblobs
	void append_vortices(VortexBlobs &NewVortBlobs);

	/// Biot-Savart relationship
	void biotsavart();

	/// Diffusion through random walk
	void diffusion_random_walk(Random &_rand, double nu, double dt);

	/// Find the total circulation
	double totalcirc();

	/// Current number of vortex blobs
	unsigned size() const;

	/// Print the location of each vortex blob
	void print_location();

	/// Print the velocity of each vortex blob
	void print_velocity();

	/// Print the circulation of each vortex blob
	void print_circulation();

	/// Write the blob info and number of vortices to file
	/** \param time Simulation time [s] */
	void write_step(double time, unsigned step);

	/// Destructor
	~VortexBlobs();
};

#endif
