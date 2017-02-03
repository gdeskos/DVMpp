#ifndef BODY_H
#define BODY_H

#include <vector>
#include "XmlHandler.hpp"

/// Body describes the solid body

class Body
{
	public:
	std::vector<double> x; ///< x-coordinate
	std::vector<double> z; ///< z-coordinate

	public:
	Body();

	Body(XmlHandler &xml);

	/// Resize all the date members
	/** \param size New size */
	void resize(unsigned size);

	/// Number of points on the boundary
	unsigned size();

	/// Read the input coordinates file
	void read_input_coord(std::string file);

	/// Print the location each boundary point
	void print_location();
};

#endif
