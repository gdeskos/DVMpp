#ifndef BODY_H
#define BODY_H

#include <vector>

/// Body describes the solid body

class Body
{
	public:
	std::vector<double> x; ///< x-coordinate
	std::vector<double> z; ///< z-coordinate

	public:
	Body();

	/// Resize all the date members
	/** \param size New size */
	void resize(unsigned size);

	/// Number of points on the boundary
	unsigned size();

	/// Print the location each boundary point
	void print_location();
};

#endif
