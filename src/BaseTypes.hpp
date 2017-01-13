#ifndef BASETYPES
#define BASETYPES

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

typedef std::vector<double> Vector;
typedef std::vector<Vector> Matrix;
typedef std::vector<unsigned>::iterator IndexIt;


class VortexSheet
{
	public:
	VortexSheet();
	void resize(unsigned size);
	unsigned size();
	void print_collocation();
	void print_unit_vectors();
	void print_gamma();

	public:
	std::vector<double> gamma;      // surface vorticity
	std::vector<double> gamma_prev; // surface vorticity in previous time step
	                                // (for load calculations)
	std::vector<double> x;          // x-coordinate point
	std::vector<double> z;          // z-coordinate point
	std::vector<double> xc;         // x-collocation point
	std::vector<double> zc;         // z-collocation point
	std::vector<double> theta;      // angle
	std::vector<double> ds;         // panel
	std::vector<double> enx;        // x-coor normal unit vector
	std::vector<double> enz;        // z-coor normal unit vector
	std::vector<double> etx;        // x-coor tangential unit vector
	std::vector<double> etz;        // z-coor tangential unit vector
};

class Body
{
	public:
	Body();
	void resize(unsigned size);
	unsigned size();
	void print_location();

	public:
	std::vector<double> x; // x-coordinate
	std::vector<double> z; // z-coordinate
};
#endif
