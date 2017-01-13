#include "Body.hpp"

#include <iostream>

Body::Body()
{
}

void Body::resize(unsigned size)
{
	x.resize(size);
	z.resize(size);
}

unsigned Body::size()
{
	if (x.size() != z.size()) {
		throw std::string("Size mismatch in Body");
	}
	return x.size();
}

void Body::print_location()
{
	for (unsigned i = 0; i < size(); i++) {
		std::cout << " Body coor (x,z) = (" << x[i] << "," << z[i] << ")"
		          << std::endl;
	}
}
