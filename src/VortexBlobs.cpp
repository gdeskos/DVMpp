#include "VortexBlobs.hpp"

#include <iostream>
#include <string>

VortexBlobs::VortexBlobs()
{
}

void VortexBlobs::resize(unsigned size)
{
	ID.resize(size);
	x.set_size(size);
	z.set_size(size);
	circ.set_size(size);
	sigma.set_size(size);
	u.set_size(size);
	w.set_size(size);
	uvs.set_size(size);
	wvs.set_size(size);
	omega.set_size(size);
}

double VortexBlobs::totalcirc()
{
	return arma::sum(circ);
}

unsigned VortexBlobs::size()
{
	if ((ID.size() != x.size()) && (ID.size() != z.size())
	    && (ID.size() != circ.size())
	    && (ID.size() != sigma.size())
	    && (ID.size() != u.size())
	    && (ID.size() != w.size())
	    && (ID.size() != uvs.size())
	    && (ID.size() != wvs.size())
	    && ID.size() != omega.size()) {
		throw std::string("Size mismatch in VortexBlobs");
	}
	return ID.size();
}

void VortexBlobs::print_location()
{
	for (unsigned i = 0; i < size(); i++) {
		std::cout << " The location of Vortex Blob " << ID[i] << " is (x,z) = ("
		          << x(i) << "," << z(i) << ")" << std::endl;
	}
}

void VortexBlobs::print_velocity()
{
	for (unsigned i = 0; i < size(); i++) {
		std::cout << "(u,w) = (" << u(i) << "," << w(i) << ")" << std::endl;
	}
}

void VortexBlobs::print_circulation()
{
	for (unsigned i = 0; i < size(); i++) {
		std::cout << " circ = " << circ(i) << std::endl;
	}
}
