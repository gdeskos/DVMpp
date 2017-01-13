#include "VortexSheet.hpp"

VortexSheet::VortexSheet()
{
}

void VortexSheet::resize(unsigned size)
{
	gamma.resize(size);
	gamma_prev.resize(size);
	x.resize(size);
	z.resize(size);
	xc.resize(size);
	zc.resize(size);
	theta.resize(size);
	ds.resize(size);
	enx.resize(size);
	enz.resize(size);
	etx.resize(size);
	etz.resize(size);
}

unsigned VortexSheet::size()
{
	if ((gamma.size() != xc.size()) && (gamma_prev.size() != xc.size())
	    && (gamma.size() != zc.size())
	    && (gamma.size() != ds.size())
	    && (gamma.size() != theta.size())
	    && (gamma.size() != enx.size())
	    && (gamma.size() != enz.size())
	    && (gamma.size() != etx.size())
	    && (gamma.size() != etz.size())) {
		throw std::string("Size mismatch in VortexSheet");
	}
	return gamma.size();
}

void VortexSheet::print_collocation()
{
	for (unsigned i = 0; i < size(); i++) {
		std::cout << "(xc,zc) = (" << xc[i] << "," << zc[i] << ")" << std::endl;
	}
}

void VortexSheet::print_unit_vectors()
{
	for (unsigned i = 0; i < size(); i++) {
		std::cout << "For panel No " << i << "the normal vectors are en = ("
		          << enx[i] << "," << enz[i] << ") and et = (" << etx[i] << ","
		          << etz[i] << ")" << std::endl;
	}
}

void VortexSheet::print_gamma()
{
	for (unsigned i = 0; i < size(); i++) {
		std::cout << gamma[i] << std::endl;
	}
}
