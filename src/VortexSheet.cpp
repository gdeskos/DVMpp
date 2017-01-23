#include "VortexSheet.hpp"

VortexSheet::VortexSheet()
{
}

VortexSheet::VortexSheet(const XmlHandler &xml)
{
	m_rho = xml.getValueAttribute("constants", "density");
	m_dt = xml.getValueAttribute("time", "dt");
}


void VortexSheet::resize(unsigned size)
{
	gamma.resize(size);
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

void VortexSheet::compute_loads()
{
	double dp;

	std::vector<double> p;
	p.resize(size());

	// not so sure here (reference pressure)
	p[0] = 0;
	for (unsigned i = 1; i < size(); i++) {

		dp = m_rho / m_dt * gamma[i] * ds[i];

		p[i] = p[0] - dp;
	}

	m_fx = 0;
	m_fz = 0;
	for (unsigned i = 0; i < size(); i++) {
		m_fx += p[i] * enx[i] * ds[i];
		m_fz += p[i] * enz[i] * ds[i];
	}

	std::cout << "Fx = " << m_fx << " Fz = " << m_fz << std::endl;
}

unsigned VortexSheet::size()
{
	if ((gamma.size() != xc.size())
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

std::tuple<double, double> VortexSheet::get_forces()
{
	return std::make_tuple(m_fx, m_fz);
}
