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
	gamma.set_size(size);
	x.set_size(size+1);
	z.set_size(size+1);
	xc.set_size(size);
	zc.set_size(size);
	theta.set_size(size);
	ds.set_size(size);
	enx.set_size(size);
	enz.set_size(size);
	etx.set_size(size);
	etz.set_size(size);
}

void VortexSheet::compute_loads(double Urel)
{
	Vector P(size());
	Vector Dp = -m_rho / m_dt * (gamma % ds);

    for (unsigned i=0;i<size()-1;i++)
    {
        P(i+1)=P(i)+Dp(i);
    }
   
    // Finding the average of the two pressures at the boundary
    for (unsigned i=0;i<size()-1;i++)
    {
        P(i)=0.5*(P(i)+P(i+1));
    }

    P(size()-1)=0.5*(P(size()-1)+P(0));
	double pmax=arma::max(P);
    double pref=0.5*m_rho*Urel*Urel;

    P += (pref-pmax)*arma::ones(size(),1); 

    m_fx = -arma::sum(P % enx % ds);
	m_fz = -arma::sum(P % enz % ds);

}

unsigned VortexSheet::size()
{
	if ((gamma.size() != xc.size())
	    && (gamma.size() != zc.size())
	    && (gamma.size() != x.size()-1)
	    && (gamma.size() != z.size()-1)
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
		std::cout << "(xc,zc) = (" << xc(i) << "," << zc(i) << ")" << std::endl;
	}
}

void VortexSheet::print_unit_vectors()
{
	for (unsigned i = 0; i < size(); i++) {
		std::cout << "For panel No " << i << "the normal vectors are en = ("
		          << enx(i) << "," << enz(i) << ") and et = (" << etx(i) << ","
		          << etz(i) << ")" << std::endl;
	}
}

void VortexSheet::print_gamma()
{
	gamma.print();
}

std::tuple<double, double> VortexSheet::get_forces()
{
	return std::make_tuple(m_fx, m_fz);
}
