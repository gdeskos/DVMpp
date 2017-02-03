#include "VortexBlobs.hpp"

#include <iostream>
#include <string>
//************************* Constructors *************************************//
VortexBlobs::VortexBlobs()
{
	// Don't put anything in here, it is never called!
}

VortexBlobs::VortexBlobs(const XmlHandler &xml)
{
	m_pi = 4.0 * atan(1.0);
	m_rpi2 = 1.0 / (2.0 * m_pi);

}


// ********************************* Public Methods *****************************//
VortexBlobs::VortexBlobs(const unsigned &N)
{
	m_pi = 4.0 * atan(1.0);
	m_rpi2 = 1.0 / (2.0 * m_pi);
	
    resize(N);
}

void VortexBlobs::resize(unsigned size)
{
	m_ID.resize(size);
	m_x.set_size(size);
	m_z.set_size(size);
	m_circ.set_size(size);
	m_sigma.set_size(size);
	m_u.set_size(size);
	m_w.set_size(size);
	m_uvs.set_size(size);
	m_wvs.set_size(size);
	m_omega.set_size(size);
}

void VortexBlobs::append_vortices(VortexBlobs& NewVortBlobs) 
{
    // Get the number of Vortices
    auto Nold=size();
    auto Nnew=NewVortBlobs.size();
    for (unsigned i=0; i<Nnew;i++)
    {
    m_ID.push_back(Nold+NewVortBlobs.m_ID[i]);
    }
    m_x.insert_rows(m_x.n_elem, NewVortBlobs.m_x);
	m_z.insert_rows(m_z.n_elem, NewVortBlobs.m_z);
	m_circ.insert_rows(m_circ.n_elem, NewVortBlobs.m_circ);
	m_sigma.insert_rows(m_sigma.n_elem, NewVortBlobs.m_sigma); 
	m_u.insert_rows(m_u.n_elem, NewVortBlobs.m_u);
	m_w.insert_rows(m_w.n_elem, NewVortBlobs.m_w);
	m_uvs.insert_rows(m_uvs.n_elem, NewVortBlobs.m_uvs);
	m_wvs.insert_rows(m_wvs.n_elem, NewVortBlobs.m_wvs);
	m_omega.insert_rows(m_omega.n_elem, NewVortBlobs.m_omega); 
}
void VortexBlobs::biotsavart()
{

#pragma omp parallel for
	for (unsigned i = 0; i < size(); i++) {

		double dx_ij, dz_ij, dK_ij, dr_ij2, threshold, rsigmasqr;

		m_u(i) = 0.0;
		m_w(i) = 0.0;

		for (unsigned j = 0; j < size(); j++) {

			if (i != j) {
				dx_ij = m_x(i) - m_x(j);
				dz_ij = m_z(i) - m_z(j);
				dr_ij2 = std::pow(dx_ij, 2) + std::pow(dz_ij, 2.0);

				threshold =10.0*std::pow(m_sigma(j), 2.0);
				rsigmasqr = 1.0 / std::pow(m_sigma(j), 2.0);

				if (dr_ij2 < threshold) {
					dK_ij = (1.0 - std::exp(-dr_ij2 * rsigmasqr)) / dr_ij2;
				} else {
					dK_ij = 1.0 / dr_ij2;
				}

				m_u(i) -= dK_ij * dz_ij * m_circ(j);
				m_w(i) += dK_ij * dx_ij * m_circ(j);
			}
		}
		m_u(i) *= m_rpi2;
		m_w(i) *= m_rpi2;
	}
}

double VortexBlobs::totalcirc()
{
	return arma::sum(m_circ);
}

unsigned VortexBlobs::size()
{
	auto ID_sz = m_ID.size();
	if ((ID_sz != m_x.size()) && (ID_sz != m_z.size())
	    && (ID_sz != m_circ.size())
	    && (ID_sz != m_sigma.size())
	    && (ID_sz != m_u.size())
	    && (ID_sz != m_w.size())
	    && (ID_sz != m_uvs.size())
	    && (ID_sz != m_wvs.size())
	    && (ID_sz != m_omega.size())) {
		throw std::string("Size mismatch in VortexBlobs");
	}
	return ID_sz;
}

void VortexBlobs::print_location()
{
	for (unsigned i = 0; i < size(); i++) {
		std::cout << " The location of Vortex Blob " << m_ID[i] << " is (x,z) = ("
		          << m_x(i) << "," << m_z(i) << ")" << std::endl;
	}
}

void VortexBlobs::print_velocity()
{
	for (unsigned i = 0; i < size(); i++) {
		std::cout << "(u,w) = (" << m_u(i) << "," << m_w(i) << ")" << std::endl;
	}
}

void VortexBlobs::print_circulation()
{
	for (unsigned i = 0; i < size(); i++) {
		std::cout << " circ = " << m_circ(i) << std::endl;
	}
}

//**************************************** Destructor *****************************************************//
VortexBlobs::~VortexBlobs()
{
    //Destructor
}

