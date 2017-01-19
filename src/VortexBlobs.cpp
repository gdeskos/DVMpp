#include "VortexBlobs.hpp"

#include <iostream>
#include <string>

VortexBlobs::VortexBlobs()
{
	// Don't put anything in here, it is never called!
}

VortexBlobs::VortexBlobs(const XmlHandler &xml)
{
	m_pi = 4.0 * atan(1.0);
	m_rpi2 = 1.0 / (2.0 * m_pi);

	m_kernel_threshold = xml.getValueAttribute("constants", "kernel_threshold");
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

				threshold =
				    m_kernel_threshold * std::pow(m_sigma(j), 2.0);
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
