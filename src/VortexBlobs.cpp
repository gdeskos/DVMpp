#include "VortexBlobs.hpp"

#include <iostream>
#include <string>
//************************* Constructors *************************************//
VortexBlobs::VortexBlobs()
{
}

VortexBlobs::VortexBlobs(const XmlHandler &xml, const std::string &stamp)
{
	// Initialise the output files with the headers
	auto outdir = xml.getStringAttribute("io", "output_dir", true);
	auto dt = xml.getValueAttribute("time", "dt");
	auto steps = xml.getValueAttribute("time", "steps");

	m_blobsfile = outdir + stamp + std::string("_vortex.dat");
	auto bf = std::ofstream(m_blobsfile);
	bf << size() << " # Number of nodes\n"
	   << dt << " # Time step\n"
	   << steps << " # Steps\n"
	   << "Time [s]\tx-position [m]\tz-position [m]\tcirculation\n";

	m_numfile = outdir + stamp + std::string("_vortex_num.dat");
	auto nf = std::ofstream(m_numfile);
	nf << "Time [s]\tNumber of vortices\n";
}

//****************************** Public Methods ******************************//
VortexBlobs::VortexBlobs(const unsigned &N)
{
	m_ID.resize(N);
	m_x.set_size(N);
	m_z.set_size(N);
	m_circ.set_size(N);
	m_sigma.set_size(N);
	m_u.set_size(N);
	m_w.set_size(N);
	m_uvs.set_size(N);
	m_wvs.set_size(N);
}

void VortexBlobs::append_vortices(VortexBlobs &NewVortBlobs)
{
	// Get the number of Vortices
	auto Nold = size();
	auto Nnew = NewVortBlobs.size();
	for (unsigned i = 0; i < Nnew; i++) {
		m_ID.push_back(Nold + NewVortBlobs.m_ID[i]);
	}
	m_x.insert_rows(m_x.n_elem, NewVortBlobs.m_x);
	m_z.insert_rows(m_z.n_elem, NewVortBlobs.m_z);
	m_circ.insert_rows(m_circ.n_elem, NewVortBlobs.m_circ);
	m_sigma.insert_rows(m_sigma.n_elem, NewVortBlobs.m_sigma);
	m_u.insert_rows(m_u.n_elem, NewVortBlobs.m_u);
	m_w.insert_rows(m_w.n_elem, NewVortBlobs.m_w);
	m_uvs.insert_rows(m_uvs.n_elem, NewVortBlobs.m_uvs);
	m_wvs.insert_rows(m_wvs.n_elem, NewVortBlobs.m_wvs);
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

				threshold = 10.0 * std::pow(m_sigma(j), 2.0);
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
		m_u(i) *= math::rpi2;
		m_w(i) *= math::rpi2;
	}
}

void VortexBlobs::diffusion_random_walk(Random &_rand, double nu, double dt)
{
	double R1, R2, rrw, thetarw;

	for (unsigned i = 0; i < size(); i++) {

		// Generate two random numbers in the range 0...1
		R1 = _rand.rand();
		R2 = _rand.rand();

		// Calculate r and theta for the random walk
		rrw = std::sqrt(4.0 * nu * dt * std::log(1.0 / R1));
		thetarw = 2.0 * math::pi * R2;

		m_x(i) += rrw * cos(thetarw);
		m_z(i) += rrw * sin(thetarw);
	}
}

double VortexBlobs::totalcirc()
{
	return arma::sum(m_circ);
}

unsigned VortexBlobs::size() const
{
	auto ID_sz = m_ID.size();
	if ((ID_sz != m_x.size()) && (ID_sz != m_z.size())
	    && (ID_sz != m_circ.size())
	    && (ID_sz != m_sigma.size())
	    && (ID_sz != m_u.size())
	    && (ID_sz != m_w.size())
	    && (ID_sz != m_uvs.size())
	    && (ID_sz != m_wvs.size())) {
		throw std::string("Size mismatch in VortexBlobs");
	}
	return ID_sz;
}

void VortexBlobs::print_location()
{
	for (unsigned i = 0; i < size(); i++) {
		std::cout << " The location of Vortex Blob " << m_ID[i]
		          << " is (x,z) = (" << m_x(i) << "," << m_z(i) << ")"
		          << std::endl;
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

void VortexBlobs::write_step(double time, unsigned step)
{
	if (m_blobsfile.empty() || m_numfile.empty()) {
		throw std::string(
		    "VortexBlobs::write_step -> must set filenames before writing");
	}

	// Make sure we open up to end of the file
	auto bf = std::ofstream(m_blobsfile, std::ios_base::app);

	for (unsigned i = 0; i < size(); ++i) {
		bf << time << "\t" << m_x(i) << "\t" << m_z(i) << "\t" << m_circ(i)
		   << "\n";
	}

	auto nf = std::ofstream(m_numfile, std::ios_base::app);
	nf << step << "\t" << size() << "\n";
}

//****************************** Destructor **********************************//
VortexBlobs::~VortexBlobs()
{
	// Destructor
}
