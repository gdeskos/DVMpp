#include "VortexSheet.hpp"

VortexSheet::VortexSheet()
{
	// Don't put anything in here, it is never called
}

VortexSheet::VortexSheet(const XmlHandler &xml)
{
	m_rho = xml.getValueAttribute("constants", "density");
	m_dt = xml.getValueAttribute("time", "dt");

	m_pi = 4.0 * atan(1.0);
	m_rpi2 = 1.0 / (2.0 * m_pi);

	auto indir = xml.getStringAttribute("io", "input_dir");
	// Make sure we have trailing separator
	if (indir.back() != '/') {
		indir += '/';
	}

	std::string domain = indir + xml.getStringAttribute("io", "domain_file");

	read_input_coord(domain);
	form_vortex_sheet();
}

void VortexSheet::read_input_coord(std::string file)
{
	// Make use of the armadillo functionality to read the file
	Matrix coords;
	bool status = coords.load(file, arma::raw_ascii);

	if (status != true) {
		throw std::string("Could not open coordinates file");
	}

	m_x = coords.col(0);
	m_z = coords.col(1);

	std::cout << "Succesfully loaded coordinate file with " << m_x.n_elem
	          << " points." << std::endl;
}

void VortexSheet::form_vortex_sheet()
{
	// Create the collocation points
	auto end = m_x.n_elem - 1;
	m_xc = 0.5 * (m_x.rows(0, end - 1) + m_x.rows(1, end));
	m_zc = 0.5 * (m_z.rows(0, end - 1) + m_z.rows(1, end));

	// ds and theta along the vortex sheet
	Vector dx = m_x.rows(1, end) - m_x.rows(0, end - 1);
	Vector dz = m_z.rows(1, end) - m_z.rows(0, end - 1);

	m_ds = arma::sqrt(arma::pow(dx, 2) + arma::pow(dz, 2));
	m_theta = arma::atan2(dz, dx);

	// Outwards facing normals and tangentials
	m_enx = dz / m_ds;
	m_enz = -dx / m_ds;
	m_etx = -m_enz;
	m_etz = m_enx;

	// Inwards facing normals and tangentials
	// m_etx = arma::cos(m_theta);
	// m_etz = arma::sin(m_theta);
	// m_enx = -m_etz;
	// m_enz = m_etx;

	// Make sure the surface vorticity is the correct size
	m_gamma.copy_size(m_xc);

	std::cout << "Created vortex sheet of size " << m_gamma.n_elem << std::endl;
}

void VortexSheet::resize(unsigned size)
{
	m_gamma.set_size(size);
	m_xc.set_size(size);
	m_zc.set_size(size);
	m_theta.set_size(size);
	m_ds.set_size(size);
	m_enx.set_size(size);
	m_enz.set_size(size);
	m_etx.set_size(size);
	m_etz.set_size(size);
}

void VortexSheet::vortexsheetbc(VortexBlobs &blobs)
{
	double c1, c2, c5, c6, c7, c8, c9;
	Matrix px, py, qx, qy;

	px.set_size(blobs.size(), size());
	qx.copy_size(px);
	py.copy_size(px);
	qy.copy_size(px);

// compute the coefficients
	for (unsigned i = 0; i < blobs.size(); i++) {

		double xi = blobs.m_x(i);
		double zi = blobs.m_z(i);

		for (unsigned j = 0; j < size(); j++) {

			double xj = m_xc(j);
			double zj = m_zc(j);
			double thetaj = m_theta(j);
			double dsj = m_ds(j);

			c1 = -(xi - xj) * cos(thetaj) - (zi - zj) * sin(thetaj);
			c2 = std::pow(xi - xj, 2.0) + std::pow(zi - zj, 2.0);
			c5 = (xi - xj) * sin(thetaj) - (zi - zj) * cos(thetaj);
			c6 = log(1.0 + dsj * ((dsj + 2 * c1) / c2));
			c7 = atan2(c5 * dsj, c2 + c1 * dsj);
			c8 =
			    (xi - xj) * sin(-2.0 * thetaj) + (zi - zj) * cos(-2.0 * thetaj);
			c9 =
			    (xi - xj) * cos(-2.0 * thetaj) - (zi - zj) * sin(-2.0 * thetaj);

			qx(i, j) = -sin(thetaj) + 0.5 * c8 * c6 / dsj
			           + (c1 * cos(thetaj) + c5 * sin(thetaj)) * c7 / dsj;
			px(i, j) = -0.5 * c6 * sin(thetaj) - c7 * cos(thetaj) - qx(i, j);

			qy(i, j) = cos(thetaj) + 0.5 * c9 * c6 / dsj
			           + (c1 * sin(thetaj) - c5 * cos(thetaj)) * c7 / dsj;
			py(i, j) = 0.5 * c6 * cos(thetaj) - c7 * sin(thetaj) - qy(i, j);
		}
	}

	// calculate the vortexsheet induced velocities
	unsigned last = size() - 1;

#pragma omp parallel for
	for (unsigned i = 0; i < blobs.size(); i++) {

		blobs.m_uvs(i) = 0.0;
		blobs.m_wvs(i) = 0.0;

		for (unsigned j = 0; j < size(); j++) {
			if (j == last) {
				blobs.m_uvs(i) += px(i, last) * m_gamma(last)
				                   + qx(i, last) * m_gamma(0);
				blobs.m_wvs(i) += py(i, last) * m_gamma(last)
				                   + qy(i, last) * m_gamma(0);
			} else {
				blobs.m_uvs(i) += px(i, j) * m_gamma(j)
				                   + qx(i, j) * m_gamma(j + 1);
				blobs.m_wvs(i) += py(i, j) * m_gamma(j)
				                   + qy(i, j) * m_gamma(j + 1);
			}
		}
		blobs.m_uvs(i) *= m_rpi2;
		blobs.m_wvs(i) *= m_rpi2;
	}
}

void VortexSheet::compute_loads(double Urel)
{
	Vector P(size());
	Vector Dp = -m_rho / m_dt * (m_gamma % m_ds);

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

    m_fx = -arma::sum(P % m_enx % m_ds);
	m_fz = -arma::sum(P % m_enz % m_ds);

}

unsigned VortexSheet::size()
{
	if ((m_gamma.size() != m_xc.size())
	    && (m_gamma.size() != m_zc.size())
	    && (m_gamma.size() != m_x.size()-1)
	    && (m_gamma.size() != m_z.size()-1)
	    && (m_gamma.size() != m_ds.size())
	    && (m_gamma.size() != m_theta.size())
	    && (m_gamma.size() != m_enx.size())
	    && (m_gamma.size() != m_enz.size())
	    && (m_gamma.size() != m_etx.size())
	    && (m_gamma.size() != m_etz.size())) {
		throw std::string("Size mismatch in VortexSheet");
	}
	return m_gamma.size();
}

void VortexSheet::print_collocation()
{
	for (unsigned i = 0; i < size(); i++) {
		std::cout << "(xc,zc) = (" << m_xc(i) << "," << m_zc(i) << ")" << std::endl;
	}
}

void VortexSheet::print_unit_vectors()
{
	for (unsigned i = 0; i < size(); i++) {
		std::cout << "For panel No " << i << "the normal vectors are en = ("
		          << m_enx(i) << "," << m_enz(i) << ") and et = (" << m_etx(i) << ","
		          << m_etz(i) << ")" << std::endl;
	}
}

void VortexSheet::print_gamma()
{
	m_gamma.print();
}

std::tuple<double, double> VortexSheet::get_forces()
{
	return std::make_tuple(m_fx, m_fz);
}
