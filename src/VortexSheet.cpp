#include "VortexSheet.hpp"

VortexSheet::VortexSheet()
{
	// Don't put anything in here, it is never called
}

VortexSheet::VortexSheet(const XmlHandler &xml)
{
	m_rho = xml.getValueAttribute("constants", "density");
	m_nu = xml.getValueAttribute("constants", "nu");
	m_dt = xml.getValueAttribute("time", "dt");

	m_Ux = xml.getValueAttribute("flow", "ux");
	m_Uz = xml.getValueAttribute("flow", "uz");

	m_maxGamma = xml.getValueAttribute("constants", "max_Gamma");
	m_maxNumPanelVort = xml.getValueAttribute("constants", "max_NumPanelVort");
	m_cutoff_exp = xml.getValueAttribute("constants", "cutoff_exp");

	m_kernel_threshold = xml.getValueAttribute("constants", "kernel_threshold");

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
	compute_influence_matrix();
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

void VortexSheet::compute_influence_matrix()
{
	//========================================================================
	// Compute influence matrix according to coefficients after Mogenthal
	// =======================================================================

	// Follow Morgenthal (2002)
	unsigned Nl = m_gamma.n_elem;

	std::cout << "Nl = " << Nl << "\n";

	m_infM.set_size(Nl + 1, Nl);

	double c1, c2, c3, c4, c5, c6, c7, c9;
	Matrix p, q;
	p.set_size(Nl, Nl);
	q.set_size(Nl, Nl);

	for (unsigned i = 0; i < Nl; i++) {
		double xci = m_xc(i);
		double zci = m_zc(i);
		double thetai = m_theta(i);

		for (unsigned j = 0; j < Nl; j++) {
			if (i == j) {
				p(i, j) = -1.0;
				q(i, j) = 1.0;
			} else {

				double xj = m_x(j);
				double zj = m_z(j);
				double thetaj = m_theta(j);

				double dsj = m_ds(j);

				c1 = -(xci - xj) * std::cos(thetaj) - (zci - zj) * std::sin(thetaj);
				c2 = std::pow(xci - xj, 2) + std::pow(zci - zj, 2);
				c3 = std::sin(thetai - thetaj);
				c4 = std::cos(thetai - thetaj);
				c5 = (xci - xj) * std::sin(thetaj) - (zci - zj) * std::cos(thetaj);
				c6 = std::log(1.0 + dsj * ((dsj + 2 * c1) / c2));
				c7 = std::atan2((c5 * dsj), (c2 + c1 * dsj));
				c9 = (xci - xj) * std::cos(thetai - 2.0 * thetaj)
				     - (zci - zj) * std::sin(thetai - 2.0 * thetaj);

				q(i, j) =
				    c4 + 0.5 * c9 * c6 / dsj - ((c1 * c3 + c4 * c5) * c7) / dsj;
				p(i, j) = 0.5 * c4 * c6 + c3 * c7 - q(i, j);
			}
		}
	}

	for (unsigned i = 0; i < Nl; i++) {
		m_infM(i, 0) = -m_rpi2 * (p(i, 0) + q(i, Nl - 1));

		for (unsigned j = 1; j < Nl; j++) {
			m_infM(i, j) = -m_rpi2 * (p(i, j) + q(i, j - 1));
		}
	}

	// Enforcing the total circulation
	for (unsigned j = 0; j < Nl; j++) {
		m_infM(Nl, j) = m_ds(j);
	}

	std::cout << "Computed influence matrix of size " << arma::size(m_infM)
	          << " after Kuette and Chow" << std::endl;

}

void VortexSheet::solvevortexsheet(VortexBlobs &blobs)
{
	unsigned Nl = m_gamma.n_elem;

	Vector brhs(Nl + 1);

	if (blobs.size() == 0) {
		brhs.rows(0, Nl - 1) = m_Ux * m_enx + m_Uz * m_enz;
		brhs(Nl) = -blobs.totalcirc();

	} else {
		Vector u(Nl);
		Vector w(Nl);

#pragma omp parallel for
		for (unsigned i = 0; i < Nl; i++) {

			double dK_ij, rsigmasqr, dx_ij, dz_ij, dr_ij2, threshold;

			u(i) = 0.0;
			w(i) = 0.0;

			for (unsigned j = 0; j < blobs.size(); j++) {

				dx_ij = m_xc(i) - blobs.m_x(j);
				dz_ij = m_zc(i) - blobs.m_z(j);
				dr_ij2 = std::pow(dx_ij, 2) + std::pow(dz_ij, 2);

				threshold =
				    m_kernel_threshold * std::pow(blobs.m_sigma(j), 2);
				rsigmasqr = 1.0 / std::pow(blobs.m_sigma(j), 2);

				if (dr_ij2 < threshold) {
					dK_ij = (1.0 - std::exp(-dr_ij2 * rsigmasqr)) / dr_ij2;
				} else {
					dK_ij = 1.0 / dr_ij2;
				}

				u(i) -= dK_ij * dz_ij * blobs.m_circ(j);
				w(i) += dK_ij * dx_ij * blobs.m_circ(j);
			}
		}

		u *= m_rpi2;
		w *= m_rpi2;

		// Not entirely convinced that this is the correct BC (see Morgenthal)
		brhs.rows(0, Nl - 1) = (m_Ux + u) % m_enx + (m_Uz + w) % m_enz;
		brhs(Nl) = -blobs.totalcirc();
	}

	// Solve system
	m_gamma = arma::solve(m_infM.t() * m_infM, m_infM.t() * brhs);
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

VortexBlobs VortexSheet::release_nascent_vortices_rw(Random &_rand)
{
	// double x, z, circ, sigma;
	auto Nvs = size();

	double R, rrw;

	// First we need to determine how many of these vortices will be created at
	// each panel. This is in accordance with Morgenthal PhD 4.3.6 Vortex release
	// algorithm (eq. 4.108)
	Vector_un PanelNewVort(Nvs);       // Panel's new vortices

	Vector PanelCirc = m_gamma % m_ds; // Panel's total circulation

	Vector AbsPanelCirc =
	    arma::abs(PanelCirc / arma::max(PanelCirc)
	              * m_maxNumPanelVort); // Absolute value of the panel vorticity

	Vector Circ_new(Nvs);
	// GD---> Should not need to do a loop here. Unfortunately for some reason
	// armadillo does not allow me to do
	// PanelNewVort=arma::round(AbsPanelCirc)+arma::ones(Nvs,1), perhaps
	// MAB can fix this.
	for (unsigned i = 0; i < Nvs; i++) {
		// Number of released vortices per panel
		PanelNewVort(i) = std::floor(AbsPanelCirc(i)) + 1;
		assert(PanelNewVort(i) > 0);

		// Circulation of each vortex we release from the ith panel
		Circ_new(i) = PanelCirc(i) / PanelNewVort(i);
	}

	// Number of the new vortices.
	// if we set the maximum number of vortices from each panel=1, then Nrv=Nsv
	unsigned Nrv = arma::sum(PanelNewVort);

	// Initialize the vortexblobs ready for release
	VortexBlobs nascentVort(Nrv);

	// Calculating the position of the newelly released vortices
	double xm, zm;
	unsigned counter = 0;

	for (unsigned i = 0; i < Nvs; i++) {
		unsigned j = 0;
		while (j < PanelNewVort(i)) {
			// Find locations at the ith panel from which the vortices should be
			// released
			xm = m_x(i) + m_ds(i) * m_etx(i) / (PanelNewVort(i) + 1);
			zm = m_z(i) + m_ds(i) * m_etz(i) / (PanelNewVort(i) + 1);

			// Here is the tricky bit !!! Morgenthal shows in figure 4.7 that
			// the particles may be released with
			// some random walk in both the normal (to the panel) and the
			// tangential directions. This is
			// wrong according to Chorin 1978. Since we are in the boundary
			// layer, the diffusion process takes place only
			// in the normal direction and not the streamwise. This also
			// according to Prandtl's boundary layer approximation.
			// I will implement the Chorin 1978 here and not the Morgenthal one.
			// In the end, this should not make
			// a huge difference.

			// Create Random walk values
			R = _rand.rand();

			// This is half normal distribution only in the outwards direction
			rrw = std::abs(std::sqrt(4.0 * m_nu * m_dt * std::log(1.0 / R)));

			// Add the released vortex
			nascentVort.m_x(counter) = xm + rrw * m_enx(i);
			nascentVort.m_z(counter) = zm + rrw * m_enz(i);
			nascentVort.m_circ(counter) = Circ_new(i);

			// Now for the cut-off kernel we implement things as suggested by
			// Mirta Perlman 1985 (JCP)
			// using sigma=ds^q where 0.5<q<1. She suggests using q=0.625! This
			// q parameter is the coefficient not the cutoff
			// parameter the inputs file.
			nascentVort.m_sigma(counter) = std::pow(m_ds(i), m_cutoff_exp);

			counter++;
			assert(counter <= Nrv);
			j++;
		}
	}
	return nascentVort;
}

void VortexSheet::reflect(VortexBlobs& vortex)
{
	Vector_un closest_panel(vortex.size());
	Vector min_dist(vortex.size());
	Vector _mirror(2);

	double x_init, z_init, x_0, z_0, x_1, z_1;

	double dx, dz, dr, min_prev;

	for (unsigned i = 0; i < vortex.size(); i++) {
		min_dist(i) = 0;
		closest_panel(i) = 0;

		if (inside_body(vortex.m_x(i), vortex.m_z(i))) {

			// Find which panel is closest to the vortex
			min_prev = 10E6;
			for (unsigned j = 0; j < size(); j++) {

				dx = m_xc(j) - vortex.m_x(i);
				dz = m_zc(j) - vortex.m_z(i);
				dr = std::sqrt(std::pow(dx, 2.0) + std::pow(dz, 2.0));

				if (dr < min_prev) {
					closest_panel(i) = j;
					min_prev = dr;
				}
			}

			// Find the mirror image vortex blob
			x_init = vortex.m_x(i);
			z_init = vortex.m_z(i);

			x_0 = m_x(closest_panel(i));
			z_0 = m_z(closest_panel(i));
			x_1 = m_x(closest_panel(i) + 1);
			z_1 = m_z(closest_panel(i) + 1);

			_mirror = mirror(x_init, z_init, x_0, z_0, x_1, z_1);

			vortex.m_x(i) = _mirror(0);
			vortex.m_z(i) = _mirror(1);
            // All other properties of the vortices 
		}
	}
}

int VortexSheet::inside_body(const double& xcoor, const double& zcoor)
{
    int cn = 0; 
    
    for (unsigned i = 0; i < size(); i++) {
		if (((m_z(i) <= zcoor) && (m_z(i + 1) > zcoor))
		    || ((m_z(i) > zcoor) && (m_z(i + 1) <= zcoor))) {

			float vt = (float)(zcoor - m_z(i)) / (m_z(i + 1) - m_z(i));
			if (xcoor < m_x(i) + vt * (m_x(i + 1) - m_x(i))) {
				++cn;
			}
		}
	}
	return (cn & 1);
}

Vector VortexSheet::mirror(const double &x_init,
                           const double &z_init,
                           const double &x_0,
                           const double &z_0,
                           const double &x_1,
                           const double &z_1)
{
	Vector p2(2);
	double dx, dz, a, b;

	dx = x_1 - x_0;
	dz = z_1 - z_0;

	a = (dx * dx - dz * dz) / (dx * dx + dz * dz);
	b = 2.0 * dx * dz / (dx * dx + dz * dz);

	p2(0) = a * (x_init - x_0) + b * (z_init - z_0) + x_0;
	p2(1) = b * (x_init - x_0) - a * (z_init - z_0) + z_0;

	return p2;
}

void VortexSheet::compute_loads(double Ur)
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
    double pref=0.5*m_rho*Ur*Ur;

    P += (pref-pmax)*arma::ones(size(),1); 

    // To find the forces we need to measure on the particle
    // we need to change sign 
    m_fx = -arma::sum(P % m_enx % m_ds);
	m_fz = -arma::sum(P % m_enz % m_ds);
    
    std::cout<<"C_D = "<<m_fx/(0.5*m_rho*1.0*Ur*Ur)<<"\t"<<"C_L = "<<m_fz/(0.5*m_rho*1.0*Ur*Ur)<<std::endl;
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
