#include "VortexSheet.hpp"

VortexSheet::VortexSheet()
{
	// Don't put anything in here, it is never called
}

VortexSheet::VortexSheet(const XmlHandler &xml, const std::string &stamp)
{
	m_rho = xml.getValueAttribute("constants", "density");
	m_nu = xml.getValueAttribute("constants", "nu");
	m_dt = xml.getValueAttribute("time", "dt");

	m_Ux = xml.getValueAttribute("flow", "ux");
	m_Uz = xml.getValueAttribute("flow", "uz");

	m_maxNumPanelVort = xml.getValueAttribute("constants", "max_NumPanelVort");
	m_cutoff_exp = xml.getValueAttribute("constants", "cutoff_exp");

	auto indir = xml.getStringAttribute("io", "input_dir");
	// Make sure we have trailing separator
	if (indir.back() != '/') {
		indir += '/';
	}

	std::string domain = indir + xml.getStringAttribute("io", "domain_file");

	read_input_coord(domain);
	form_vortex_sheet();
	compute_influence_matrix();

	// Initialise the output file
	auto outdir = xml.getStringAttribute("io", "output_dir", true);
	m_gammafile.open(outdir + stamp + std::string("_gamma.dat"));
	m_gammafile << m_xc.n_elem << " # Number of collocation points\n"
	            << m_dt << " # Time Step\n"
	            << xml.getValueAttribute("time", "steps") << " # Steps\n"
	            << "Time [s] Gamma - Vortex sheet strength\n";

	m_forcefile.open(outdir + stamp + std::string("_loads.dat"));
	m_forcefile << "Time [s]\tF_x [-]\tF_z [s]\n";
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

	// Assume for the moment that we are looking at a cylinder - get the
	// diameter for the drag coefficient reference area.
	m_A = arma::max(m_x) - arma::min(m_x);
	std::cout << "A = " << m_A << " (assuming cylinder)\n";
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

				c1 = -(xci - xj) * std::cos(thetaj)
				     - (zci - zj) * std::sin(thetaj);
				c2 = std::pow(xci - xj, 2) + std::pow(zci - zj, 2);
				c3 = std::sin(thetai - thetaj);
				c4 = std::cos(thetai - thetaj);
				c5 = (xci - xj) * std::sin(thetaj)
				     - (zci - zj) * std::cos(thetaj);
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
		m_infM(i, 0) = -math::rpi2 * (p(i, 0) + q(i, Nl - 1));

		for (unsigned j = 1; j < Nl; j++) {
			m_infM(i, j) = -math::rpi2 * (p(i, j) + q(i, j - 1));
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

				threshold = 10 * std::pow(blobs.m_sigma(j), 2);
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

		u *= math::rpi2;
		w *= math::rpi2;

		// Not entirely convinced that this is the correct BC (see Morgenthal)
		brhs.rows(0, Nl - 1) = (m_Ux + u) % m_enx + (m_Uz + w) % m_enz;
		brhs(Nl) = -blobs.totalcirc();
	}

	// Solve system
	m_gamma = arma::solve(m_infM.t() * m_infM, m_infM.t() * brhs);
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
				blobs.m_uvs(i) +=
				    px(i, last) * m_gamma(last) + qx(i, last) * m_gamma(0);
				blobs.m_wvs(i) +=
				    py(i, last) * m_gamma(last) + qy(i, last) * m_gamma(0);
			} else {
				blobs.m_uvs(i) +=
				    px(i, j) * m_gamma(j) + qx(i, j) * m_gamma(j + 1);
				blobs.m_wvs(i) +=
				    py(i, j) * m_gamma(j) + qy(i, j) * m_gamma(j + 1);
			}
		}
		blobs.m_uvs(i) *= math::rpi2;
		blobs.m_wvs(i) *= math::rpi2;
	}
}

VortexBlobs VortexSheet::release_nascent_vortices_rw(Random &_rand)
{
	// double x, z, circ, sigma;
	auto Nvs = size();

	double R, rrw;

	// First we need to determine how many of these vortices will be created at
	// each panel. This is in accordance with Morgenthal PhD 4.3.6 Vortex
	// release
	// algorithm (eq. 4.108)
	Vector PanelCirc = m_gamma % m_ds; // Panel's total circulation

	// Absolute value of the panel vorticity
	Vector AbsPanelCirc =
	    arma::abs(PanelCirc / arma::max(PanelCirc) * m_maxNumPanelVort);

	Vector PanelNewVort = arma::floor(AbsPanelCirc) + arma::ones(Nvs);
	// Make sure the new panel vortices are non-negative
	if (!arma::all(PanelNewVort > 1.0e-15)) {
		throw std::string("VortexSheet::release_nascent_vortices_rw -> "
		                  "negative PanelNewVort");
	}

	Vector Circ_new = PanelCirc / PanelNewVort;

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

void VortexSheet::write_step(double time)
{
	for (unsigned i = 0; i < m_xc.n_elem; ++i) {
		m_gammafile << time << "\t " << m_xc(i) << "\t" << m_zc(i) << "\t"
		            << m_gamma(i) << "\t" << m_ds(i) << "\n";
	}

	m_forcefile << time << "\t" << m_fx << "\t" << m_fz << "\n";
}

void VortexSheet::reflect(VortexBlobs &vortex)
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
			// All other properties of the vortices remain the same
		}
	}
}

bool VortexSheet::inside_body(double xcoor, double zcoor)
{
	// Using this algorithm, cn == 1 when the point is inside the boundary
	int cn = 0;

	for (unsigned i = 0; i < size(); i++) {
		if (((m_z(i) <= zcoor) && (m_z(i + 1) > zcoor))
		    || ((m_z(i) > zcoor) && (m_z(i + 1) <= zcoor))) {

			auto vt = (zcoor - m_z(i)) / (m_z(i + 1) - m_z(i));
			if (xcoor < m_x(i) + vt * (m_x(i + 1) - m_x(i))) {
				++cn;
			}
		}
	}

	return (cn == 1);
}

Vector VortexSheet::mirror(double x_init,
                           double z_init,
                           double x_0,
                           double z_0,
                           double x_1,
                           double z_1)
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
	Vector P(size(), arma::fill::zeros);
	Vector Dp = -m_rho / m_dt * (m_gamma % m_ds);

	for (unsigned i = 0; i < size() - 1; i++) {
		P(i + 1) = P(i) + Dp(i);
	}

	// Finding the average of the two pressures at the boundary
	for (unsigned i = 0; i < size() - 1; i++) {
		P(i) = 0.5 * (P(i) + P(i + 1));
	}

	P(size() - 1) = 0.5 * (P(size() - 1) + P(0));
	double pmax = arma::max(P);
	double pref = 0.5 * m_rho * Ur * Ur;

	P += (pref - pmax) * arma::ones(size(), 1);

	// To find the forces we need to measure on the particle
	// we need to change sign
	m_fx = -arma::sum(P % m_enx % m_ds);
	m_fz = -arma::sum(P % m_enz % m_ds);

	std::cout << "C_D = " << m_fx / (0.5 * m_rho * m_A * Ur * Ur) << "\t"
	          << "C_L = " << m_fz / (0.5 * m_rho * m_A * Ur * Ur) << std::endl;
}

unsigned VortexSheet::size()
{
	if ((m_gamma.size() != m_xc.size()) && (m_gamma.size() != m_zc.size())
	    && (m_gamma.size() != m_x.size() - 1)
	    && (m_gamma.size() != m_z.size() - 1)
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
		std::cout << "(xc,zc) = (" << m_xc(i) << "," << m_zc(i) << ")"
		          << std::endl;
	}
}

void VortexSheet::print_unit_vectors()
{
	for (unsigned i = 0; i < size(); i++) {
		std::cout << "For panel No " << i << "the normal vectors are en = ("
		          << m_enx(i) << "," << m_enz(i) << ") and et = (" << m_etx(i)
		          << "," << m_etz(i) << ")" << std::endl;
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

void VortexSheet::probe_velocities(const VortexBlobs &blobs, Probe &probe)
{
	double rsigmasqr;
	double dx_ij, dz_ij, dr_ij2, threshold, xkernel;
	double dK_ij, zkernel;
	double x_i, z_i, u_i, w_i;

	// Compute the velocity vector at the probe points

	for (unsigned i = 0; i < probe.size(); i++) {
		probe.m_u(i) = 0.0;
		probe.m_w(i) = 0.0;
	}

	for (unsigned i = 0; i < probe.size(); i++) {
		x_i = probe.m_x(i);
		z_i = probe.m_z(i);
		u_i = 0.0;
		w_i = 0.0;
		for (unsigned j = 1; j < blobs.size(); j++) {
			dx_ij = x_i - blobs.m_x(j);
			dz_ij = z_i - blobs.m_z(j);
			dr_ij2 = std::pow(dx_ij, 2) + std::pow(dz_ij, 2);

			threshold = 10 * std::pow(blobs.m_sigma(j), 2);
			rsigmasqr = 1.0 / std::pow(blobs.m_sigma(j), 2);

			if (dr_ij2 < threshold) {
				dK_ij = (1.0 - std::exp(-dr_ij2 * rsigmasqr)) / dr_ij2;
			} else {
				dK_ij = 1.0 / dr_ij2;
			}

			xkernel = dK_ij * dz_ij;
			zkernel = dK_ij * dx_ij;
			u_i = -xkernel * blobs.m_circ(j);
			w_i = +zkernel * blobs.m_circ(j);
			probe.m_u(i) = probe.m_u(i) + u_i;
			probe.m_w(i) = probe.m_w(i) + w_i;
		}
	}

	double c1, c2, c5, c6, c7, c8, c9;
	Matrix px, py, qx, qy;
	px.set_size(probe.size(), size());
	qx.set_size(probe.size(), size());
	py.set_size(probe.size(), size());
	qy.set_size(probe.size(), size());

	for (unsigned i = 0; i < probe.size(); i++) {

		probe.m_uvs(i) = 0.0;
		probe.m_wvs(i) = 0.0;
	}

	for (unsigned i = 0; i < probe.size(); i++) {
		for (unsigned j = 1; j < size(); j++) {
			c1 = -(probe.m_x(i) - m_xc(j)) * cos(m_theta(j))
			     - (probe.m_z(i) - m_zc(j)) * sin(m_theta(j));
			c2 = std::pow((probe.m_x(i) - m_xc(j)), 2)
			     + std::pow((probe.m_z(i) - m_zc(j)), 2);
			c5 = (probe.m_x(i) - m_xc(j)) * sin(m_theta(j))
			     - (probe.m_z(i) - m_zc(j)) * cos(m_theta(j));
			c6 = log(1.0 + m_ds(j) * ((m_ds(j) + 2 * c1) / c2));
			c7 = atan2((c5 * m_ds(j)), (c2 + c1 * m_ds(j)));
			c8 = (probe.m_x(i) - m_xc(j)) * sin(-2.0 * m_theta(j))
			     + (probe.m_z(i) - m_zc(j)) * cos(-2.0 * m_theta(j));
			c9 = (probe.m_x(i) - m_xc(j)) * cos(-2.0 * m_theta(j))
			     + (probe.m_z(i) - m_zc(j)) * sin(-2.0 * m_theta(j));

			qx(i, j) =
			    -sin(m_theta(j)) + 0.5 * c8 * c6 / m_ds(j)
			    + (c1 * cos(m_theta(j)) + c5 * sin(m_theta(j))) * c7 / m_ds(j);
			px(i, j) =
			    -0.5 * c6 * sin(m_theta(j)) - c7 * cos(m_theta(j)) - qx(i, j);

			qy(i, j) =
			    cos(m_theta(j)) + 0.5 * c9 * c6 / m_ds(j)
			    + (c1 * cos(m_theta(j)) - c5 * sin(m_theta(j))) * c7 / m_ds(j);
			py(i, j) =
			    -0.5 * c6 * cos(m_theta(j)) - c7 * sin(m_theta(j)) - qy(i, j);

			if (j == size() - 1) {

				probe.m_uvs(i) = (px(i, size() - 1) * m_gamma(size() - 1)
				                  + qx(i, size() - 1) * m_gamma(0));
				probe.m_wvs(i) = (py(i, size() - 1) * m_gamma(size() - 1)
				                  - qy(i, size() - 1) * m_gamma(0));

			} else {
				probe.m_uvs(i) =
				    (px(i, j) * m_gamma(j) + qx(i, j) * m_gamma(j + 1));
				probe.m_wvs(i) =
				    (py(i, j) * m_gamma(j) + qy(i, j) * m_gamma(j + 1));
			}
		}
	}

	for (unsigned i = 0; i < probe.size(); i++) {
		probe.m_u(i) = math::rpi2 * probe.m_u(i) + probe.m_uvs(i) + m_Ux;
		probe.m_w(i) = math::rpi2 * probe.m_w(i) + probe.m_wvs(i) + m_Uz;
	}
}
