#include "DVMBase.hpp"
#include <armadillo>
#include <cmath>
#include <iomanip>
#include <time.h>

DVMBase::DVMBase()
{
	m_pi = 4.0 * atan(1.0);
	m_step = 0;
	m_rpi2 = 1.0 / (2.0 * m_pi);
}

void DVMBase::init(pugi::xml_document &xml_doc)
{
	std::string in_dir = xml_doc.child("IC2DDVM")
	                         .child("io")
	                         .child("input_dir")
	                         .attribute("string")
	                         .value();
	std::string out_dir = xml_doc.child("IC2DDVM")
	                          .child("io")
	                          .child("output_dir")
	                          .attribute("string")
	                          .value();
	std::string domain_file = xml_doc.child("IC2DDVM")
	                              .child("io")
	                              .child("domain_file")
	                              .attribute("string")
	                              .value();

	double rho = atof(xml_doc.child("IC2DDVM")
	                      .child("constants")
	                      .child("density")
	                      .attribute("val")
	                      .value());
	if (rho <= 0) {
		throw std::string(
		    "<IC2DDVM><constants><density> must be larger than zero");
	}

	double nu = atof(xml_doc.child("IC2DDVM")
	                     .child("constants")
	                     .child("nu")
	                     .attribute("val")
	                     .value());
	if (nu <= 0) {
		throw std::string("<IC2DDVM><constants><nu> must be larger than zero");
	}

	double max_gamma = atof(xml_doc.child("IC2DDVM")
	                            .child("constants")
	                            .child("max_gamma")
	                            .attribute("val")
	                            .value());
	if (max_gamma <= 0) {
		throw std::string(
		    "<IC2DDVM><constants><max_gamma> must be larger than zero");
	}

	double kernel = atof(xml_doc.child("IC2DDVM")
	                         .child("constants")
	                         .child("kernel_threshold")
	                         .attribute("val")
	                         .value());
	if (kernel <= 0) {
		throw std::string(
		    "<IC2DDVM><constants><kernel_threshold> must be larger than zero");
	}

	double cutoff = atof(xml_doc.child("IC2DDVM")
	                         .child("constants")
	                         .child("sigma_cutoff")
	                         .attribute("val")
	                         .value());
	if (cutoff <= 0) {
		throw std::string(
		    "<IC2DDVM><constants><sigma_cutoff> must be larger than zero");
	}
	if (cutoff >= 1) {
		throw std::string(
		    "<IC2DDVM><constants><sigma_cutoff> must be less than one");
	}

	double dt = atof(xml_doc.child("IC2DDVM")
	                     .child("time")
	                     .child("dt")
	                     .attribute("val")
	                     .value());
	if (dt <= 0) {
		throw std::string("<IC2DDVM><time><dt> must be larger than zero");
	}

	unsigned steps = atof(xml_doc.child("IC2DDVM")
	                          .child("time")
	                          .child("steps")
	                          .attribute("val")
	                          .value());
	if (steps <= 0) {
		throw std::string("<IC2DDVM><time><steps> must be larger than zero");
	}

	double ux = atof(xml_doc.child("IC2DDVM")
	                     .child("flow")
	                     .child("ux")
	                     .attribute("val")
	                     .value());
	double uz = atof(xml_doc.child("IC2DDVM")
	                     .child("flow")
	                     .child("uz")
	                     .attribute("val")
	                     .value());

	// Generalise this to arrays
	unsigned probe_x = atof(xml_doc.child("IC2DDVM")
	                            .child("probe")
	                            .child("x")
	                            .attribute("val")
	                            .value());
	unsigned probe_z = atof(xml_doc.child("IC2DDVM")
	                            .child("probe")
	                            .child("z")
	                            .attribute("val")
	                            .value());

	// Directories and file names
	m_in_dir = in_dir;
	m_out_dir = out_dir;
	m_domain_file = domain_file;

	// Parameters
	m_maxGamma = max_gamma;
	m_kernel_threshold = kernel;
	m_sigma_cutoff = cutoff;

	m_dt = dt;
	m_rho = rho;
	m_steps = steps;

	m_nu = nu;
	m_Ux = ux;
	m_Uz = uz;

	// Generalise this to arrays
	m_probe.resize(1);

	// Define probe point
	m_probe.x[0] = probe_x;
	m_probe.z[0] = probe_z;
}

void DVMBase::read_input_coord()
{
	// Read the body coordinates
	std::string file = m_in_dir + m_domain_file;
	std::string line;
	std::ifstream coor_file(file.c_str());
	int index = 0;
	double tmp1, tmp2;

	if (coor_file.is_open()) {
		while (coor_file.good()) {
			std::getline(coor_file, line);
			std::istringstream buffer(line);
			buffer >> tmp1 >> tmp2;
			m_body.x.push_back(tmp1);
			m_body.z.push_back(tmp2);
		}
		m_body.x.pop_back();
		m_body.z.pop_back();
		std::cout << "Succesfully loaded coordinate file with " << m_body.size()
		          << " points." << std::endl;
		coor_file.close();
		m_n = m_body.size();
	} else {
		std::string error_msg;
		error_msg = "Unable to open coordinate file from " + file;
		throw error_msg;
	}

	m_body.print_location();
}

void DVMBase::init_outputs()
{

	// Init all io files
	time_t rawtime;
	struct tm *ptm;
	time(&rawtime);
	ptm = gmtime(&rawtime);

	std::ostringstream os;
	os << ptm->tm_year + 1900 << "_";
	if (ptm->tm_mon + 1 < 10) {
		os << "0";
	}
	os << ptm->tm_mon + 1 << "_";
	if (ptm->tm_mday < 10) {
		os << "0";
	}
	os << ptm->tm_mday << "_";
	if (ptm->tm_hour + 1 < 10) {
		os << "0";
	}
	os << ptm->tm_hour + 1 << "_";
	if (ptm->tm_min < 10) {
		os << "0";
	}
	os << ptm->tm_min << "_";
	if (ptm->tm_sec < 10) {
		os << "0";
	}
	os << ptm->tm_sec;

	std::cout << "File timestamp is " << os.str() << std::endl;

	dev_dvm.open((m_out_dir + os.str() + std::string("_vortex.dat")).c_str());
	dev_dvm << m_vortex.size() << " # Number of Nodes" << std::endl;
	dev_dvm << m_dt << " # Time Step" << std::endl;
	dev_dvm << m_steps << " # Steps " << std::endl;
	dev_dvm << "Time[s]"
	        << " "
	        << " x-position [m]"
	        << " "
	        << "z-position [m]"
	        << " "
	        << "circulation" << std::endl;

	dev_Num.open(
	    (m_out_dir + os.str() + std::string("_vortex_num.dat")).c_str());
	dev_Num << "Time [s]"
	        << " "
	        << "Number of vortices" << std::endl;

	dev_gamma.open((m_out_dir + os.str() + std::string("_gamma.dat")).c_str());
	dev_gamma << m_vortsheet.size() << " # Number of collocation points"
	          << std::endl;
	dev_gamma << m_dt << " # Time Step" << std::endl;
	dev_gamma << m_steps << " # Steps " << std::endl;
	dev_gamma << "Time [s]"
	          << " "
	          << "Gamma - Vortex sheet strength" << std::endl;

	dev_loads.open((m_out_dir + os.str() + std::string("_loads.dat")).c_str());
	dev_loads << "Time [s]"
	          << " "
	          << "\tFx [N]"
	          << "\tFz [N]" << std::endl;

	dev_probe.open((m_out_dir + os.str() + std::string("_probe.dat")).c_str());
	dev_probe << m_dt << " # Time Step" << std::endl;
	dev_probe << m_steps << " # Steps " << std::endl;
	dev_probe << "Time [s]"
	          << "\t"
	          << "u [m/s]"
	          << "\t"
	          << "w{m/s" << std::endl;
}

void DVMBase::form_vortex_sheet()
{
	// Initialize vortex sheet
	m_vortsheet.resize(m_n - 1);
	double dx, dz, theta;

	// Define the characteristics of the vortex sheet
	for (unsigned i = 0; i < m_vortsheet.size(); i++) {

		// collocation points on the vortex sheet
		m_vortsheet.xc[i] = 0.5 * (m_body.x[i] + m_body.x[i + 1]);
		m_vortsheet.zc[i] = 0.5 * (m_body.z[i] + m_body.z[i + 1]);

		// ds and theta along the vortex sheet
		dx = m_body.x[i + 1] - m_body.x[i];
		dz = m_body.z[i + 1] - m_body.z[i];
		theta = atan2(dz, dx);

		m_vortsheet.ds[i] = std::sqrt(std::pow(dx, 2) + std::pow(dz, 2));
		m_vortsheet.theta[i] = theta;

		// normals and tangentials

		// Inwards facing
		// m_vortsheet.etx[i] = cos(theta);
		// m_vortsheet.etz[i] = sin(theta);
		// m_vortsheet.enx[i] = -m_vortsheet.etz[i];
		// m_vortsheet.enz[i] = m_vortsheet.etx[i];

		// Outwards facing
		m_vortsheet.enx[i] = dz / m_vortsheet.ds[i];
		m_vortsheet.enz[i] = -dx / m_vortsheet.ds[i];
		m_vortsheet.etx[i] = -m_vortsheet.enz[i];
		m_vortsheet.etz[i] = m_vortsheet.enx[i];
	}

	// m_Gamma_abs.resize(m_vortsheet.size());

	m_vortsheet.print_collocation();
	m_vortsheet.print_unit_vectors();

	std::cout << "Created vortex sheet of size " << m_vortsheet.size()
	          << std::endl;

	// throw std::string("Stop here");
}

void DVMBase::compute_influence_matrix()
{

	//========================================================================
	// Compute influence matrix according to coefficients after Mogenthal
	// =======================================================================

	std::cout << "m_n = " << m_n << std::endl;
	std::cout << "m_body.size() = " << m_body.size() << std::endl;
	std::cout << "m_vortsheet.size() = " << m_vortsheet.size() << std::endl;

	// Follow Morgenthal (2002)

	unsigned Nl = m_vortsheet.size();

	m_infM.resize(Nl + 1, std::vector<double>(Nl));

	std::cout << "m_infM.size(1) " << m_infM.size() << std::endl;
	std::cout << "m_infM.size(2) " << m_infM[0].size() << std::endl;

	// throw std::string("Stop here");

	double c1, c2, c3, c4, c5, c6, c7, c8, c9;
	Matrix p, q;
	p.resize(Nl, std::vector<double>(Nl));
	q.resize(Nl, std::vector<double>(Nl));

	for (unsigned i = 0; i < Nl; i++) {

		double xci = m_vortsheet.xc[i];
		double zci = m_vortsheet.zc[i];
		double thetai = m_vortsheet.theta[i];

		for (unsigned j = 0; j < Nl; j++) {
			if (i == j) {
				p[i][j] = -1.0;
				q[i][j] = 1.0;
			} else {

				double xj = m_body.x[j];
				double zj = m_body.z[j];
				double thetaj = m_vortsheet.theta[j];

				double dsj = m_vortsheet.ds[j];

				c1 = -(xci - xj) * cos(thetaj) - (zci - zj) * sin(thetaj);
				c2 = std::pow(xci - xj, 2) + std::pow(zci - zj, 2);
				c3 = sin(thetai - thetaj);
				c4 = cos(thetai - thetaj);
				c5 = (xci - xj) * sin(thetaj) - (zci - zj) * cos(thetaj);
				c6 = log(1.0 + dsj * ((dsj + 2 * c1) / c2));
				c7 = atan2((c5 * dsj), (c2 + c1 * dsj));
				c8 = (xci - xj) * sin(thetai - 2.0 * thetaj)
				     + (zci - zj) * cos(thetai - 2.0 * thetaj);
				c9 = (xci - xj) * cos(thetai - 2.0 * thetaj)
				     - (zci - zj) * sin(thetai - 2.0 * thetaj);

				q[i][j] =
				    c4 + 0.5 * c9 * c6 / dsj - ((c1 * c3 + c4 * c5) * c7) / dsj;
				p[i][j] = 0.5 * c4 * c6 + c3 * c7 - q[i][j];
			}
		}
	}

	for (unsigned i = 0; i < Nl; i++) {

		m_infM[i][0] = -m_rpi2 * (p[i][0] + q[i][Nl - 1]);

		for (unsigned j = 1; j < Nl; j++) {
			m_infM[i][j] = -m_rpi2 * (p[i][j] + q[i][j - 1]);
		}
	}

	// Enforcing the total circulation
	for (unsigned j = 0; j < Nl; j++) {
		m_infM[Nl][j] = m_vortsheet.ds[j];
	}

	// Print the matrix

	//===============================================
	for (unsigned i = 0; i < m_infM.size(); i++) {
		for (unsigned j = 0; j < m_infM[0].size(); j++) {
			std::cout << std::setprecision(3) << m_infM[i][j] << '\t';
		}
		std::cout << "\n";
	}
	//==============================================

	// throw std::string("Stop here");

	std::cout << "Computed influence matrix of size " << m_infM.size() << ","
	          << m_infM[0].size() << " after Kuette and Chow" << std::endl;
}

double DVMBase::get_time()
{
	return m_time;
}

unsigned DVMBase::get_size()
{
	return m_vortex.size();
}

unsigned DVMBase::get_vs_size()
{
	return m_vortsheet.size();
}

unsigned DVMBase::get_steps()
{
	return m_steps;
}

void DVMBase::increment_step()
{
	m_step++;
	m_time = m_step * m_dt;
}

void DVMBase::write_outputs()
{
	for (unsigned i = 0; i < m_vortex.size(); i++) {
		dev_dvm << m_time << " " << m_vortex.x[i] << " " << m_vortex.z[i] << " "
		        << m_vortex.circ[i] << std::endl;
	}

	dev_Num << m_step << " " << m_vortex.size() << std::endl;

	for (unsigned i = 0; i < m_vortsheet.size(); i++) {
		dev_gamma << m_time << "\t " << m_vortsheet.xc[i] << "\t "
		          << m_vortsheet.zc[i] << " " << m_vortsheet.gamma[i] << "\t"
		          << m_vortsheet.ds[i] << std::endl;
	}

	dev_loads << m_step << " " << m_fx << "\t" << m_fz << std::endl;

	for (unsigned i = 0; i < m_probe.size(); i++) {
		dev_probe << m_time << " " << m_probe.u[i] << " " << m_probe.w[i]
		          << std::endl;
	}
}

void DVMBase::biotsavart()
{

#pragma omp parallel for
	for (unsigned i = 0; i < m_vortex.size(); i++) {

		double dx_ij, dz_ij, dK_ij, dr_ij2, threshold, rsigmasqr;

		m_vortex.u[i] = 0.0;
		m_vortex.w[i] = 0.0;

		for (unsigned j = 0; j < m_vortex.size(); j++) {

			if (i != j) {
				dx_ij = m_vortex.x[i] - m_vortex.x[j];
				dz_ij = m_vortex.z[i] - m_vortex.z[j];
				dr_ij2 = std::pow(dx_ij, 2) + std::pow(dz_ij, 2.0);

				threshold =
				    m_kernel_threshold * std::pow(m_vortex.sigma[j], 2.0);
				rsigmasqr = 1.0 / std::pow(m_vortex.sigma[j], 2.0);

				if (dr_ij2 < threshold) {
					dK_ij = (1.0 - std::exp(-dr_ij2 * rsigmasqr)) / dr_ij2;
				} else {
					dK_ij = 1.0 / dr_ij2;
				}

				m_vortex.u[i] -= dK_ij * dz_ij * m_vortex.circ[j];
				m_vortex.w[i] += dK_ij * dx_ij * m_vortex.circ[j];
			}
		}
		m_vortex.u[i] *= m_rpi2;
		m_vortex.w[i] *= m_rpi2;
	}
}

void DVMBase::convect(unsigned order)
{
	// Convects according to an order of magnitude
	if (order == 1) {
		biotsavart();
		vortexsheetbc();

		for (unsigned i = 0; i < m_vortex.size(); i++) {
			m_vortex.x[i] += (m_vortex.u[i] + m_vortex.uvs[i] + m_Ux) * m_dt;
			m_vortex.z[i] += (m_vortex.w[i] + m_vortex.wvs[i] + m_Uz) * m_dt;

			// std::cout << "Vortex [" << i << "] x = " << m_vortex.x[i] << " z
			// = " << m_vortex.z[i] << std::endl;

			// std::cout << "Ux = " << m_Ux << " Uz = " << m_Uz << std::endl;
		}
	} else if (order == 2) {

		throw std::string("DVMBase::convet - only first order convection "
		                  "supported at present");

		// This requires the full problem to be solved at a higher level -
		// probably best done in main.cpp?
		/* std::vector<double> x_old, z_old,u_old, w_old, u_tmp, w_tmp;
		 x_old.resize(m_vortex.size());
		 z_old.resize(m_vortex.size());
		 u_old.resize(m_vortex.size());
		 w_old.resize(m_vortex.size());
		 u_tmp.resize(m_vortex.size());
		 w_tmp.resize(m_vortex.size());

		 biotsavart();
		 vortexsheetbc();
		 for(unsigned i=0;i<m_vortex.size();i++)
		 {
		     x_old[i]=m_vortex.x[i]; // Old values k step//
		     z_old[i]=m_vortex.z[i];//  Old values k step//
		     u_old[i]=m_vortex.u[i]+m_vortex.uvs[i]+m_Ux;
		     w_old[i]=m_vortex.w[i]+m_vortex.wvs[i]+m_Uz;
		     m_vortex.x[i]=x_old[i]+2.0/3.0*u_old[i]*dt; // Temporary values
		     m_vortex.z[i]=z_old[i]+2.0/3.0*w_old[i]*dt; // Temporary values

		 }

		 biotsavart();
		 vortexsheetbc();
		 for(unsigned i=0;i<m_vortex.size();i++)
		 {
		     u_tmp[i]=m_vortex.u[i]+m_vortex.uvs[i]+m_Ux;
		     w_tmp[i]=m_vortex.w[i]+m_vortex.wvs[i]+m_Uz;
		     m_vortex.x[i]=x_old[i] + dt*(3.0/4.0*u_old[i]+1.0/4.0*u_tmp[i]);
		     m_vortex.z[i]=z_old[i] + dt*(3.0/4.0*w_old[i]+1.0/4.0*w_tmp[i]);
		 }*/
	} else {
		throw std::string("DVMBase::convect - only first and second order "
		                  "convection supported");
	}
}

void DVMBase::solvevortexsheet()
{

	unsigned Nl = m_vortsheet.size();

	std::vector<double> brhs;
	brhs.resize(Nl + 1);

	if (m_vortex.size() == 0) {

		for (unsigned i = 0; i < Nl; i++) {
			brhs[i] = (m_Ux * m_vortsheet.enx[i] + m_Uz * m_vortsheet.enz[i]);
		}
		brhs[Nl] = -m_vortex.totalcirc();
		;

	} else {
		std::vector<double> u, w;
		u.resize(m_vortsheet.size());
		w.resize(m_vortsheet.size());

#pragma omp parallel for
		for (unsigned i = 0; i < m_vortsheet.size(); i++) {

			double dK_ij, rsigmasqr, dx_ij, dz_ij, dr_ij2, threshold;

			u[i] = 0.0;
			w[i] = 0.0;

			for (unsigned j = 0; j < m_vortex.size(); j++) {

				dx_ij = m_vortsheet.xc[i] - m_vortex.x[j];
				dz_ij = m_vortsheet.zc[i] - m_vortex.z[j];
				dr_ij2 = std::pow(dx_ij, 2) + std::pow(dz_ij, 2);

				threshold = m_kernel_threshold * std::pow(m_vortex.sigma[j], 2);
				rsigmasqr = 1.0 / std::pow(m_vortex.sigma[j], 2);

				if (dr_ij2 < threshold) {
					dK_ij = (1.0 - std::exp(-dr_ij2 * rsigmasqr)) / dr_ij2;
				} else {
					dK_ij = 1.0 / dr_ij2;
				}

				u[i] -= dK_ij * dz_ij * m_vortex.circ[j];
				w[i] += dK_ij * dx_ij * m_vortex.circ[j];
			}
		}

		for (unsigned i = 0; i < Nl; i++) {
			u[i] *= m_rpi2;
			w[i] *= m_rpi2;
		}

		// Not entirely convinced that this is the correct BC (see Morgenthal)
		for (unsigned i = 0; i < Nl; i++) {
			brhs[i] = ((m_Ux + u[i]) * m_vortsheet.enx[i]
			           + (m_Uz + w[i]) * m_vortsheet.enz[i]);
		}
		brhs[Nl] = -m_vortex.totalcirc();
	}

	// Use Armadillo to solve the overdetermined system
	arma::mat MM, B, X;

	MM.zeros(m_infM.size(), m_infM[0].size());
	B.zeros(brhs.size(), 1);

	// Copy Files into Armadillo format
	// for(unsigned i=0; i<m_n+1; i++){ //JS edit
	for (unsigned i = 0; i < B.n_rows; i++) {
		B(i, 0) = brhs[i];
	}

	for (unsigned i = 0; i < MM.n_rows; i++) {
		for (unsigned j = 0; j < MM.n_cols; j++) {
			MM(i, j) = m_infM[i][j];
		}
	}

	// Solve system
	X = arma::solve(MM.t() * MM, MM.t() * B);

	// Copy back to the vortexsheet - removed special treatment here
	// Why was this included? gamma_0 simply solution of Matrix
	for (unsigned i = 0; i < m_vortsheet.size(); i++) {
		m_vortsheet.gamma[i] = X(i, 0);
	}

	/*

	// Print the RHS
	for (unsigned i=0; i < B.n_rows; i++){
	    std::cout << "RHS " << i << "\t" << B(i,0) << std::endl;
	}
	std::cout << std::endl;

	// Print the matrix
	//===============================================
	for (unsigned i=0; i< MM.n_rows; i++){
	    for(unsigned j=0; j< MM.n_cols; j++){
	        std::cout<< std::setprecision(3) << MM(i,j)<<'\t';
	    }
	    std::cout << std::endl;
	}
	std::cout << std::endl;
	//==============================================

	// Test that the output is correct and satisfies Ax = b
	std::vector<double> brhs_check;
	brhs_check.resize(Nl+1);
	for (unsigned i=0; i<Nl+1; i++){
	    brhs_check[i] = 0;
	    for (unsigned j=0; j<Nl; j++){
	        brhs_check[i] += m_infM[i][j]*X(j);
	    }
	}



	// Print the solution
	for (unsigned i=0; i < X.n_rows; i++){
	    std::cout << "X = " << i << "\t" << X(i,0) << std::endl;
	}

	for (unsigned i=0; i < brhs_check.size(); i++){
	    std::cout << "RHS check = " << i << "\t" << brhs_check[i] << std::endl;
	}

	// Print the circulation
	for (unsigned i=0; i < m_vortsheet.size(); i++){
	    std::cout << "gamma = " << i << "\t" << m_vortsheet.gamma[i] <<
	std::endl;
	}
	std::cout << std::endl;

	throw std::string("Stop here");
	*/
}

void DVMBase::vortexsheetbc()
{
	double c1, c2, c3, c4, c5, c6, c7, c8, c9;
	Matrix px, py, qx, qy;

	px.resize(m_vortex.size(), std::vector<double>(m_vortsheet.size()));
	qx.resize(m_vortex.size(), std::vector<double>(m_vortsheet.size()));
	py.resize(m_vortex.size(), std::vector<double>(m_vortsheet.size()));
	qy.resize(m_vortex.size(), std::vector<double>(m_vortsheet.size()));

// compute the coefficients
#pragma omp parallel for
	for (unsigned i = 0; i < m_vortex.size(); i++) {

		double xi = m_vortex.x[i];
		double zi = m_vortex.z[i];

		for (unsigned j = 0; j < m_vortsheet.size(); j++) {

			double xj = m_vortsheet.xc[j];
			double zj = m_vortsheet.zc[j];
			double thetaj = m_vortsheet.theta[j];
			double dsj = m_vortsheet.ds[j];

			c1 = -(xi - xj) * cos(thetaj) - (zi - zj) * sin(thetaj);
			c2 = std::pow(xi - xj, 2.0) + std::pow(zi - zj, 2.0);
			c5 = (xi - xj) * sin(thetaj) - (zi - zj) * cos(thetaj);
			c6 = log(1.0 + dsj * ((dsj + 2 * c1) / c2));
			c7 = atan2(c5 * dsj, c2 + c1 * dsj);
			c8 =
			    (xi - xj) * sin(-2.0 * thetaj) + (zi - zj) * cos(-2.0 * thetaj);
			c9 =
			    (xi - xj) * cos(-2.0 * thetaj) - (zi - zj) * sin(-2.0 * thetaj);

			qx[i][j] = -sin(thetaj) + 0.5 * c8 * c6 / dsj
			           + (c1 * cos(thetaj) + c5 * sin(thetaj)) * c7 / dsj;
			px[i][j] = -0.5 * c6 * sin(thetaj) - c7 * cos(thetaj) - qx[i][j];

			qy[i][j] = cos(thetaj) + 0.5 * c9 * c6 / dsj
			           + (c1 * sin(thetaj) - c5 * cos(thetaj)) * c7 / dsj;
			py[i][j] = 0.5 * c6 * cos(thetaj) - c7 * sin(thetaj) - qy[i][j];
		}
	}

	// calculate the vortexsheet induced velocities
	unsigned last = m_vortsheet.size() - 1;

#pragma omp parallel for
	for (unsigned i = 0; i < m_vortex.size(); i++) {

		m_vortex.uvs[i] = 0.0;
		m_vortex.wvs[i] = 0.0;

		for (unsigned j = 0; j < m_vortsheet.size(); j++) {
			if (j == last) {
				m_vortex.uvs[i] += px[i][last] * m_vortsheet.gamma[last]
				                   + qx[i][last] * m_vortsheet.gamma[0];
				m_vortex.wvs[i] += py[i][last] * m_vortsheet.gamma[last]
				                   + qy[i][last] * m_vortsheet.gamma[0];
			} else {
				m_vortex.uvs[i] += px[i][j] * m_vortsheet.gamma[j]
				                   + qx[i][j] * m_vortsheet.gamma[j + 1];
				m_vortex.wvs[i] += py[i][j] * m_vortsheet.gamma[j]
				                   + qy[i][j] * m_vortsheet.gamma[j + 1];
			}
		}
		m_vortex.uvs[i] *= m_rpi2;
		m_vortex.wvs[i] *= m_rpi2;
	}
}

void DVMBase::diffrw()
{
	double R1, R2, rrw, thetarw;

	unsigned short seed = time(NULL) * 1000000;
	seed48(&seed);

	for (unsigned i = 0; i < m_vortex.size(); i++) {

		// Generate two random numbers in the range 0...1
		R1 = drand48();
		R2 = drand48();

		// Calculate r and theta for the random walk
		rrw = std::sqrt(4.0 * m_nu * m_dt * std::log(1.0 / R1));
		thetarw = 2.0 * m_pi * R2;

		m_vortex.x[i] += rrw * cos(thetarw);
		m_vortex.z[i] += rrw * sin(thetarw);
	}
}

void DVMBase::diffuse_vs_rw()
{
	double x, z, circ, sigma;

	double R1, R2, rrw, thetarw;

	unsigned short seed = time(NULL) * 1000000;
	seed48(&seed);

	unsigned nv_panel;
	double nv_gamma;

	for (unsigned i = 0; i < m_vortsheet.size(); i++) {

		R1 = drand48();
		R2 = drand48();

		rrw = std::sqrt(4.0 * m_nu * m_dt * std::log(1.0 / R1));
		thetarw = 2.0 * m_pi * R2;

		// The position of the released vortex blox after diffusion using random
		// walk
		x = m_vortsheet.xc[i] + rrw * cos(thetarw);
		z = m_vortsheet.zc[i] + rrw * sin(thetarw);

		// Obtain the circulation (see Morgenthal p. 36)
		if (i == m_vortsheet.size() - 1) {
			circ = 0.5 * (m_vortsheet.gamma[i] * m_vortsheet.ds[i]
			              + m_vortsheet.gamma[0] * m_vortsheet.ds[0]);
		} else {
			circ = 0.5 * (m_vortsheet.gamma[i] * m_vortsheet.ds[i]
			              + m_vortsheet.gamma[i + 1] * m_vortsheet.ds[i + 1]);
		}
		sigma = std::pow(m_vortsheet.ds[i], m_sigma_cutoff);

		// Add the released vortex
		m_vortex.ID.push_back(m_vortex.size() + 1);
		m_vortex.x.push_back(x);
		m_vortex.z.push_back(z);
		m_vortex.circ.push_back(circ);
		m_vortex.sigma.push_back(sigma);

		// Why are these initially zero???
		m_vortex.u.push_back(0);
		m_vortex.w.push_back(0);
		m_vortex.uvs.push_back(0);
		m_vortex.wvs.push_back(0);
		m_vortex.omega.push_back(0);
	}
}

void DVMBase::reflect()
{
	std::vector<unsigned> closest_panel;
	std::vector<double> min_dist;

	closest_panel.resize(m_vortex.size());
	min_dist.resize(m_vortex.size());

	std::vector<double> _mirror;
	_mirror.resize(2);
	double x_init, z_init, x_0, z_0, x_1, z_1;

	double dx, dz, dr, min_prev;

	for (unsigned i = 0; i < m_vortex.size(); i++) {
		min_dist[i] = 0;
		closest_panel[i] = 0;

		if (inside_body(m_vortex.x[i], m_vortex.z[i])) {

			// Find which panel is closest to the vortex
			min_prev = 10E6;
			for (unsigned j = 0; j < m_vortsheet.size(); j++) {

				dx = m_vortsheet.xc[j] - m_vortex.x[i];
				dz = m_vortsheet.zc[j] - m_vortex.z[i];
				dr = std::sqrt(std::pow(dx, 2.0) + std::pow(dz, 2.0));

				if (dr < min_prev) {
					closest_panel[i] = j;
					min_prev = dr;
				}
			}

			// Find the mirror image vortex blob
			x_init = m_vortex.x[i];
			z_init = m_vortex.z[i];

			x_0 = m_body.x[closest_panel[i]];
			z_0 = m_body.z[closest_panel[i]];
			x_1 = m_body.x[closest_panel[i] + 1];
			z_1 = m_body.z[closest_panel[i] + 1];

			_mirror = mirror(x_init, z_init, x_0, z_0, x_1, z_1);

			// std::cout << "Mirrored vortex from " << m_vortex.x[i] << "," <<
			// m_vortex.z[i];

			m_vortex.x[i] = _mirror[0];
			m_vortex.z[i] = _mirror[1];
			m_vortex.circ[i] = m_vortex.circ[i]; // This does not do anything
			                                     // ??? - should this be
			                                     // different ???
			// std::cout<<m_vortex.x[i]<<"\t"<<m_vortex.z[i]<<"\n";

			// std::cout << " to "  << m_vortex.x[i] << "," << m_vortex.z[i] <<
			// std::endl;
		}
	}
}

std::vector<double> DVMBase::mirror(double x_init,
                                    double z_init,
                                    double x_0,
                                    double z_0,
                                    double x_1,
                                    double z_1)
{
	std::vector<double> p2;
	p2.resize(2);
	double dx, dz, a, b;

	dx = x_1 - x_0;
	dz = z_1 - z_0;

	a = (dx * dx - dz * dz) / (dx * dx + dz * dz);
	b = 2.0 * dx * dz / (dx * dx + dz * dz);

	p2[0] = a * (x_init - x_0) + b * (z_init - z_0) + x_0;
	p2[1] = b * (x_init - x_0) - a * (z_init - z_0) + z_0;

	return p2;
}

// Check whether a discrete vortex is inside the body
// (Crossing rule method)
int DVMBase::inside_body(double x, double z)
{
	int cn = 0;

	// loop through all the edges of the polygon
	for (unsigned i = 0; i < m_body.size() - 1; i++) {
		if (((m_body.z[i] <= z) && (m_body.z[i + 1] > z))
		    || ((m_body.z[i] > z) && (m_body.z[i + 1] <= z))) {
			float vt =
			    (float)(z - m_body.z[i]) / (m_body.z[i + 1] - m_body.z[i]);
			if (x < m_body.x[i] + vt * (m_body.x[i + 1] - m_body.x[i])) {
				++cn;
			}
		}
	}
	return (cn & 1);
}

void DVMBase::compute_loads()
{

	double dp;

	std::vector<double> p;
	p.resize(m_vortsheet.size());

	// not so sure here (reference pressure)
	p[0] = 0;
	for (unsigned i = 1; i < m_vortsheet.size(); i++) {

		dp = m_rho / m_dt * m_vortsheet.gamma[i] * m_vortsheet.ds[i];

		p[i] = p[0] - dp;
	}

	m_fx = 0;
	m_fz = 0;
	for (unsigned i = 0; i < m_vortsheet.size(); i++) {
		m_fx += p[i] * m_vortsheet.enx[i] * m_vortsheet.ds[i];
		m_fz += p[i] * m_vortsheet.enz[i] * m_vortsheet.ds[i];
	}

	std::cout << "Fx = " << m_fx << " Fz = " << m_fz << std::endl;
}

void DVMBase::save_vort()
{
	for (unsigned i = 0; i < m_vortsheet.size(); i++) {
		m_vortsheet.gamma_prev[i] = m_vortsheet.gamma[i];
	}
}

void DVMBase::probe_velocities()
{
	double rsigmasqr;
	double rpi2 = 1.0 / (2.0 * m_pi);
	double dx_ij, dz_ij, dr_ij2, threshold, xkernel;
	double dK_ij, zkernel;
	double x_i, z_i, u_i, w_i, c_i;

	// Compute the velocity vector at the probe points

	for (unsigned i = 0; i < m_probe.size(); i++) {
		m_probe.u[i] = 0.0;
		m_probe.w[i] = 0.0;
	}

	for (unsigned i = 0; i < m_probe.size(); i++) {
		x_i = m_probe.x[i];
		z_i = m_probe.z[i];
		u_i = 0.0;
		w_i = 0.0;
		for (unsigned j = 1; j < m_vortex.size(); j++) {
			dx_ij = x_i - m_vortex.x[j];
			dz_ij = z_i - m_vortex.z[j];
			dr_ij2 = std::pow(dx_ij, 2) + std::pow(dz_ij, 2);

			threshold = m_kernel_threshold * std::pow(m_vortex.sigma[j], 2);
			rsigmasqr = 1.0 / std::pow(m_vortex.sigma[j], 2);

			if (dr_ij2 < threshold) {
				dK_ij = (1.0 - std::exp(-dr_ij2 * rsigmasqr)) / dr_ij2;
			} else {
				dK_ij = 1.0 / dr_ij2;
			}

			xkernel = dK_ij * dz_ij;
			zkernel = dK_ij * dx_ij;
			u_i = -xkernel * m_vortex.circ[j];
			w_i = +zkernel * m_vortex.circ[j];
			m_probe.u[i] = m_probe.u[i] + u_i;
			m_probe.w[i] = m_probe.w[i] + w_i;
		}
	}

	double c1, c2, c3, c4, c5, c6, c7, c8, c9;
	Matrix px, py, qx, qy;
	px.resize(m_probe.size(), std::vector<double>(m_vortsheet.size()));
	qx.resize(m_probe.size(), std::vector<double>(m_vortsheet.size()));
	py.resize(m_probe.size(), std::vector<double>(m_vortsheet.size()));
	qy.resize(m_probe.size(), std::vector<double>(m_vortsheet.size()));

	for (unsigned i = 0; i < m_probe.size(); i++) {

		m_probe.uvs[i] = 0.0;
		m_probe.wvs[i] = 0.0;
	}

	for (unsigned i = 0; i < m_probe.size(); i++) {
		for (unsigned j = 1; j < m_vortsheet.size(); j++) {
			c1 = -(m_probe.x[i] - m_vortsheet.xc[j]) * cos(m_vortsheet.theta[j])
			     - (m_probe.z[i] - m_vortsheet.zc[j])
			           * sin(m_vortsheet.theta[j]);
			c2 = std::pow((m_probe.x[i] - m_vortsheet.xc[j]), 2)
			     + std::pow((m_probe.z[i] - m_vortsheet.zc[j]), 2);
			c5 = (m_probe.x[i] - m_vortsheet.xc[j]) * sin(m_vortsheet.theta[j])
			     - (m_probe.z[i] - m_vortsheet.zc[j])
			           * cos(m_vortsheet.theta[j]);
			c6 = log(1.0
			         + m_vortsheet.ds[j] * ((m_vortsheet.ds[j] + 2 * c1) / c2));
			c7 = atan2((c5 * m_vortsheet.ds[j]), (c2 + c1 * m_vortsheet.ds[j]));
			c8 = (m_probe.x[i] - m_vortsheet.xc[j])
			         * sin(-2.0 * m_vortsheet.theta[j])
			     + (m_probe.z[i] - m_vortsheet.zc[j])
			           * cos(-2.0 * m_vortsheet.theta[j]);
			c9 = (m_probe.x[i] - m_vortsheet.xc[j])
			         * cos(-2.0 * m_vortsheet.theta[j])
			     + (m_probe.z[i] - m_vortsheet.zc[j])
			           * sin(-2.0 * m_vortsheet.theta[j]);

			qx[i][j] = -sin(m_vortsheet.theta[j])
			           + 0.5 * c8 * c6 / m_vortsheet.ds[j]
			           + (c1 * cos(m_vortsheet.theta[j])
			              + c5 * sin(m_vortsheet.theta[j]))
			                 * c7 / m_vortsheet.ds[j];
			px[i][j] = -0.5 * c6 * sin(m_vortsheet.theta[j])
			           - c7 * cos(m_vortsheet.theta[j]) - qx[i][j];

			qy[i][j] = cos(m_vortsheet.theta[j])
			           + 0.5 * c9 * c6 / m_vortsheet.ds[j]
			           + (c1 * cos(m_vortsheet.theta[j])
			              - c5 * sin(m_vortsheet.theta[j]))
			                 * c7 / m_vortsheet.ds[j];
			py[i][j] = -0.5 * c6 * cos(m_vortsheet.theta[j])
			           - c7 * sin(m_vortsheet.theta[j]) - qy[i][j];

			if (j == m_vortsheet.size() - 1) {

				m_probe.uvs[i] =
				    (px[i][m_vortsheet.size() - 1]
				         * m_vortsheet.gamma[m_vortsheet.size() - 1]
				     + qx[i][m_vortsheet.size() - 1] * m_vortsheet.gamma[0]);
				m_probe.wvs[i] =
				    (py[i][m_vortsheet.size() - 1]
				         * m_vortsheet.gamma[m_vortsheet.size() - 1]
				     - qy[i][m_vortsheet.size() - 1] * m_vortsheet.gamma[0]);

			} else {
				m_probe.uvs[i] = (px[i][j] * m_vortsheet.gamma[j]
				                  + qx[i][j] * m_vortsheet.gamma[j + 1]);
				m_probe.wvs[i] = (py[i][j] * m_vortsheet.gamma[j]
				                  + qy[i][j] * m_vortsheet.gamma[j + 1]);
			}
		}
	}

	for (unsigned i = 0; i < m_probe.size(); i++) {
		m_probe.u[i] = rpi2 * m_probe.u[i] + m_probe.uvs[i] + m_Ux;
		m_probe.w[i] = rpi2 * m_probe.w[i] + m_probe.wvs[i] + m_Uz;
	}
}
