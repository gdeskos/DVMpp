#include "DVMBase.hpp"
#include <armadillo>
#include <cmath>
#include <iomanip>
#include <time.h>
#include <iostream>
#include <cassert>

DVMBase::DVMBase(XmlHandler &xml) : m_vortex(xml), m_vortsheet(xml)
{
	m_pi = 4.0 * atan(1.0);
	m_step = 0;
	m_rpi2 = 1.0 / (2.0 * m_pi);
}

void DVMBase::init(XmlHandler &xml, std::string timestamp)
{
	m_timestamp = timestamp;

	// Some helpers so we can do less typing
	auto getVal = [&xml](const char *l, const char *p) {
		return xml.getValueAttribute(l, p);
	};

	auto getStr = [&xml](const char *l, const char *p) {
		return xml.getStringAttribute(l, p);
	};

	auto checksep = [](std::string str) {
		return (str.back() == '/') ? str : str + '/';
	};

	// Get all the various inputs. Valid inputs are checked by the xml handler
	m_in_dir = checksep(getStr("io", "input_dir"));
	m_out_dir = checksep(getStr("io", "output_dir"));
	m_domain_file = getStr("io", "domain_file");

	m_rho = getVal("constants", "density");
	m_nu = getVal("constants", "nu");
	m_maxNumPanelVort = getVal("constants", "max_NumPanelVort");
	m_cutoff_exp = getVal("constants", "cutoff_exp");
	// Need to find a way in the xml handler to deal with this too
	if (m_cutoff_exp >= 1 ) {
		throw std::string(
		    "<IC2DDVM><constants><cutoff_exponent> must be less than one and greater than 0.5");
	}

	m_dt = getVal("time", "dt");
	m_steps = getVal("time", "steps");

	m_Ux = getVal("flow", "ux");
	m_Uz = getVal("flow", "uz");

	m_probe.m_x = xml.getList("probe", "x");
	m_probe.m_z = xml.getList("probe", "z");

	auto seed = xml.getIntAttribute("constants", "seed");
	m_rand.seed(seed);

	auto scheme = getStr("time", "scheme");
	if (scheme.compare("euler") == 0) {
		m_scheme = Scheme::Euler;
	} else if (scheme.compare("RK2") == 0) {
		m_scheme = Scheme::RK2;
	} else if (scheme.compare("RK4") == 0) {
		m_scheme = Scheme::RK4;
	} // Invalid cases dealt with by the xml handler
}

void DVMBase::solve()
{
	// First timestep
	form_vortex_sheet();
	compute_influence_matrix();

	// Timeloop
	for (unsigned j = 1; j <= get_steps(); j++) {

		increment_step();
		compute_step();

		// Output
		write_outputs();

		// Screen output
		std::cout << "Simulation time          = " << get_time() << "\tStep "
		          << j << "/" << get_steps() << std::endl;
		std::cout << "Number of vortex blobs   = " << get_size() << std::endl;
	}
}

void DVMBase::compute_step()
{
	// Inviscid Substep
	solvevortexsheet();

    double Urel;
    Urel=std::sqrt(m_Ux*m_Ux+m_Uz*m_Uz);
	m_vortsheet.compute_loads(Urel);
    //This should not be here! This is just a test
	double fx, fz;
	std::tie(fx, fz) = m_vortsheet.get_forces();
    std::cout<<"C_D = "<<fx/(0.5*m_rho*1.0*Urel*Urel)<<"\t"<<"C_L = "<<fz/(0.5*m_rho*1.0*Urel*Urel)<<std::endl;

	convect();
    
    // The diffusion substep is split into two steps
    // A diffussion problem with only a flux of vorticity in the boundaries dgamma/dn=a 
    VortexBlobs NewVortices=m_vortsheet.release_nascent_vortices_rw(m_rand); 
    std::cout<<NewVortices.size()<<std::endl;
    m_vortex.append_vortices(NewVortices); 

	diffrw();                    // A diffussion problem in an infinite domain

	// Housekeeping
    // For a large time step vortices may cross the boundary due to random walk! We reflect them back
	reflect();
}

void DVMBase::read_input_coord()
{
	// Read the body coordinates
	std::string file = m_in_dir + m_domain_file;
	std::string line;
	std::ifstream coor_file(file.c_str());
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
	dev_dvm.open(
	    (m_out_dir + m_timestamp + std::string("_vortex.dat")).c_str());
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
	    (m_out_dir + m_timestamp + std::string("_vortex_num.dat")).c_str());
	dev_Num << "Time [s]"
	        << " "
	        << "Number of vortices" << std::endl;

	dev_gamma.open(
	    (m_out_dir + m_timestamp + std::string("_gamma.dat")).c_str());
	dev_gamma << m_vortsheet.size() << " # Number of collocation points"
	          << std::endl;
	dev_gamma << m_dt << " # Time Step" << std::endl;
	dev_gamma << m_steps << " # Steps " << std::endl;
	dev_gamma << "Time [s]"
	          << " "
	          << "Gamma - Vortex sheet strength" << std::endl;

	dev_loads.open(
	    (m_out_dir + m_timestamp + std::string("_loads.dat")).c_str());
	dev_loads << "Time [s]"
	          << " "
	          << "\tC_x [-]"
	          << "\tC_z [-]" << std::endl;

	dev_probe.open(
	    (m_out_dir + m_timestamp + std::string("_probe.dat")).c_str());
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
        
        // Assign the body points to the vortsheet, remember this is m_vortsheet.size()+1
		m_vortsheet.x[i]=m_body.x[i];
		m_vortsheet.z[i]=m_body.z[i];
		m_vortsheet.x[i+1]=m_body.x[i+1];
		m_vortsheet.z[i+1]=m_body.z[i+1];
        
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

	m_infM.set_size(Nl + 1, Nl);

	std::cout << "m_infM.size " << m_infM.size() << std::endl;

	double c1, c2, c3, c4, c5, c6, c7, c9;
	Matrix p, q;
	p.set_size(Nl, Nl);
	q.set_size(Nl, Nl);

	for (unsigned i = 0; i < Nl; i++) {

		double xci = m_vortsheet.xc[i];
		double zci = m_vortsheet.zc[i];
		double thetai = m_vortsheet.theta[i];

		for (unsigned j = 0; j < Nl; j++) {
			if (i == j) {
				p(i, j) = -1.0;
				q(i, j) = 1.0;
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
				c9 = (xci - xj) * cos(thetai - 2.0 * thetaj)
				     - (zci - zj) * sin(thetai - 2.0 * thetaj);

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
		m_infM(Nl, j) = m_vortsheet.ds[j];
	}

	// Print the matrix
	m_infM.print();

	std::cout << "Computed influence matrix of size " << m_infM.size()
	          << " after Kuette and Chow" << std::endl;
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
		dev_dvm << m_time << " " << m_vortex.m_x(i) << " " << m_vortex.m_z(i) << " "
		        << m_vortex.m_circ(i) << std::endl;
	}

	dev_Num << m_step << " " << m_vortex.size() << std::endl;

	for (unsigned i = 0; i < m_vortsheet.size(); i++) {
		dev_gamma << m_time << "\t " << m_vortsheet.xc[i] << "\t "
		          << m_vortsheet.zc[i] << " " << m_vortsheet.gamma[i] << "\t"
		          << m_vortsheet.ds[i] << std::endl;
	}

	double fx, fz;
	std::tie(fx, fz) = m_vortsheet.get_forces();
	//dev_loads << m_step << " " << fx << "\t" << fz << std::endl;
    dev_loads<<m_step<<"  "<<fx/(0.5*m_rho*1.0*m_Ux*m_Ux)<<"\t"<<fz/(0.5*m_rho*1.0*m_Ux*m_Ux)<<std::endl;

	for (unsigned i = 0; i < m_probe.size(); i++) {
		dev_probe << m_time << " " << m_probe.m_u[i] << " " << m_probe.m_w[i]
		          << std::endl;
	}
}

void DVMBase::convect()
{
	switch (m_scheme) {
	case Scheme::Euler:
		m_vortex.biotsavart();
		vortexsheetbc();

		m_vortex.m_x += (m_vortex.m_u + m_vortex.m_uvs + m_Ux) * m_dt;
		m_vortex.m_z += (m_vortex.m_w + m_vortex.m_wvs + m_Uz) * m_dt;
		break;
	default:
		throw std::string("nothing else implemented yet");
	}
}

void DVMBase::solvevortexsheet()
{

	unsigned Nl = m_vortsheet.size();

	Vector brhs(Nl + 1);

	if (m_vortex.size() == 0) {
		brhs.rows(0, Nl - 1) = m_Ux * m_vortsheet.enx + m_Uz * m_vortsheet.enz;
		brhs(Nl) = -m_vortex.totalcirc();

	} else {
		Vector u(Nl);
		Vector w(Nl);

#pragma omp parallel for
		for (unsigned i = 0; i < m_vortsheet.size(); i++) {

			double dK_ij, rsigmasqr, dx_ij, dz_ij, dr_ij2, threshold;

			u(i) = 0.0;
			w(i) = 0.0;

			for (unsigned j = 0; j < m_vortex.size(); j++) {

				dx_ij = m_vortsheet.xc[i] - m_vortex.m_x[j];
				dz_ij = m_vortsheet.zc[i] - m_vortex.m_z[j];
				dr_ij2 = std::pow(dx_ij, 2) + std::pow(dz_ij, 2);

				threshold =
				    10*std::pow(m_vortex.m_sigma[j], 2);
				rsigmasqr = 1.0 / std::pow(m_vortex.m_sigma[j], 2);

				if (dr_ij2 < threshold) {
					dK_ij = (1.0 - std::exp(-dr_ij2 * rsigmasqr)) / dr_ij2;
				} else {
					dK_ij = 1.0 / dr_ij2;
				}

				u(i) -= dK_ij * dz_ij * m_vortex.m_circ[j];
				w(i) += dK_ij * dx_ij * m_vortex.m_circ[j];
			}
		}

		u *= m_rpi2;
		w *= m_rpi2;

		// Not entirely convinced that this is the correct BC (see Morgenthal)
		brhs.rows(0, Nl - 1) =
		    (m_Ux + u) % m_vortsheet.enx + (m_Uz + w) % m_vortsheet.enz;
		brhs(Nl) = -m_vortex.totalcirc();
	}

	// Solve system
	m_vortsheet.gamma = arma::solve(m_infM.t() * m_infM, m_infM.t() * brhs);
}

void DVMBase::vortexsheetbc()
{
	double c1, c2, c5, c6, c7, c8, c9;
	Matrix px, py, qx, qy;

	px.set_size(m_vortex.size(), m_vortsheet.size());
	qx.set_size(m_vortex.size(), m_vortsheet.size());
	py.set_size(m_vortex.size(), m_vortsheet.size());
	qy.set_size(m_vortex.size(), m_vortsheet.size());

// compute the coefficients
	for (unsigned i = 0; i < m_vortex.size(); i++) {

		double xi = m_vortex.m_x[i];
		double zi = m_vortex.m_z[i];

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

			qx(i, j) = -sin(thetaj) + 0.5 * c8 * c6 / dsj
			           + (c1 * cos(thetaj) + c5 * sin(thetaj)) * c7 / dsj;
			px(i, j) = -0.5 * c6 * sin(thetaj) - c7 * cos(thetaj) - qx(i, j);

			qy(i, j) = cos(thetaj) + 0.5 * c9 * c6 / dsj
			           + (c1 * sin(thetaj) - c5 * cos(thetaj)) * c7 / dsj;
			py(i, j) = 0.5 * c6 * cos(thetaj) - c7 * sin(thetaj) - qy(i, j);
		}
	}

	// calculate the vortexsheet induced velocities
	unsigned last = m_vortsheet.size() - 1;

#pragma omp parallel for
	for (unsigned i = 0; i < m_vortex.size(); i++) {

		m_vortex.m_uvs[i] = 0.0;
		m_vortex.m_wvs[i] = 0.0;

		for (unsigned j = 0; j < m_vortsheet.size(); j++) {
			if (j == last) {
				m_vortex.m_uvs[i] += px(i, last) * m_vortsheet.gamma[last]
				                   + qx(i, last) * m_vortsheet.gamma[0];
				m_vortex.m_wvs[i] += py(i, last) * m_vortsheet.gamma[last]
				                   + qy(i, last) * m_vortsheet.gamma[0];
			} else {
				m_vortex.m_uvs[i] += px(i, j) * m_vortsheet.gamma[j]
				                   + qx(i, j) * m_vortsheet.gamma[j + 1];
				m_vortex.m_wvs[i] += py(i, j) * m_vortsheet.gamma[j]
				                   + qy(i, j) * m_vortsheet.gamma[j + 1];
			}
		}
		m_vortex.m_uvs[i] *= m_rpi2;
		m_vortex.m_wvs[i] *= m_rpi2;
	}
}

void DVMBase::diffrw()
{
	double R1, R2, rrw, thetarw;

	for (unsigned i = 0; i < m_vortex.size(); i++) {

		// Generate two random numbers in the range 0...1
		R1 = m_rand.rand();
		R2 = m_rand.rand();

		// Calculate r and theta for the random walk
		rrw = std::sqrt(4.0 * m_nu * m_dt * std::log(1.0 / R1));
		thetarw = 2.0 * m_pi * R2;

		m_vortex.m_x[i] += rrw * cos(thetarw);
		m_vortex.m_z[i] += rrw * sin(thetarw);
	}
}

void DVMBase::diffuse_vs_rw() // I will change its name to releave vortices 
{
	//double x, z, circ, sigma;
	auto Nvs = m_vortsheet.size();

	double R, rrw;

    // First we need to determine how many of these vortices will be created at each 
    // panel. This is in accordance with Morgenthal PhD 4.3.6 Vortex release algorithm (eq. 4.108)
    Vector_un PanelNewVort(Nvs);                            // Panel's new vortices
    Vector PanelCirc = m_vortsheet.gamma%m_vortsheet.ds;    // Panel's total circulation
    Vector AbsPanelCirc = arma::abs(PanelCirc/arma::max(PanelCirc)*m_maxNumPanelVort); // Absolute value of the panel vorticity   
    Vector Circ_new(Nvs);
    //GD---> Should not need to do a loop here. Unfortunately for some reason armadillo does not allow me to do 
    //       PanelNewVort=arma::round(AbsPanelCirc)+arma::ones(Nvs,1), perhaps MAB can fix this. 
    for (unsigned i=0;i<Nvs;i++)
    {
    PanelNewVort[i] = std::floor(AbsPanelCirc[i])+1; // Number of released vortices per panel
    assert(PanelNewVort[i]>0);    
    Circ_new[i]=PanelCirc[i]/PanelNewVort[i]; // Circulation of each vortex we release from the ith panel 
    }
    
    unsigned Nrv=arma::sum(PanelNewVort); //Number of the new vortices. 
    // if we set the maximum number of vortices from each panel=1, then Nrv=Nsv
    Vector x_vec(Nrv);
	Vector z_vec(Nrv);
	Vector circ_vec(Nrv);
	Vector sigma_vec(Nrv); 
    unsigned counter=0; 

    //Calculating the position of the newelly released vortices
    double xstart,zstart,xm,zm;

    for (unsigned i=0; i<Nvs;i++){
        xstart=m_vortsheet.x[i];zstart=m_vortsheet.z[i]; // Coordinates for the start point of the panel
        unsigned j=0;
        while (j<PanelNewVort[i])
        {
            //Find locations at the ith panel from which the vortices should be released
            xm=xstart+m_vortsheet.ds[i]*m_vortsheet.etx[i]/(PanelNewVort[i]+1);
            zm=zstart+m_vortsheet.ds[i]*m_vortsheet.etz[i]/(PanelNewVort[i]+1);
            //Here is the tricky bit !!! Morgenthal shows in figure 4.7 that the particles may be released with
            //some random walk in both the normal (to the panel) and the tangential directions. This is
            //wrong according to Chorin 1978. Since we are in the boundary layer, the diffusion process takes place only
            //in the normal direction and not the streamwise. This also according to Prandtl's boundary layer approximation.
            //I will implement the Chorin 1978 here and not the Morgenthal one. In the end, this should not make 
            //a huge difference.
        
            //Create Random walk values
            R=m_rand.rand();
            rrw=std::abs(std::sqrt(4.0*m_nu*m_dt*std::log(1.0/R))); //This is half normal distribution only in the ourwards direction
	        // Add the released vortex
            m_vortex.m_ID.push_back(m_vortex.m_ID.size() + 1); //add the id
            x_vec(counter)=xm+rrw*m_vortsheet.enx[i];
            z_vec(counter)=zm+rrw*m_vortsheet.enz[i];
            circ_vec(counter)=Circ_new[i];
            //Now for the cut-off kernel we implement things as suggested by Mirta Perlman 1985 (JCP)
            //using sigma=ds^q where 0.5<q<1. She suggests using q=0.625! This q parameter is the coefficient not the cutoff
            //parameter the inputs file.
            sigma_vec(counter)=std::pow(m_vortsheet.ds[i],m_cutoff_exp);
            counter++;
            assert(counter<=Nrv);
            j++;
        }   
            
    }
	        m_vortex.m_x.insert_rows(m_vortex.m_x.n_elem, x_vec);
	        m_vortex.m_z.insert_rows(m_vortex.m_z.n_elem, z_vec);
	        m_vortex.m_circ.insert_rows(m_vortex.m_circ.n_elem, circ_vec);
	        m_vortex.m_sigma.insert_rows(m_vortex.m_sigma.n_elem, sigma_vec);
            
            //And then we initialize everything else to zero
	        m_vortex.m_u.insert_rows(m_vortex.m_u.n_elem, Nrv, true);
	        m_vortex.m_w.insert_rows(m_vortex.m_w.n_elem, Nrv, true);
	        m_vortex.m_uvs.insert_rows(m_vortex.m_uvs.n_elem, Nrv, true);
	        m_vortex.m_wvs.insert_rows(m_vortex.m_wvs.n_elem, Nrv, true);
	        m_vortex.m_omega.insert_rows(m_vortex.m_omega.n_elem, Nrv, true); 
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

		if (inside_body(m_vortex.m_x(i), m_vortex.m_z(i))) {

			// Find which panel is closest to the vortex
			min_prev = 10E6;
			for (unsigned j = 0; j < m_vortsheet.size(); j++) {

				dx = m_vortsheet.xc[j] - m_vortex.m_x(i);
				dz = m_vortsheet.zc[j] - m_vortex.m_z(i);
				dr = std::sqrt(std::pow(dx, 2.0) + std::pow(dz, 2.0));

				if (dr < min_prev) {
					closest_panel[i] = j;
					min_prev = dr;
				}
			}

			// Find the mirror image vortex blob
			x_init = m_vortex.m_x(i);
			z_init = m_vortex.m_z(i);

			x_0 = m_body.x[closest_panel[i]];
			z_0 = m_body.z[closest_panel[i]];
			x_1 = m_body.x[closest_panel[i] + 1];
			z_1 = m_body.z[closest_panel[i] + 1];

			_mirror = mirror(x_init, z_init, x_0, z_0, x_1, z_1);

			// std::cout << "Mirrored vortex from " << m_vortex.x[i] << "," <<
			// m_vortex.z[i];

			m_vortex.m_x(i) = _mirror[0];
			m_vortex.m_z(i) = _mirror[1];
			m_vortex.m_circ(i) = m_vortex.m_circ(i); // This does not do anything
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

void DVMBase::probe_velocities()
{
	double rsigmasqr;
	double rpi2 = 1.0 / (2.0 * m_pi);
	double dx_ij, dz_ij, dr_ij2, threshold, xkernel;
	double dK_ij, zkernel;
	double x_i, z_i, u_i, w_i;

	// Compute the velocity vector at the probe points

	for (unsigned i = 0; i < m_probe.size(); i++) {
		m_probe.m_u[i] = 0.0;
		m_probe.m_w[i] = 0.0;
	}

	for (unsigned i = 0; i < m_probe.size(); i++) {
		x_i = m_probe.m_x[i];
		z_i = m_probe.m_z[i];
		u_i = 0.0;
		w_i = 0.0;
		for (unsigned j = 1; j < m_vortex.size(); j++) {
			dx_ij = x_i - m_vortex.m_x(j);
			dz_ij = z_i - m_vortex.m_z(j);
			dr_ij2 = std::pow(dx_ij, 2) + std::pow(dz_ij, 2);

			threshold = 10 * std::pow(m_vortex.m_sigma(j), 2);
			rsigmasqr = 1.0 / std::pow(m_vortex.m_sigma(j), 2);

			if (dr_ij2 < threshold) {
				dK_ij = (1.0 - std::exp(-dr_ij2 * rsigmasqr)) / dr_ij2;
			} else {
				dK_ij = 1.0 / dr_ij2;
			}

			xkernel = dK_ij * dz_ij;
			zkernel = dK_ij * dx_ij;
			u_i = -xkernel * m_vortex.m_circ(j);
			w_i = +zkernel * m_vortex.m_circ(j);
			m_probe.m_u[i] = m_probe.m_u[i] + u_i;
			m_probe.m_w[i] = m_probe.m_w[i] + w_i;
		}
	}

	double c1, c2, c5, c6, c7, c8, c9;
	Matrix px, py, qx, qy;
	px.set_size(m_probe.size(), m_vortsheet.size());
	qx.set_size(m_probe.size(), m_vortsheet.size());
	py.set_size(m_probe.size(), m_vortsheet.size());
	qy.set_size(m_probe.size(), m_vortsheet.size());

	for (unsigned i = 0; i < m_probe.size(); i++) {

		m_probe.m_uvs[i] = 0.0;
		m_probe.m_wvs[i] = 0.0;
	}

	for (unsigned i = 0; i < m_probe.size(); i++) {
		for (unsigned j = 1; j < m_vortsheet.size(); j++) {
			c1 = -(m_probe.m_x[i] - m_vortsheet.xc[j]) * cos(m_vortsheet.theta[j])
			     - (m_probe.m_z[i] - m_vortsheet.zc[j])
			           * sin(m_vortsheet.theta[j]);
			c2 = std::pow((m_probe.m_x[i] - m_vortsheet.xc[j]), 2)
			     + std::pow((m_probe.m_z[i] - m_vortsheet.zc[j]), 2);
			c5 = (m_probe.m_x[i] - m_vortsheet.xc[j]) * sin(m_vortsheet.theta[j])
			     - (m_probe.m_z[i] - m_vortsheet.zc[j])
			           * cos(m_vortsheet.theta[j]);
			c6 = log(1.0
			         + m_vortsheet.ds[j] * ((m_vortsheet.ds[j] + 2 * c1) / c2));
			c7 = atan2((c5 * m_vortsheet.ds[j]), (c2 + c1 * m_vortsheet.ds[j]));
			c8 = (m_probe.m_x[i] - m_vortsheet.xc[j])
			         * sin(-2.0 * m_vortsheet.theta[j])
			     + (m_probe.m_z[i] - m_vortsheet.zc[j])
			           * cos(-2.0 * m_vortsheet.theta[j]);
			c9 = (m_probe.m_x[i] - m_vortsheet.xc[j])
			         * cos(-2.0 * m_vortsheet.theta[j])
			     + (m_probe.m_z[i] - m_vortsheet.zc[j])
			           * sin(-2.0 * m_vortsheet.theta[j]);

			qx(i, j) = -sin(m_vortsheet.theta[j])
			           + 0.5 * c8 * c6 / m_vortsheet.ds[j]
			           + (c1 * cos(m_vortsheet.theta[j])
			              + c5 * sin(m_vortsheet.theta[j]))
			                 * c7 / m_vortsheet.ds[j];
			px(i, j) = -0.5 * c6 * sin(m_vortsheet.theta[j])
			           - c7 * cos(m_vortsheet.theta[j]) - qx(i, j);

			qy(i, j) = cos(m_vortsheet.theta[j])
			           + 0.5 * c9 * c6 / m_vortsheet.ds[j]
			           + (c1 * cos(m_vortsheet.theta[j])
			              - c5 * sin(m_vortsheet.theta[j]))
			                 * c7 / m_vortsheet.ds[j];
			py(i, j) = -0.5 * c6 * cos(m_vortsheet.theta[j])
			           - c7 * sin(m_vortsheet.theta[j]) - qy(i, j);

			if (j == m_vortsheet.size() - 1) {

				m_probe.m_uvs[i] =
				    (px(i, m_vortsheet.size() - 1)
				         * m_vortsheet.gamma[m_vortsheet.size() - 1]
				     + qx(i, m_vortsheet.size() - 1) * m_vortsheet.gamma[0]);
				m_probe.m_wvs[i] =
				    (py(i, m_vortsheet.size() - 1)
				         * m_vortsheet.gamma[m_vortsheet.size() - 1]
				     - qy(i, m_vortsheet.size() - 1) * m_vortsheet.gamma[0]);

			} else {
				m_probe.m_uvs[i] = (px(i, j) * m_vortsheet.gamma[j]
				                  + qx(i, j) * m_vortsheet.gamma[j + 1]);
				m_probe.m_wvs[i] = (py(i, j) * m_vortsheet.gamma[j]
				                  + qy(i, j) * m_vortsheet.gamma[j + 1]);
			}
		}
	}

	for (unsigned i = 0; i < m_probe.size(); i++) {
		m_probe.m_u[i] = rpi2 * m_probe.m_u[i] + m_probe.m_uvs[i] + m_Ux;
		m_probe.m_w[i] = rpi2 * m_probe.m_w[i] + m_probe.m_wvs[i] + m_Uz;
	}
}
