#include "DVMBase.hpp"
#include <armadillo>
#include <cmath>
#include <iomanip>
#include <time.h>
#include <iostream>
#include <cassert>

DVMBase::DVMBase(XmlHandler &xml) : m_vortex(xml), m_vortsheet(xml), m_body(xml)
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
	m_maxGamma = getVal("constants", "max_Gamma");
	m_maxNumPanelVort = getVal("constants", "max_NumPanelVort");
	m_cutoff_exp = getVal("constants", "cutoff_exp");
	// Need to find a way in the xml handler to deal with this too
	if (m_cutoff_exp >= 1) {
		throw std::string(
		    "<IC2DDVM><constants><cutoff_exponent> must be less than one and greater than 0.5");
	}

	m_dt = getVal("time", "dt");
	m_steps = getVal("time", "steps");

	m_Ux = getVal("flow", "ux");
	m_Uz = getVal("flow", "uz");
	m_Ur = std::sqrt(m_Ux * m_Ux + m_Uz * m_Uz);

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
    
    std::string surfacecross = "REFLECT";//getStr("algorithms", "surface_crossing");
    if (surfacecross.compare("DELETE") == 0) {
		m_surfcross = SurfaceCross::DELETE;
	} else if (surfacecross.compare("ABSORB") == 0) {
		m_surfcross = SurfaceCross::ABSORB;
	} else if (surfacecross.compare("REFLECT") == 0) {
		m_surfcross = SurfaceCross::REFLECT;
	} // Invalid cases dealt with by the xml handler
	
}

void DVMBase::solve()
{
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
	m_vortsheet.solvevortexsheet(m_vortex);

	m_vortsheet.compute_loads(m_Ur);
	
	convect();
    
    // The diffusion substep is split into two steps
    // A diffussion problem with only a flux of vorticity in the boundaries dgamma/dn=a 
    VortexBlobs NewVortices=m_vortsheet.release_nascent_vortices_rw(m_rand); 
    m_vortex.append_vortices(NewVortices); 
	
    m_vortex.diffusion_random_walk(m_rand,m_nu,m_dt); // A diffussion problem in an infinite domain
    
    // If a large time step is used some vortices may cross the boundary due to random walk! 
    // Care is taken of these vortices by 1) Deleting them, 2) Absorbing them, 3) Reflecting them back to the flow
    m_vortsheet.reflect(m_vortex);		
	
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
	          << "\tF_x [-]"
	          << "\tF_z [-]" << std::endl;

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
		dev_gamma << m_time << "\t " << m_vortsheet.m_xc(i) << "\t "
		          << m_vortsheet.m_zc(i) << " " << m_vortsheet.m_gamma(i) << "\t"
		          << m_vortsheet.m_ds(i) << std::endl;
	}

	double fx, fz;
	std::tie(fx, fz) = m_vortsheet.get_forces();
	dev_loads << m_step << " " << fx << "\t" << fz << std::endl;
    //dev_loads<<m_step<<"  "<<fx/(0.5*m_rho*1.0*m_Ux*m_Ux)<<"\t"<<fz/(0.5*m_rho*1.0*m_Ux*m_Ux)<<std::endl;

	for (unsigned i = 0; i < m_probe.size(); i++) {
		dev_probe << m_time << " " << m_probe.m_u(i) << " " << m_probe.m_w(i)
		          << std::endl;
	}
}

void DVMBase::convect()
{
	switch (m_scheme) {
	case Scheme::Euler:
		m_vortex.biotsavart();
		m_vortsheet.vortexsheetbc(m_vortex);

		m_vortex.m_x += (m_vortex.m_u + m_vortex.m_uvs + m_Ux) * m_dt;
		m_vortex.m_z += (m_vortex.m_w + m_vortex.m_wvs + m_Uz) * m_dt;
		break;
	default:
		throw std::string("nothing else implemented yet");
	}
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
		m_probe.m_u(i) = 0.0;
		m_probe.m_w(i) = 0.0;
	}

	for (unsigned i = 0; i < m_probe.size(); i++) {
		x_i = m_probe.m_x(i);
		z_i = m_probe.m_z(i);
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
			m_probe.m_u(i) = m_probe.m_u(i) + u_i;
			m_probe.m_w(i) = m_probe.m_w(i) + w_i;
		}
	}

	double c1, c2, c5, c6, c7, c8, c9;
	Matrix px, py, qx, qy;
	px.set_size(m_probe.size(), m_vortsheet.size());
	qx.set_size(m_probe.size(), m_vortsheet.size());
	py.set_size(m_probe.size(), m_vortsheet.size());
	qy.set_size(m_probe.size(), m_vortsheet.size());

	for (unsigned i = 0; i < m_probe.size(); i++) {

		m_probe.m_uvs(i) = 0.0;
		m_probe.m_wvs(i) = 0.0;
	}

	for (unsigned i = 0; i < m_probe.size(); i++) {
		for (unsigned j = 1; j < m_vortsheet.size(); j++) {
			c1 = -(m_probe.m_x(i) - m_vortsheet.m_xc(j)) * cos(m_vortsheet.m_theta(j))
			     - (m_probe.m_z(i) - m_vortsheet.m_zc(j))
			           * sin(m_vortsheet.m_theta(j));
			c2 = std::pow((m_probe.m_x(i) - m_vortsheet.m_xc(j)), 2)
			     + std::pow((m_probe.m_z(i) - m_vortsheet.m_zc(j)), 2);
			c5 = (m_probe.m_x(i) - m_vortsheet.m_xc(j)) * sin(m_vortsheet.m_theta(j))
			     - (m_probe.m_z(i) - m_vortsheet.m_zc(j))
			           * cos(m_vortsheet.m_theta(j));
			c6 = log(1.0
			         + m_vortsheet.m_ds(j) * ((m_vortsheet.m_ds(j) + 2 * c1) / c2));
			c7 = atan2((c5 * m_vortsheet.m_ds(j)), (c2 + c1 * m_vortsheet.m_ds(j)));
			c8 = (m_probe.m_x(i) - m_vortsheet.m_xc(j))
			         * sin(-2.0 * m_vortsheet.m_theta(j))
			     + (m_probe.m_z(i) - m_vortsheet.m_zc(j))
			           * cos(-2.0 * m_vortsheet.m_theta(j));
			c9 = (m_probe.m_x(i) - m_vortsheet.m_xc(j))
			         * cos(-2.0 * m_vortsheet.m_theta(j))
			     + (m_probe.m_z(i) - m_vortsheet.m_zc(j))
			           * sin(-2.0 * m_vortsheet.m_theta(j));

			qx(i, j) = -sin(m_vortsheet.m_theta(j))
			           + 0.5 * c8 * c6 / m_vortsheet.m_ds(j)
			           + (c1 * cos(m_vortsheet.m_theta(j))
			              + c5 * sin(m_vortsheet.m_theta(j)))
			                 * c7 / m_vortsheet.m_ds(j);
			px(i, j) = -0.5 * c6 * sin(m_vortsheet.m_theta(j))
			           - c7 * cos(m_vortsheet.m_theta(j)) - qx(i, j);

			qy(i, j) = cos(m_vortsheet.m_theta(j))
			           + 0.5 * c9 * c6 / m_vortsheet.m_ds(j)
			           + (c1 * cos(m_vortsheet.m_theta(j))
			              - c5 * sin(m_vortsheet.m_theta(j)))
			                 * c7 / m_vortsheet.m_ds(j);
			py(i, j) = -0.5 * c6 * cos(m_vortsheet.m_theta(j))
			           - c7 * sin(m_vortsheet.m_theta(j)) - qy(i, j);

			if (j == m_vortsheet.size() - 1) {

				m_probe.m_uvs(i) =
				    (px(i, m_vortsheet.size() - 1)
				         * m_vortsheet.m_gamma(m_vortsheet.size() - 1)
				     + qx(i, m_vortsheet.size() - 1) * m_vortsheet.m_gamma(0));
				m_probe.m_wvs(i) =
				    (py(i, m_vortsheet.size() - 1)
				         * m_vortsheet.m_gamma(m_vortsheet.size() - 1)
				     - qy(i, m_vortsheet.size() - 1) * m_vortsheet.m_gamma(0));

			} else {
				m_probe.m_uvs(i) = (px(i, j) * m_vortsheet.m_gamma(j)
				                  + qx(i, j) * m_vortsheet.m_gamma(j + 1));
				m_probe.m_wvs(i) = (py(i, j) * m_vortsheet.m_gamma(j)
				                  + qy(i, j) * m_vortsheet.m_gamma(j + 1));
			}
		}
	}

	for (unsigned i = 0; i < m_probe.size(); i++) {
		m_probe.m_u(i) = rpi2 * m_probe.m_u(i) + m_probe.m_uvs(i) + m_Ux;
		m_probe.m_w(i) = rpi2 * m_probe.m_w(i) + m_probe.m_wvs(i) + m_Uz;
	}
}
