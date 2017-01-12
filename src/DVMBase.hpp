#ifndef DVMBASE
#define DVMBASE

#include "BaseTypes.hpp"
#include "pugiconfig.hpp"
#include "pugixml.hpp"
#include <algorithm>
#include <armadillo>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <string>

class DVMBase
{
	public:
	DVMBase();

	// initialise the object
	void init(pugi::xml_document &);

	// read paramters
	void init_dipole(double zdistance,
	                 double xdistance,
	                 unsigned nvb1,
	                 unsigned nvb2,
	                 double radius1,
	                 double radius2);
	void init_outputs();
	void read_input_coord();

	// form matrices / and computational entinties
	void form_vortex_sheet();
	void compute_influence_matrix();

	// Convect stuff
	void biotsavart();
	void convect(unsigned order);

	// Vortex sheet stuff
	void solvevortexsheet();
	void vortexsheetbc();
	void imagebc();

	// Diffusion stuff
	void diffrw();
	void diffpse(double nu, double dt);
	void diffuse_vs_rw();
	void write_outputs();

	// Housekeeping stuff
	void delete_vort();
	void reflect();
	void absorb();
	int inside_body(double x, double z);
	std::vector<double> mirror(double x_init,
	                           double z_init,
	                           double x_0,
	                           double z_0,
	                           double x_1,
	                           double z_1);

	// Diagnostics
	void save_vort();
	void compute_loads();
	void compute_strouhal();
	void probe_velocities();
	double get_time();
	unsigned get_size();
	unsigned get_vs_size();
	unsigned get_steps();

	void increment_step();

	public:
	VortexBlobs m_vortex;
	VortexSheet m_vortsheet;
	Body m_body;
	VortexBlobs m_probe;

	public:
	Matrix m_infM;
	std::ofstream dev_dvm, dev_Num, dev_gamma, dev_loads, dev_probe;
	// std::vector<double> xi_rw, eta_rw; // Why were these here???
	std::vector<double> uvs, wvs, uI, wI;
	// std::vector<double> m_Gamma_abs;
	double m_dt, m_nu, m_Ux, m_Uz, m_n, m_np;
	double m_kernel_threshold, m_sigma_cutoff;
	double m_time;
	unsigned m_step, m_steps;
	double m_GammaDel, m_maxGamma;
	double m_St;

	// Reading and writing
	std::string m_in_dir;
	std::string m_out_dir;
	std::string m_domain_file;

	// some constants
	double m_pi, m_rpi2;

	// Loads related variables
	double m_fx, m_fz, m_rho;
};

#endif
