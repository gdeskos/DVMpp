#ifndef DVMBASE
#define DVMBASE

#include "BaseTypes.hpp"
#include <iostream>
#include <sstream>
#include <string.h>
#include <fstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <armadillo>

class DVMBase
{
public:
     DVMBase();
     // read paramters 
     void init_dipole(double zdistance, double xdistance, unsigned nvb1, unsigned nvb2,double radius1, double radius2);
     void init_outputs();
     void read_input_coord();
     
     // form matrices / and computational entinties
     void form_vortex_sheet();
     void compute_influence_matrix();
     
     // Convect stuff
     void biotsavart();
     void convect(unsigned order, double dt);
     
     // Vortex sheet stuff
     void solvevortexsheet();
     void vortexsheetbc();
     void imagebc();

     // Diffusion stuff
     void diffrw(double nu,double dt);
     void diffpse(double nu, double dt);
     void diffuse_vs_rw();
     void write_outputs();

     // Housekeeping stuff
     void delete_vort();
     void reflect();
     void absorb();
     int inside_body(double x, double z);
    std::vector<double> mirror(double x_init, double z_init, double x_0, double z_0, double x_1, double z_1);

    // Diagnostics
    void save_vort();
    void compute_loads();
    void compute_strouhal();
    void probe_velocities();

public:
     VortexBlobs m_vortex;
     VortexSheet m_vortsheet;
     Body  m_body;
     VortexBlobs m_probe; 

public:
     Matrix m_infM;
     std::ofstream  dev_dvm, dev_Num, dev_gamma, dev_loads, dev_probe;
     std::vector<double> xi_rw, eta_rw;
     std::vector<double> uvs, wvs, uI, wI;
     std::vector<double> m_gamma_prev;
     std::vector<double> m_Gamma_abs;     
     double m_dt, m_nu, m_Ux, m_Uz, m_n, m_np;
     double m_time, m_step, m_steps;   
     double m_GammaDel, m_maxGamma;
     double m_St;

     // some constants
     double m_pi;

     // Loads related variables
     double m_CD, m_rho;
};


#endif

