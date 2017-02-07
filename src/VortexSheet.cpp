#include "VortexSheet.hpp"

VortexSheet::VortexSheet()
{
}

VortexSheet::VortexSheet(const XmlHandler &xml)
{
	m_rho = xml.getValueAttribute("constants", "density");
	m_nu = xml.getValueAttribute("constants", "nu");
	m_dt = xml.getValueAttribute("time", "dt");
	m_maxGamma= xml.getValueAttribute("constants", "max_Gamma");
	m_maxNumPanelVort= xml.getValueAttribute("constants", "max_NumPanelVort");
	m_cutoff_exp = xml.getValueAttribute("constants", "cutoff_exp");

}


void VortexSheet::resize(unsigned size)
{
	gamma.set_size(size);
	x.set_size(size+1);
	z.set_size(size+1);
	xc.set_size(size);
	zc.set_size(size);
	theta.set_size(size);
	ds.set_size(size);
	enx.set_size(size);
	enz.set_size(size);
	etx.set_size(size);
	etz.set_size(size);
}

VortexBlobs VortexSheet::release_nascent_vortices_rw(Random& _rand)  
{
	//double x, z, circ, sigma;
	auto Nvs = size();

	double R, rrw;

    // First we need to determine how many of these vortices will be created at each 
    // panel. This is in accordance with Morgenthal PhD 4.3.6 Vortex release algorithm (eq. 4.108)
    Vector_un PanelNewVort(Nvs);                      // Panel's new vortices
    Vector PanelCirc = gamma%ds;    // Panel's total circulation
    Vector AbsPanelCirc = arma::abs(PanelCirc/arma::max(PanelCirc)*m_maxNumPanelVort); // Absolute value of the panel vorticity   
    Vector Circ_new(Nvs);
    //GD---> Should not need to do a loop here. Unfortunately for some reason armadillo does not allow me to do 
    //       PanelNewVort=arma::round(AbsPanelCirc)+arma::ones(Nvs,1), perhaps MAB can fix this. 
    for (unsigned i=0;i<Nvs;i++)
    {
    PanelNewVort(i) = std::floor(AbsPanelCirc(i))+1; // Number of released vortices per panel
    assert(PanelNewVort(i)>0);    
    Circ_new(i)=PanelCirc(i)/PanelNewVort(i); // Circulation of each vortex we release from the ith panel 
    }
    
    unsigned Nrv=arma::sum(PanelNewVort); //Number of the new vortices. 
    // if we set the maximum number of vortices from each panel=1, then Nrv=Nsv
    
    // Initialize the vortexblobs ready for release 
    VortexBlobs nascentVort(Nrv);

    //Calculating the position of the newelly released vortices
    double xm,zm;
    unsigned counter=0; 

    for (unsigned i=0; i<Nvs;i++){
        unsigned j=0;
        while (j<PanelNewVort(i))
        {
            //Find locations at the ith panel from which the vortices should be released
            xm=x(i)+ds(i)*etx(i)/(PanelNewVort(i)+1);
            zm=z(i)+ds(i)*etz(i)/(PanelNewVort(i)+1);
            //Here is the tricky bit !!! Morgenthal shows in figure 4.7 that the particles may be released with
            //some random walk in both the normal (to the panel) and the tangential directions. This is
            //wrong according to Chorin 1978. Since we are in the boundary layer, the diffusion process takes place only
            //in the normal direction and not the streamwise. This also according to Prandtl's boundary layer approximation.
            //I will implement the Chorin 1978 here and not the Morgenthal one. In the end, this should not make 
            //a huge difference.
        
            //Create Random walk values
            R=_rand.rand();
            rrw=std::abs(std::sqrt(4.0*m_nu*m_dt*std::log(1.0/R))); //This is half normal distribution only in the ourwards direction
	        // Add the released vortex
            nascentVort.m_x(counter)=xm+rrw*enx(i);
            nascentVort.m_z(counter)=zm+rrw*enz(i);
            nascentVort.m_circ(counter)=Circ_new(i);
            //Now for the cut-off kernel we implement things as suggested by Mirta Perlman 1985 (JCP)
            //using sigma=ds^q where 0.5<q<1. She suggests using q=0.625! This q parameter is the coefficient not the cutoff
            //parameter the inputs file.
            nascentVort.m_sigma(counter)=std::pow(ds(i),m_cutoff_exp);
            counter++;
            assert(counter<=Nrv);
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

				dx = xc(j) - vortex.m_x(i);
				dz = zc(j) - vortex.m_z(i);
				dr = std::sqrt(std::pow(dx, 2.0) + std::pow(dz, 2.0));

				if (dr < min_prev) {
					closest_panel(i) = j;
					min_prev = dr;
				}
			}

			// Find the mirror image vortex blob
			x_init = vortex.m_x(i);
			z_init = vortex.m_z(i);

			x_0 = x(closest_panel(i));
			z_0 = z(closest_panel(i));
			x_1 = x(closest_panel(i) + 1);
			z_1 = z(closest_panel(i) + 1);

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
		if (((z(i) <= zcoor) && (z(i + 1) > zcoor))
		    || ((z(i) > zcoor) && (z(i + 1) <= zcoor))) {
			float vt =
			    (float)(zcoor - z(i)) / (z(i + 1) - z(i));
			if (xcoor < x(i) + vt * (x(i + 1) - x(i))) {
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
	Vector Dp = -m_rho / m_dt * (gamma % ds);

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
    m_fx = -arma::sum(P % enx % ds);
	m_fz = -arma::sum(P % enz % ds);
    
    std::cout<<"C_D = "<<m_fx/(0.5*m_rho*1.0*Ur*Ur)<<"\t"<<"C_L = "<<m_fz/(0.5*m_rho*1.0*Ur*Ur)<<std::endl;

}

unsigned VortexSheet::size()
{
	if ((gamma.size() != xc.size())
	    && (gamma.size() != zc.size())
	    && (gamma.size() != x.size()-1)
	    && (gamma.size() != z.size()-1)
	    && (gamma.size() != ds.size())
	    && (gamma.size() != theta.size())
	    && (gamma.size() != enx.size())
	    && (gamma.size() != enz.size())
	    && (gamma.size() != etx.size())
	    && (gamma.size() != etz.size())) {
		throw std::string("Size mismatch in VortexSheet");
	}
	return gamma.size();
}

void VortexSheet::print_collocation()
{
	for (unsigned i = 0; i < size(); i++) {
		std::cout << "(xc,zc) = (" << xc(i) << "," << zc(i) << ")" << std::endl;
	}
}

void VortexSheet::print_unit_vectors()
{
	for (unsigned i = 0; i < size(); i++) {
		std::cout << "For panel No " << i << "the normal vectors are en = ("
		          << enx(i) << "," << enz(i) << ") and et = (" << etx(i) << ","
		          << etz(i) << ")" << std::endl;
	}
}

void VortexSheet::print_gamma()
{
	gamma.print();
}

std::tuple<double, double> VortexSheet::get_forces()
{
	return std::make_tuple(m_fx, m_fz);
}
