#include "BaseTypes.hpp"



//---------------- Vortex Blobs ---------------------//

VortexBlobs::VortexBlobs()
{
}

void VortexBlobs::resize(unsigned size)
{
    ID.resize(size);
    x.resize(size); 
    z.resize(size);
    circ.resize(size);
    sigma.resize(size);
    u.resize(size);
    w.resize(size);
    uvs.resize(size);
    wvs.resize(size); 
    omega.resize(size);
}	

double VortexBlobs::totalcirc()
{ 
    double sum=0.0;
    for(unsigned i=0; i<size(); i++){
        sum=sum+circ[i];
    }
    return sum;
}


unsigned VortexBlobs::size()
{ 
	if ((ID.size() != x.size())&&(ID.size()!=z.size())&&(ID.size()!=circ.size())&&(ID.size()!=sigma.size()) &&(ID.size()!=u.size())&&(ID.size()!=w.size())&&  (ID.size()!=uvs.size())&&(ID.size()!=wvs.size())&& ID.size()!=omega.size()) 
    { throw std::string("Size mismatch in VortexBlobs");} 
	return ID.size();
}


void VortexBlobs::print_location()
{ 
	for (unsigned i=0;i<size(); i++){ 
		std::cout <<" The location of Vortex Blob "<<ID[i]<<" is (x,z) = (" << x[i] << "," << z[i] << ")" << std::endl;
	}
}

void VortexBlobs::print_velocity()
{ 
	for (unsigned i=0;i<size(); i++){ 
		std::cout << "(u,w) = (" << u[i] << "," << w[i] << ")" << std::endl;
	}
}

void VortexBlobs::print_circulation()
{ 
	for (unsigned i=0;i<size(); i++){ 
		std::cout << " circ = " << circ[i] << std::endl;
	}
}



//-----------------Vortex Sheet------------------//

VortexSheet::VortexSheet()
{
}

void VortexSheet::resize(unsigned size)
{ 
	gamma.resize(size);
    x.resize(size);
    z.resize(size);
    xc.resize(size); 
	zc.resize(size);
    theta.resize(size);
    ds.resize(size);
    enx.resize(size);
    enz.resize(size);
    etx.resize(size);
    etz.resize(size);
}	

unsigned VortexSheet::size()
{ 
	if ((gamma.size() != xc.size())&&(gamma.size()!=zc.size())&&(gamma.size()!=ds.size())&&(gamma.size()!=theta.size())&&(gamma.size()!=enx.size()) &&(gamma.size()!=enz.size())&&(gamma.size()!=etx.size())&& (gamma.size()!=etz.size())) 
    { throw std::string("Size mismatch in VortexSheet");} 
	return gamma.size();
}

void VortexSheet::print_collocation()
{ 
	for (unsigned i=0;i<size(); i++){ 
		std::cout <<"(xc,zc) = (" << xc[i] << "," << zc[i] << ")" << std::endl;
	}
}

void VortexSheet::print_unit_vectors()
{ 
	for (unsigned i=0;i<size(); i++){ 
		std::cout << "For panel No "<<i<<"the normal vectors are en = (" << enx[i] << "," << enz[i] << ") and et = (" <<etx[i]<<","<<etz[i]<<")"<<std::endl;
        }
}

void VortexSheet::print_gamma()
{ 
	for (unsigned i=0;i<size(); i++){ 
		std::cout << gamma[i] << std::endl;
	}
}

// ---------------  Body ---------------- //

Body::Body()
{
}

void Body::resize(unsigned size)
{
	x.resize(size); 
	z.resize(size);
}	

unsigned Body::size()
{ 
	if (x.size() != z.size() )
    { throw std::string("Size mismatch in VortexBlobs");} 
	return x.size();
}


void Body::print_location()
{ 
	for (unsigned i=0;i<size(); i++){ 
		std::cout <<" Body coor (x,z) = (" << x[i] << "," << z[i] << ")" << std::endl;
	}
}

