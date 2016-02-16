#ifndef BASETYPES
#define BASETYPES

#include <vector>
#include <cstring>
#include <iostream>
#include <cstdlib>


typedef std::vector<double> Vector;
typedef std::vector<Vector> Matrix;
typedef std::vector<unsigned>::iterator IndexIt;	


class VortexBlobs
{
public:
	VortexBlobs();
	void resize(unsigned size);	
	unsigned size();
	void print_location();
	void print_velocity();
    	void print_circulation();
    	void print_cutoff();
    	void print_vorticity();
    	void print_all();
    	double totalcirc();

public: 
    
    std::vector<unsigned> ID; // Vortex blob ID
    std::vector<double> x; //x-coordinate   
    std::vector<double> z; //z-coordinate
    std::vector<double> circ; //circulation
    std::vector<double> sigma; //vortex cut-off
    std::vector<double> u; // local x-velocity
    std::vector<double> w; // local z-velocity
    std::vector<double> uvs; 
    std::vector<double> wvs;
    std::vector<double> omega; //Local vorticity
};

class VortexSheet
{
public:
    VortexSheet();
    void resize(unsigned size);
    unsigned size();
    void print_collocation();
    void print_unit_vectors();
    void print_gamma();

public:
    std::vector<double> gamma; // surface vorticity
    std::vector<double> x;    // x-coordinate point
    std::vector<double> z;    // z-coordinate point
    std::vector<double> xc;  // x-collocation point
    std::vector<double> zc;  // z-collocation point
    std::vector<double> theta; //angle
    std::vector<double> ds;  // panel
    std::vector<double> enx; // x-coor normal unit vector  
    std::vector<double> enz; // z-coor normal unit vector
    std::vector<double> etx; // x-coor tangential unit vector
    std::vector<double> etz; // z-coor tangential unit vector 
};

class Body
{
public:
    Body();
    void resize(unsigned size);
    unsigned size();
    void print_location();

public:
    std::vector<double> x; // x-coordinate
    std::vector<double> z; // z-coordinate
};
#endif

