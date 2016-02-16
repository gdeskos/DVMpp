// The main application file
#include <iostream>
#include <time.h>
#include "BaseTypes.hpp"
#include "DVMBase.hpp"
#include "pugiconfig.hpp"
#include "pugixml.hpp"
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
    
try{
    
    if(argc != 2){ throw std::string("Must supply input XML file name");}
    
    std::string input_xml_file(argv[1]);
    
    pugi::xml_document xml_doc;
    pugi::xml_parse_result result = xml_doc.load_file(input_xml_file.c_str());
    
    std::string experiment_name, domain_type;
    //FileNames file_names;
    
    if(!result){
        std::ostringstream os;
        os << "An error occured whilst reading the input XML file: " << result.description();
        throw std::string(os.str());
    }
    
    std::cout << "==============================================================" << std::endl;
    std::cout << "Running the ICL DVM Test code Cpp version."  << std::endl;
    std::cout << "==============================================================" << std::endl;

    std::cout << "Successfully loaded XML input file" << std::endl;

//    std::string in_dir = xml_doc.child("IC2DDVM").child("io").child("input_dir").attribute("string").value();
//    std::string out_dir = xml_doc.child("IC2DDVM").child("io").child("output_dir").attribute("stringX").value();

    
    double rho = atof(xml_doc.child("IC2DDVM").child("constants").child("density").attribute("val").value());
    if(rho <=0){ throw std::string("<IC2DDVM><constants><density> must be larger than zero");}
   
    double nu = atof(xml_doc.child("IC2DDVM").child("constants").child("nu").attribute("val").value());
    if(nu <=0){ throw std::string("<IC2DDVM><constants><nu> must be larger than zero");}
    
    double max_gamma = atof(xml_doc.child("IC2DDVM").child("constants").child("max_gamma").attribute("val").value());
    if(max_gamma <=0){ throw std::string("<IC2DDVM><constants><max_gamma> must be larger than zero");}
    
    
    double dt = atof(xml_doc.child("IC2DDVM").child("time").child("dt").attribute("val").value());
    if(dt <=0){ throw std::string("<IC2DDVM><time><dt> must be larger than zero");}
    
    unsigned steps = atof(xml_doc.child("IC2DDVM").child("time").child("steps").attribute("val").value());
    if(steps <=0){ throw std::string("<IC2DDVM><time><steps> must be larger than zero");}
  
    unsigned ux = atof(xml_doc.child("IC2DDVM").child("flow").child("ux").attribute("val").value());
    unsigned uz = atof(xml_doc.child("IC2DDVM").child("flow").child("uz").attribute("val").value());
    
    // Generalise this to arrays
    unsigned probe_x = atof(xml_doc.child("IC2DDVM").child("probe").child("x").attribute("val").value());
    unsigned probe_z = atof(xml_doc.child("IC2DDVM").child("probe").child("z").attribute("val").value());
    
    
    DVMBase bl;
    
    //Parameters
    bl.m_maxGamma=max_gamma;
    bl.m_dt=dt;
    bl.m_rho=rho;
    bl.m_steps=steps;
    
    bl.m_nu=nu;
    bl.m_Ux=ux;
    bl.m_Uz=uz;
   // bl.m_body; ?? What is this command doing here?

    bl.m_probe.resize(1);
 
    // Define probe point
    bl.m_probe.x[0]=probe_x;
    bl.m_probe.z[0]=probe_z;

    /*
    
    bl.read_input_coord();

    bl.form_vortex_sheet(); 
    bl.m_Gamma_abs.resize(bl.m_vortsheet.size());
    bl.compute_influence_matrix();
 
    bl.init_outputs(); 
    
    clock_t start;
    double cpu_time;

    */
    /*
    
    start = clock();
    // Start the time loop;
    for(unsigned j=1;j<=bl.m_steps;j++)
    {
        bl.m_step=j; 
        bl.m_time=bl.m_dt*j; 

        if(bl.m_vortex.size()==0){
            
            bl.solvevortexsheet();
            bl.save_vort();
        
        }else{
        
        // Inviscid Substep
        bl.solvevortexsheet();
        //bl.compute_loads();
        bl.save_vort();
  
        bl.convect(1,bl.m_dt);
        bl.diffrw(bl.m_nu,bl.m_dt); 

        }
        

        // Viscous Substep
        bl.diffuse_vs_rw();
       
        // Housekeeping
        bl.reflect();
   
        // Compute certain values
        bl.probe_velocities();

        // Output
        bl.write_outputs(); 
        
        // Screen output
        cout<<"Simulation message : time = "<<bl.m_time<<std::endl;
        cout<<"Number of vortex blobs : NVB = "<<bl.m_vortex.size()<<std::endl;
    }
       
    // output the cpu time on screen 
    cpu_time = (clock() - start) / (double) CLOCKS_PER_SEC;
    cout<<"The CPU time is : "<< cpu_time<<endl;
    
    */
}catch (char* str){
    std::cout << "Exception thrown: " << str << std::endl;
}catch (std::string str){
    std::cout << "Exception thrown: " << str << std::endl;
}catch ( ...) {
    std::cout << "Unhandled exception type" << std::endl;
}
    
    std::cout << "Completed DVM computations succesfully" << std::endl;
    
return 0;
    
    
}


