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

    
    DVMBase bl;
    
    //Parameters
    bl.m_maxGamma=0.05;
    bl.m_dt=0.5;
    bl.m_rho=1.0;
    double T=100;
    bl.m_steps=T/bl.m_dt;
    bl.m_nu=0.000001;
    bl.m_Ux=0.01;
    bl.m_Uz=0;
    bl.m_body;

    bl.m_probe.resize(1);
    // Define probe point
    bl.m_probe.x[0]=20;
    bl.m_probe.z[0]=1;

    bl.read_input_coord();

    bl.form_vortex_sheet(); 
    bl.m_Gamma_abs.resize(bl.m_vortsheet.size());
    bl.compute_influence_matrix();
 
    //cin.get();

    bl.init_outputs(); 
    
    clock_t start;
    double cpu_time;

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
    
    
}catch (char* str){
    std::cout << "Exception thrown: " << str << std::endl;
}catch (std::string str){
    std::cout << "Exception thrown: " << str << std::endl;
}catch ( ...) {
    std::cout << "Unhandled exception type" << std::endl;
}
    
return 0;
    
    
}


