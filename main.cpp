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

    
    DVMBase dvm;
    
    dvm.init(xml_doc);
    dvm.read_input_coord();
    dvm.init_outputs();
    
    
    dvm.form_vortex_sheet();
    dvm.compute_influence_matrix();
    
    
    
    clock_t start;
    double cpu_time;

    start = clock();

    // Start the time loop;
    for(unsigned j=1;j<=dvm.get_steps();j++){
        
        dvm.increment_step();
        
        if(dvm.get_vs_size()==0){ // Is this the correct get_size command?
            dvm.solvevortexsheet();
            dvm.save_vort();
        }else{
        
            // Inviscid Substep
            dvm.solvevortexsheet();
            dvm.save_vort();
  
            //bl.compute_loads();
            
            dvm.convect(1); // first order second order scheme
            dvm.diffrw();
        }
        
        
        // Viscous Substep
        dvm.diffuse_vs_rw(); // a number of question here - not entirely clear
       
    
        // Housekeeping
        dvm.reflect();
   
        
        // Compute certain values
        //bl.probe_velocities();

        // Output
        dvm.write_outputs();
        
        
         
        // Screen output
        cout<<"Simulation time          = " << dvm.get_time() << "\tStep " << j << "/" << dvm.get_steps() << std::endl;
        cout<<"Number of vortex blobs   = " << dvm.get_size() << std::endl;
    }
       
    // output the cpu time on screen 
    cpu_time = (clock() - start) / (double) CLOCKS_PER_SEC;
    cout << "Used CPU time is : "<< cpu_time << endl;
    
    
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


