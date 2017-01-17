// The main application file
#include "BaseTypes.hpp"
#include "DVMBase.hpp"
#include "XmlHandler.hpp"
#include "pugiconfig.hpp"
#include "pugixml.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <time.h>

using namespace std;

int main(int argc, char *argv[])
{

	try {

		if (argc != 2) {
			throw std::string("Must supply input XML file name");
		}

		std::string input_xml_file(argv[1]);

		XmlHandler xml(input_xml_file);

		std::string experiment_name, domain_type;
		// FileNames file_names;

		std::cout
		    << "=============================================================="
		    << std::endl;
		std::cout << "Running the ICL DVM Test code Cpp version." << std::endl;
		std::cout
		    << "=============================================================="
		    << std::endl;

		std::cout << "Successfully loaded XML input file" << std::endl;

		// Deal with time - this is horrible!
		time_t rawtime;
		struct tm *ptm;
		time(&rawtime);
		ptm = gmtime(&rawtime);

		std::ostringstream os;
		os << ptm->tm_year + 1900 << "_";
		if (ptm->tm_mon + 1 < 10) {
			os << "0";
		}
		os << ptm->tm_mon + 1 << "_";
		if (ptm->tm_mday < 10) {
			os << "0";
		}
		os << ptm->tm_mday << "_";
		if (ptm->tm_hour + 1 < 10) {
			os << "0";
		}
		os << ptm->tm_hour + 1 << "_";
		if (ptm->tm_min < 10) {
			os << "0";
		}
		os << ptm->tm_min << "_";
		if (ptm->tm_sec < 10) {
			os << "0";
		}
		os << ptm->tm_sec;

		std::string stamp = os.str();

		std::cout << "File timestamp is " << os.str() << std::endl;

		DVMBase dvm(xml);

		dvm.init(xml, stamp);
		dvm.read_input_coord();
		dvm.init_outputs();

		dvm.form_vortex_sheet();
		dvm.compute_influence_matrix();

		clock_t start;
		double cpu_time;

		start = clock();

		// Start the time loop;
		for (unsigned j = 1; j <= dvm.get_steps(); j++) {

			dvm.increment_step();

			if (dvm.get_vs_size()
			    == 0) { // Is this the correct get_size command?
				dvm.solvevortexsheet();
				dvm.save_vort();
			} else {

				// Inviscid Substep
				dvm.solvevortexsheet();

				dvm.compute_loads();

				dvm.save_vort();

				dvm.convect(1); // first order second order scheme
				dvm.diffrw();
			}

			// Viscous Substep
			dvm.diffuse_vs_rw(); // a number of question here - not entirely
			                     // clear

			// Housekeeping
			dvm.reflect();

			// Compute certain values
			// bl.probe_velocities();

			// Output
			dvm.write_outputs();

			// Screen output
			cout << "Simulation time          = " << dvm.get_time() << "\tStep "
			     << j << "/" << dvm.get_steps() << std::endl;
			cout << "Number of vortex blobs   = " << dvm.get_size()
			     << std::endl;
		}

		// output the cpu time on screen
		cpu_time = (clock() - start) / (double)CLOCKS_PER_SEC;
		cout << "Used CPU time is : " << cpu_time << endl;

		auto outdir = xml.getStringAttribute("io", "output_dir");
		xml.save(outdir + stamp + "_xml_in.xml");

	} catch (char *str) {
		std::cout << "Exception thrown: " << str << std::endl;
	} catch (std::string str) {
		std::cout << "Exception thrown: " << str << std::endl;
	} catch (...) {
		std::cout << "Unhandled exception type" << std::endl;
	}

	std::cout << "Completed DVM computations succesfully" << std::endl;

	return 0;
}
