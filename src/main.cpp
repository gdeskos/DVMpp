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

		clock_t start;
		double cpu_time;

		start = clock();

		/// Solve the problem
		dvm.solve();

		// output the cpu time on screen
		cpu_time = (clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << "Used CPU time is : " << cpu_time << std::endl;

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
