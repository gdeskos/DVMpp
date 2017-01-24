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
#include <boost/program_options.hpp>

int main(int argc, char *argv[])
{

	try {
		// Deal with command line options
		namespace po = boost::program_options;
		po::options_description desc("DVM++ commandline options");
		desc.add_options()("help,h", "prints this help message")(
		    "input-file,f", po::value<std::string>(), "input xml file")(
		    "example-xml", "prints an example xml to stdout or file");

		po::positional_options_description p;
		p.add("input-file", 1);

		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv)
		              .options(desc)
		              .positional(p)
		              .run(),
		          vm);
		po::notify(vm);

		if (vm.empty()) {
			std::cout << desc << "\n";
			return 1;
		}

		if (vm.count("example-xml")) {
			XmlHandler xml;

			// The input file option will capture the positional argument for
			// the name of the example xml so we can just use it!
			std::string file = "-";
			if (vm.count("input-file")) {
				file = vm["input-file"].as<std::string>();
			}
			xml.writeExample(file);

			return 1;
		}

		XmlHandler xml(vm["input-file"].as<std::string>());

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
