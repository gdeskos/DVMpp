// The main application file
#include "BaseTypes.hpp"
#include "DVMBase.hpp"
#include "XmlHandler.hpp"
#include "pugiconfig.hpp"
#include "pugixml.hpp"
#include <boost/program_options.hpp>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

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
		auto time = std::time(nullptr);
		auto loctime = *std::localtime(&time);

		std::ostringstream oss;
		oss << std::put_time(&loctime, "%Y_%m_%d_%H_%M_%S");
		std::string stamp = oss.str();
		std::cout << "File timestamp is " << stamp << std::endl;

		DVMBase dvm(xml, stamp);

		auto outdir = xml.getStringAttribute("io", "output_dir");
		outdir = (outdir.back() == '/') ? outdir : outdir + '/';
		xml.save(outdir + stamp + "_xml_in.xml");
		auto start = std::chrono::system_clock::now();

		/// Solve the problem
		dvm.solve();

		// output the cpu time on screen
		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed = end - start;
		std::cout << "\nRuntime : " << elapsed.count() << "s\n";

		std::cout << "Timestamp : " << stamp << "\n";

		std::cout << "Completed DVM computations succesfully" << std::endl;

	} catch (char *str) {
		std::cout << "Exception thrown: " << str << std::endl;
	} catch (std::string str) {
		std::cout << "Exception thrown: " << str << std::endl;
	} catch (...) {
		std::cout << "Unhandled exception type" << std::endl;
	}

	return 0;
}
