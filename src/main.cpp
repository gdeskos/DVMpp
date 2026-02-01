// The main application file
#include "BaseTypes.hpp"
#include "DVMBase.hpp"
#include "Exceptions.hpp"
#include "XmlHandler.hpp"
#include "pugiconfig.hpp"
#include "pugixml.hpp"
#include <boost/program_options.hpp>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace fs = std::filesystem;

int main(int argc, char *argv[])
{
    try {
        // Deal with command line options
        namespace po = boost::program_options;
        po::options_description desc("DVMpp commandline options");
        desc.add_options()
            ("help,h", "prints this help message")
            ("input-file,f", po::value<std::string>(), "input xml file")
            ("example-xml", "prints an example xml to stdout or file");

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

            return 0;
        }

        std::cout
            << "=============================================================="
            << std::endl;
        std::cout << "Running DVMpp - 2D Discrete Vortex Method" << std::endl;
        std::cout
            << "=============================================================="
            << std::endl;

        XmlHandler xml(vm["input-file"].as<std::string>());

        // Deal with time - create timestamp for output files
        auto time = std::time(nullptr);
        auto loctime = *std::localtime(&time);

        std::ostringstream oss;
        oss << std::put_time(&loctime, "%Y_%m_%d_%H_%M_%S");
        std::string stamp = oss.str();
        std::cout << "File timestamp is " << stamp << std::endl;

        DVMBase dvm(xml, stamp);

        auto outdir = xml.getStringAttribute("io", "output_dir");
        fs::path outpath(outdir);
        if (outpath.string().back() != '/') {
            outpath = fs::path(outdir + "/");
        }
        xml.save((outpath / (stamp + "_xml_in.xml")).string());

        auto start = std::chrono::system_clock::now();

        // Solve the problem
        dvm.solve();

        std::cout << "\nTimestamp : " << stamp << "\n";

        // output the cpu time on screen
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Runtime : " << elapsed.count() << "s\n";

        std::cout << "Completed DVM computations successfully" << std::endl;

    } catch (const dvm::DVMException& e) {
        std::cerr << "DVMpp Error: " << e.what() << std::endl;
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown exception occurred" << std::endl;
        return 1;
    }

    return 0;
}
