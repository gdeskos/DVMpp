#include "Body.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

Body::Body()
{
	// Don't put anything in here - it is not called
}

Body::Body(XmlHandler &xml)
{
	auto indir = xml.getStringAttribute("io", "input_dir");
	// Make sure we have trailing separator
	if (indir.back() != '/') {
		indir += '/';
	}

	std::string domain = indir + xml.getStringAttribute("io", "domain_file");

	read_input_coord(domain);
}

void Body::resize(unsigned size)
{
	x.resize(size);
	z.resize(size);
}

void Body::read_input_coord(std::string file)
{
	// Read the body coordinates
	std::string line;
	std::ifstream coor_file(file.c_str());
	double tmp1, tmp2;

	if (coor_file.is_open()) {
		while (coor_file.good()) {
			std::getline(coor_file, line);
			std::istringstream buffer(line);
			buffer >> tmp1 >> tmp2;
			x.push_back(tmp1);
			z.push_back(tmp2);
		}
		x.pop_back();
		z.pop_back();
		std::cout << "Succesfully loaded coordinate file with " << size()
		          << " points." << std::endl;
		coor_file.close();
	} else {
		std::string error_msg;
		error_msg = "Unable to open coordinate file from " + file;
		throw error_msg;
	}

	print_location();
}

unsigned Body::size()
{
	if (x.size() != z.size()) {
		throw std::string("Size mismatch in Body");
	}
	return x.size();
}

void Body::print_location()
{
	for (unsigned i = 0; i < size(); i++) {
		std::cout << " Body coor (x,z) = (" << x[i] << "," << z[i] << ")"
		          << std::endl;
	}
}
