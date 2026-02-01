#include "XmlHandler.hpp"
#include "Exceptions.hpp"

#include <algorithm>
#include <iostream>
#include <sstream>

XmlHandler::XmlHandler()
{
    buildOptionMap();
}

XmlHandler::XmlHandler(std::string xml_file)
{
    pugi::xml_parse_result result = m_xml.load_file(xml_file.c_str());

    if (!result) {
        throw dvm::XMLParseException::cannotLoad(xml_file, result.description());
    }

    std::cout << "Loaded xml input: " << xml_file << "\n";
    buildOptionMap();
}

pugi::xml_node XmlHandler::getNode(const char *level,
                                   const char *parameter) const
{
    pugi::xml_node node = m_xml.child("DVMpp").child(level).child(parameter);
    if (!node) {
        // The node hasn't been found. There are two options:
        // (1) user has forgotten it -> provide a hint
        // (2) programmer is trying to add a new option -> must add option to map
        auto key = form_key(level, parameter);
        auto val = m_option_map.find(key);
        if (val != m_option_map.end()) {
            // case (1)
            std::ostringstream err;
            err << "xml error: <" << level << "><" << parameter
                << "> missing.\nAttribute: " << val->second.attr
                << "\nInfo: " << val->second.info
                << "\nOptions: " << val->second.reqs;
            throw dvm::XMLParseException(err.str());
        } else {
            // case (2)
            std::ostringstream err;
            err << "xml error: Are you trying to add an option without "
                   "updating XmlHandler::m_option_map?\nRequested option: <"
                << level << "><" << parameter << ">";
            throw dvm::XMLParseException(err.str());
        }
    }

    // Make sure the node is valid by checking the option map
    auto e = m_option_map.find(form_key(level, parameter));
    if (e == m_option_map.end()) {
        throw dvm::XMLParseException("<" + std::string(level) + "><" +
                                     std::string(parameter) + "> not found in the option map.");
    }

    return node;
}

std::string XmlHandler::getStringAttribute(const char *level,
                                           const char *parameter, bool checksep) const
{
    std::string attr = getNode(level, parameter).attribute("string").as_string();
    check_string(level, parameter, attr);

    if (checksep) {
        if (attr.back() != '/') {
            attr += '/';
        }
    }

    return attr;
}

std::vector<double> XmlHandler::getList(const char *level,
                                        const char *parameter) const
{
    std::string vals = getNode(level, parameter).attribute("string").as_string();
    check_string(level, parameter, vals);

    double val;
    std::vector<double> vec;
    std::istringstream s(vals);
    while (s >> val) {
        vec.push_back(val);

        if (s.peek() == ' ' || s.peek() == ',') {
            s.ignore();
        }
    }

    return vec;
}

double XmlHandler::getValueAttribute(const char *level,
                                     const char *parameter) const
{
    auto val = getNode(level, parameter).attribute("val").as_double();
    check_value(level, parameter, val);

    return val;
}

int XmlHandler::getIntAttribute(const char *level, const char *parameter) const
{
    return getNode(level, parameter).attribute("val").as_int();
}

bool XmlHandler::getBoolAttribute(const char *level,
                                  const char *parameter) const
{
    // pugi bool conversion will return false for attribute not starting with
    // '1', 't', 'T', 'y', 'Y'. This is probably safe for the moment...
    return getNode(level, parameter).attribute("string").as_bool();
}

void XmlHandler::save(std::string filename) const
{
    m_xml.save_file(filename.c_str());
}

void XmlHandler::check_string(const char* level, const char *parameter, std::string &attr) const
{
    // Get the required values, get node has already ensured we have got something valid.
    auto req = m_option_map.find(form_key(level, parameter))->second.reqs;

    // Make sure we have a required value
    if (req == "none") {
        return;
    }

    if (req.find(attr) == std::string::npos) {
        // can't find the given attr in the required string
        throw dvm::XMLParseException::invalidOption(level, parameter, req);
    }
}

void XmlHandler::check_value(const char* level, const char *parameter, double attr) const
{
    // Get the required values, get node has already ensured we have got something valid.
    auto req = m_option_map.find(form_key(level, parameter))->second.reqs;

    if (req == "none") {
        // do nothing - there is no required value
    } else if (req == ">0") {
        if (attr <= 0) {
            throw dvm::XMLParseException::invalidOption(level, parameter, req);
        }
    } else {
        throw dvm::XMLParseException::invalidOption(level, parameter, req);
    }
}

std::string XmlHandler::form_key(const char *level, const char *parameter) const
{
    return std::string(level) + "_" + std::string(parameter);
}

void XmlHandler::buildOptionMap()
{
    m_option_map.insert({"algorithms_surface_crossing", {"algorithms", "surface_crossing", "string", "surface-crossing algorithm", "DELETE ABSORB REFLECT"}});
    m_option_map.insert({"algorithms_velocity_kernel", {"algorithms", "velocity_kernel", "string", "velocity computation kernel", "direct fmm"}});
    m_option_map.insert({"algorithms_fmm_expansion_order", {"algorithms", "fmm_expansion_order", "val", "FMM expansion order", ">0"}});
    m_option_map.insert({"algorithms_fmm_theta", {"algorithms", "fmm_theta", "val", "FMM opening angle criterion", "none"}});

    m_option_map.insert({"constants_density", {"constants", "density", "val", "specific gravity of fluid", ">0"}});
    m_option_map.insert({"constants_nu", {"constants", "nu", "val", "kinematic viscosity [m^2 / s]", ">0"}});
    m_option_map.insert({"constants_max_NumPanelVort", {"constants", "max_NumPanelVort", "val", "maximum number of panel vortices", ">0"}});
    m_option_map.insert({"constants_cutoff_exp", {"constants", "cutoff_exp", "val", "cutoff exp", "none"}});
    m_option_map.insert({"constants_seed", {"constants", "seed", "val", "seed for rng. <0 random, >=0 seed", "none"}});

    m_option_map.insert({"io_input_dir", {"io", "input_dir", "string", "input directory", "none"}});
    m_option_map.insert({"io_output_dir", {"io", "output_dir", "string", "output directory", "none"}});
    m_option_map.insert({"io_domain_file", {"io", "domain_file", "string", "body geometry file", "none"}});

    m_option_map.insert({"flow_ux", {"flow", "ux", "val", "flow in the x-direction", "none"}});
    m_option_map.insert({"flow_uz", {"flow", "uz", "val", "flow in the z-direction", "none"}});

    m_option_map.insert({"probe_x", {"probe", "x", "val", "x-coordinate of the probe", "none"}});
    m_option_map.insert({"probe_z", {"probe", "z", "val", "z-coordinate of the probe", "none"}});

    m_option_map.insert({"time_scheme", {"time", "scheme", "string", "time-stepping scheme", "euler RK3"}});
    m_option_map.insert({"time_dt", {"time", "dt", "val", "time step [s]", ">0"}});
    m_option_map.insert({"time_steps", {"time", "steps", "val", "total number of timesteps", ">0"}});
}

void XmlHandler::writeExample(std::string file) const
{
    // create the xml document and add the top element
    pugi::xml_document xml_out;
    pugi::xml_node DVM = xml_out.append_child("DVMpp");

    // Go through the map to build the document tree
    for (auto e = m_option_map.begin(); e != m_option_map.end(); ++e) {
        auto l = e->second.level.c_str();
        auto p = e->second.parameter.c_str();
        auto attr = e->second.attr.c_str();
        auto info = e->second.info;
        auto req = e->second.reqs;

        // Check if the level already exists.
        pugi::xml_node node = DVM.child(l);
        if (!node) {
            // Create the level
            node = DVM.append_child(l);
        }

        // Add the parameter node and set the attribute
        auto param_node = node.append_child(p);
        param_node.append_attribute(attr) =
            (info + ". Options: " + req).c_str();
    }

    if (file == "-") {
        xml_out.save(std::cout);
    } else {
        std::cout << "Saving example xml to " << file << std::endl;
        xml_out.save_file(file.c_str());
    }
}
