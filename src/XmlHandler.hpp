// Wrapper for pugixml to give us nice error checking etc
#ifndef XMLHANDLER_H
#define XMLHANDLER_H

#include "pugixml.hpp"
#include "pugiconfig.hpp"
#include <string>
#include <vector>
#include <unordered_map>

///\file XmlHandler.hpp Wrapper for pugixml to provide error checking and nice accessors

class XmlHandler
{
	private:
	/// The parsed xml document
	pugi::xml_document m_xml;

	/// Map values
	struct opts {
		std::string level;     ///< Highest level element
		std::string parameter; ///< Parameter
		std::string attr;      ///< Attribute type
		std::string info;      ///< Description of element
		std::string reqs;      ///< Any required options. None otherwise
	};

	/// A map of the possible elements and possible values
	/** Keys are string of "level_parameter" as defined in the getXXAttribute
	 * methods below.
	 *
	 * An exception is throw if an element is provided that isn't in the map.
	 * This means that we should end up with a self documenting input file */
	std::unordered_map<std::string, opts> m_option_map;

	/// Returns the node given for the parent element level and child element
	/// parameter.
	/** This assumes one layer of nesting in the xml. Should be sufficient for
	 * the moment and having everything in here make changing things easier at a
	 * later date
	 *
	 * If the requested parameter does not exist an exception is thrown */
	pugi::xml_node getNode(const char *level, const char *parameter) const;

	/// Check that the provided string option is valid
	/** Throws exception if invalid option given with a hint */
	void check_string(const char* level, const char *parameter, std::string attr) const;

	/// Check that the provided value option is valid
	/** Throws excetion if invalid option given with a hint */
	void check_value(const char* level, const char *parameter, double attr) const;

	/// Form a key from the map given tho level and the parameter
	std::string form_key(const char* level, const char *parameter) const;

	/// Builds the option map - called by the constructor
	/** Just a long series of adding to the map statements */
	void buildOptionMap();

	public:
	/// Constructor to just build the option map
	XmlHandler();

	/// Constructor parses the xml_file and throws if there is a problem
	XmlHandler(std::string xml_file);

	/// Returns a string attribute
	std::string getStringAttribute(const char *level,
	                               const char *parameter) const;

	/// Returns a vector of the list of values in the string.
	/** Can be either space or comma separated */
	std::vector<double> getList(const char *level,
	                                  const char *parameter) const;

	/// Returns a double attribute
	double getValueAttribute(const char *level, const char *parameter) const;

	/// Returns and integer attribute
	int getIntAttribute(const char *level, const char *parameter) const;

	/// Returns a bool attribute
	bool getBoolAttribute(const char *level, const char *parameter) const;

	/// Saves the xml document to filename
	void save(std::string filename) const;

	/// Writes an example xml file using the option map
	void writeExample(std::string file) const;
};
#endif // XMLHANDLER_H

