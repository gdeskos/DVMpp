#ifndef DVM_EXCEPTIONS_HPP
#define DVM_EXCEPTIONS_HPP

/** \file Exceptions.hpp
 * \brief Exception hierarchy for DVMpp */

#include <stdexcept>
#include <string>

namespace dvm {

/// Base exception class for all DVMpp errors
class DVMException : public std::runtime_error {
public:
    explicit DVMException(const std::string& message)
        : std::runtime_error(message) {}

    explicit DVMException(const char* message)
        : std::runtime_error(message) {}
};

/// Exception for size mismatches between data structures
class SizeMismatchException : public DVMException {
public:
    explicit SizeMismatchException(const std::string& context)
        : DVMException("Size mismatch in " + context) {}
};

/// Exception for file I/O errors
class FileIOException : public DVMException {
public:
    explicit FileIOException(const std::string& message)
        : DVMException(message) {}

    static FileIOException cannotOpen(const std::string& filename) {
        return FileIOException("Cannot open file: " + filename);
    }

    static FileIOException cannotWrite(const std::string& filename) {
        return FileIOException("Cannot write to file: " + filename);
    }

    static FileIOException cannotRead(const std::string& filename) {
        return FileIOException("Cannot read file: " + filename);
    }
};

/// Exception for XML parsing errors
class XMLParseException : public DVMException {
public:
    explicit XMLParseException(const std::string& message)
        : DVMException(message) {}

    static XMLParseException cannotLoad(const std::string& filename, const std::string& reason) {
        return XMLParseException("Cannot load XML file '" + filename + "': " + reason);
    }

    static XMLParseException missingParameter(const std::string& level, const std::string& parameter) {
        return XMLParseException("Missing XML parameter: <" + level + "><" + parameter + ">");
    }

    static XMLParseException invalidOption(const std::string& level, const std::string& parameter,
                                           const std::string& options) {
        return XMLParseException("Invalid option for <" + level + "><" + parameter + ">. "
                                "Valid options: " + options);
    }
};

/// Exception for invalid algorithm state or parameters
class AlgorithmException : public DVMException {
public:
    explicit AlgorithmException(const std::string& message)
        : DVMException(message) {}
};

/// Exception for invalid vortex operations
class VortexException : public DVMException {
public:
    explicit VortexException(const std::string& message)
        : DVMException(message) {}
};

} // namespace dvm

#endif // DVM_EXCEPTIONS_HPP
