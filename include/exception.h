#ifndef VEGAS_EXCEPTION_H
#define VEGAS_EXCEPTION_H

#include <stdexcept>
#include <string>

namespace vegas {

/// @brief Base exception class for VEGAS errors
class VEGASException : public std::runtime_error {
public:
    explicit VEGASException(const std::string& message)
        : std::runtime_error(message) {}
};

/// @brief Exception for file I/O errors (cannot open, read, write, etc.)
class FileIOException : public VEGASException {
public:
    explicit FileIOException(const std::string& message)
        : VEGASException("File I/O error: " + message) {}
};

/// @brief Exception for invalid input data (malformed files, wrong format, etc.)
class InvalidInputException : public VEGASException {
public:
    explicit InvalidInputException(const std::string& message)
        : VEGASException("Invalid input: " + message) {}
};

/// @brief Exception for numeric conversion errors (parsing numbers, out of range, etc.)
class NumericConversionException : public VEGASException {
public:
    explicit NumericConversionException(const std::string& message)
        : VEGASException("Numeric conversion error: " + message) {}
};

/// @brief Exception for configuration errors (JSON parsing, inconsistent parameters, etc.)
class ConfigurationException : public VEGASException {
public:
    explicit ConfigurationException(const std::string& message)
        : VEGASException("Configuration error: " + message) {}
};

/// @brief Exception for simulation errors (invalid spin states, energy calculation, etc.)
class SimulationException : public VEGASException {
public:
    explicit SimulationException(const std::string& message)
        : VEGASException("Simulation error: " + message) {}
};

/// @brief Exception for HDF5-related errors (file creation, dataset writing, etc.)
class HDF5Exception : public VEGASException {
public:
    explicit HDF5Exception(const std::string& message)
        : VEGASException("HDF5 error: " + message) {}
};

} // namespace vegas

#endif // VEGAS_EXCEPTION_H