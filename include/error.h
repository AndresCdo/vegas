#ifndef ERROR_H
#define ERROR_H

#include "exception.h"
#include <string>

/// @brief Throw an exception with error message (replaces old exit behavior)
/// @param message Error message to display
inline void EXIT(std::string message)
{
    throw vegas::VEGASException(message);
}

#endif // ERROR_H