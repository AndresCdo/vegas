#ifndef VEGAS_CONFIG_PARSER_H
#define VEGAS_CONFIG_PARSER_H

#include "params.h"
#include "json/json.h"
#include <string>
#include <vector>
#include <optional>

namespace vegas {

struct SimulationConfig {
    std::string sampleFile;
    std::string outputFile;
    std::string initialStateFile;
    std::vector<std::string> anisotropyFiles;
    std::vector<Real> temperatures;
    std::vector<Real> fields;
    Index mcs = DEFAULT_MCS;
    Index seed = 0;
    Real kb = 1.0;
    bool hasInitialState = false;
    bool hasAnisotropy = false;
    
    void applyOverrides(
        const std::optional<Index>& mcsOverride,
        const std::optional<Index>& seedOverride,
        const std::optional<Real>& kbOverride,
        const std::optional<std::string>& outputOverride,
        const std::optional<std::string>& sampleOverride,
        const std::optional<std::string>& initialStateOverride
    );
};

class ConfigParser {
public:
    static SimulationConfig parse(const std::string& jsonFile);
    
private:
    static std::vector<Real> parseValueRange(const Json::Value& json, 
                                              const std::string& key,
                                              bool& isUnique,
                                              Real& singleValue);
    static void validate(const SimulationConfig& config);
    static std::vector<Real> expandSingleValue(Real value, size_t count);
};

} // namespace vegas

#endif // VEGAS_CONFIG_PARSER_H
