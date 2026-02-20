#include "../include/config_parser.h"
#include "../include/exception.h"

#include <fstream>
#include <ctime>

namespace vegas {

SimulationConfig ConfigParser::parse(const std::string& jsonFile) {
    SimulationConfig config;
    
    // Parse JSON file
    Json::Value root;
    Json::Reader reader;
    std::ifstream file(jsonFile, std::ifstream::binary);
    
    if (!reader.parse(file, root, false)) {
        throw vegas::FileIOException("Cannot parse configuration file: " + jsonFile);
    }
    
    // Extract basic parameters
    config.sampleFile = root["sample"].asString();
    config.mcs = root.get("mcs", DEFAULT_MCS).asInt();
    config.kb = root.get("kb", 1.0).asDouble();
    config.outputFile = root.get("out", config.sampleFile + ".h5").asString();
    config.seed = root.get("seed", static_cast<Index>(std::time(nullptr))).asUInt();
    
    // Parse temperature
    bool uniqueT = false;
    Real singleT = 0.0;
    config.temperatures = parseValueRange(root, "temperature", uniqueT, singleT);
    
    // Parse field
    bool uniqueH = false;
    Real singleH = 0.0;
    config.fields = parseValueRange(root, "field", uniqueH, singleH);
    
    // Handle single values vs arrays
    if (uniqueT && uniqueH) {
        // Both are single values - already in vectors
    } else if (uniqueT && !uniqueH) {
        // Single temperature, multiple fields
        config.temperatures = expandSingleValue(singleT, config.fields.size());
    } else if (!uniqueT && uniqueH) {
        // Multiple temperatures, single field
        config.fields = expandSingleValue(singleH, config.temperatures.size());
    }
    
    // Parse initial state
    config.initialStateFile = root.get("initialstate", "").asString();
    config.hasInitialState = root.isMember("initialstate");
    
    // Parse anisotropy files
    config.hasAnisotropy = root.isMember("anisotropy");
    if (config.hasAnisotropy) {
        if (root.get("anisotropy", 0.0).type() == 4) {
            // Single string
            config.anisotropyFiles.push_back(root.get("anisotropy", "").asString());
        } else if (root.get("anisotropy", 0.0).type() == 6) {
            // Array of strings
            const Json::Value files = root["anisotropy"];
            for (Index i = 0; i < files.size(); ++i) {
                config.anisotropyFiles.push_back(files[i].asString());
            }
        }
    }
    
    validate(config);
    return config;
}

std::vector<Real> ConfigParser::parseValueRange(const Json::Value& root,
                                                  const std::string& key,
                                                  bool& isUnique,
                                                  Real& singleValue) {
    std::vector<Real> values;
    isUnique = true;
    
    if (root.get(key, 0.0).type() == 1) {
        // Integer
        singleValue = static_cast<Real>(root.get(key, 0.0).asInt());
        values.push_back(singleValue);
    } else if (root.get(key, 0.0).type() == 2) {
        // Unsigned integer
        singleValue = static_cast<Real>(root.get(key, 0.0).asUInt());
        values.push_back(singleValue);
    } else if (root.get(key, 0.0).type() == 3) {
        // Double
        singleValue = root.get(key, 0.0).asDouble();
        values.push_back(singleValue);
    } else if (root.get(key, 0.0).type() == 6) {
        // Array of values
        isUnique = false;
        const Json::Value arr = root[key];
        for (Index i = 0; i < arr.size(); ++i) {
            values.push_back(arr[i].asDouble());
        }
    } else if (root.get(key, 0.0).type() == 7) {
        // Dictionary with range specification
        isUnique = false;
        const Json::Value range = root[key];
        
        Real start = range.get("start", key == "temperature" ? DEFAULT_TEMP_START : DEFAULT_FIELD_START).asDouble();
        Real finalVal = range.get("final", key == "temperature" ? DEFAULT_TEMP_FINAL : DEFAULT_FIELD_FINAL).asDouble();
        bool cycle = range.get("cycle", false).asBool();
        
        Index points = range.get("points", 5).asUInt();
        Real delta = (finalVal - start) / static_cast<Real>(points - 1);
        
        if (range.isMember("points") && range.isMember("delta")) {
            throw vegas::ConfigurationException(key + " section has both 'points' and 'delta' specified");
        }
        
        if (!range.isMember("points") && range.isMember("delta")) {
            delta = range.get("delta", key == "temperature" ? DEFAULT_TEMP_DELTA : DEFAULT_FIELD_DELTA).asDouble();
            if (delta == 0.0 || ((finalVal - start) / delta + 1) <= 0) {
                throw vegas::ConfigurationException("Invalid delta value in " + key + " section");
            }
            points = static_cast<Index>((finalVal - start) / delta + 1);
        }
        
        Real val = start;
        for (Index i = 0; i < points; ++i) {
            values.push_back(val);
            val += delta;
        }
        
        // Handle cycle
        if (cycle) {
            val -= 2 * delta;
            for (Index i = 0; i < points - 1; ++i) {
                values.push_back(val);
                val -= delta;
            }
        }
    } else {
        // Default value
        singleValue = 0.0;
        values.push_back(0.0);
    }
    
    return values;
}

std::vector<Real> ConfigParser::expandSingleValue(Real value, size_t count) {
    return std::vector<Real>(count, value);
}

void ConfigParser::validate(const SimulationConfig& config) {
    // Check sample file
    if (config.sampleFile.empty()) {
        throw vegas::InvalidInputException("Sample file not specified");
    }
    
    std::ifstream infile(config.sampleFile);
    if (!infile.good()) {
        throw vegas::FileIOException("Cannot open sample file: " + config.sampleFile);
    }
    
    // Check MCS
    if (config.mcs < 10) {
        throw vegas::InvalidInputException("MCS must be at least 10");
    }
    
    // Check temperature/field vector sizes match
    if (config.temperatures.size() != config.fields.size()) {
        throw vegas::ConfigurationException(
            "Temperature count (" + std::to_string(config.temperatures.size()) + 
            ") does not match field count (" + std::to_string(config.fields.size()) + ")");
    }
    
    // Check initial state file if specified
    if (config.hasInitialState && !config.initialStateFile.empty()) {
        std::ifstream stateFile(config.initialStateFile);
        if (!stateFile.good()) {
            throw vegas::FileIOException("Cannot open initial state file: " + config.initialStateFile);
        }
    }
    
    // Check anisotropy files if specified
    for (const auto& file : config.anisotropyFiles) {
        std::ifstream anisoFile(file);
        if (!anisoFile.good()) {
            throw vegas::FileIOException("Cannot open anisotropy file: " + file);
        }
    }
}

} // namespace vegas
