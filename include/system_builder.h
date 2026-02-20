#ifndef VEGAS_SYSTEM_BUILDER_H
#define VEGAS_SYSTEM_BUILDER_H

#include "system.h"
#include "config_parser.h"

namespace vegas {

class SystemBuilder {
public:
    SystemBuilder() = default;
    
    SystemBuilder& setSampleFile(const std::string& file);
    SystemBuilder& setOutputFile(const std::string& file);
    SystemBuilder& setMcs(Index mcs);
    SystemBuilder& setSeed(Index seed);
    SystemBuilder& setKb(Real kb);
    SystemBuilder& setTemperatures(const std::vector<Real>& temps);
    SystemBuilder& setFields(const std::vector<Real>& fields);
    SystemBuilder& setInitialState(const std::string& file);
    SystemBuilder& setAnisotropyFiles(const std::vector<std::string>& files);
    SystemBuilder& withConfig(const SimulationConfig& config);
    
    System build();
    
private:
    std::string sampleFile_;
    std::string outputFile_;
    std::string initialStateFile_;
    std::vector<std::string> anisotropyFiles_;
    std::vector<Real> temperatures_;
    std::vector<Real> fields_;
    Index mcs_ = DEFAULT_MCS;
    Index seed_ = 0;
    Real kb_ = 1.0;
    bool hasInitialState_ = false;
    bool hasAnisotropy_ = false;
};

} // namespace vegas

#endif // VEGAS_SYSTEM_BUILDER_H
