#include "../include/system_builder.h"
#include "../include/exception.h"

#include <fstream>

namespace vegas {

SystemBuilder& SystemBuilder::setSampleFile(const std::string& file) {
    sampleFile_ = file;
    return *this;
}

SystemBuilder& SystemBuilder::setOutputFile(const std::string& file) {
    outputFile_ = file;
    return *this;
}

SystemBuilder& SystemBuilder::setMcs(Index mcs) {
    mcs_ = mcs;
    return *this;
}

SystemBuilder& SystemBuilder::setSeed(Index seed) {
    seed_ = seed;
    return *this;
}

SystemBuilder& SystemBuilder::setKb(Real kb) {
    kb_ = kb;
    return *this;
}

SystemBuilder& SystemBuilder::setTemperatures(const std::vector<Real>& temps) {
    temperatures_ = temps;
    return *this;
}

SystemBuilder& SystemBuilder::setFields(const std::vector<Real>& fields) {
    fields_ = fields;
    return *this;
}

SystemBuilder& SystemBuilder::setInitialState(const std::string& file) {
    initialStateFile_ = file;
    hasInitialState_ = !file.empty();
    return *this;
}

SystemBuilder& SystemBuilder::setAnisotropyFiles(const std::vector<std::string>& files) {
    anisotropyFiles_ = files;
    hasAnisotropy_ = !files.empty();
    return *this;
}

SystemBuilder& SystemBuilder::withConfig(const SimulationConfig& config) {
    sampleFile_ = config.sampleFile;
    outputFile_ = config.outputFile;
    mcs_ = config.mcs;
    seed_ = config.seed;
    kb_ = config.kb;
    temperatures_ = config.temperatures;
    fields_ = config.fields;
    initialStateFile_ = config.initialStateFile;
    hasInitialState_ = config.hasInitialState;
    anisotropyFiles_ = config.anisotropyFiles;
    hasAnisotropy_ = config.hasAnisotropy;
    return *this;
}

System SystemBuilder::build() {
    // Validate required fields
    if (sampleFile_.empty()) {
        throw vegas::InvalidInputException("Sample file not specified");
    }
    
    if (temperatures_.empty() || fields_.empty()) {
        throw vegas::InvalidInputException("Temperatures and fields must be specified");
    }
    
    // Use sample filename as default output
    std::string outputFile = outputFile_.empty() ? sampleFile_ + ".h5" : outputFile_;
    
    // Create the system
    System system(sampleFile_, temperatures_, fields_, mcs_, seed_, outputFile, kb_);
    
    // Set initial state
    if (hasInitialState_ && !initialStateFile_.empty()) {
        system.setState(initialStateFile_);
    } else {
        system.randomizeSpins();
    }
    
    // Set anisotropy
    if (hasAnisotropy_ && !anisotropyFiles_.empty()) {
        system.setAnisotropies(anisotropyFiles_);
    }
    
    return system;
}

} // namespace vegas
