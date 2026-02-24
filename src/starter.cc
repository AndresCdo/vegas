#include "../include/starter.h"
#include "../include/system.h"
#include "../include/config_parser.h"
#include "../include/system_builder.h"
#include <vector>
#include <string>
#include <stdexcept>
#include <unistd.h>
#include <climits>
#include <cstdlib>

namespace STARTER {

    void CHECK(std::string sample,
               Index mcs,
               const std::vector<Real>& temps,
               const std::vector<Real>& fields)
    {
        CHECKFILE(sample);
        
        std::ifstream infile(sample);
        if (!infile.good())
            throw vegas::FileIOException("The sample file can't open or doesn't exist !!!");

        if (mcs < 10)
            throw vegas::InvalidInputException("The number of MCS must be greater than 10 !!!");

        if (temps.size() != fields.size())
            throw vegas::ConfigurationException("The amount of temperatures and fields are differents ( " +
                std::to_string(
                temps.size()) + " != " + std::to_string(fields.size()) + " )!!!");
    }

    void CHECKFILE(const std::string& filename)
    {
        // Check for empty filename
        if (filename.empty()) {
            throw vegas::InvalidInputException("Filename is empty!");
        }
        
        // Resolve absolute path
        char resolved[PATH_MAX];
        if (realpath(filename.c_str(), resolved) == nullptr) {
            // File may not exist yet (output files), but we still need to check path safety
            // For now, allow if no '..' components escape current directory
            // Simple check: ensure filename doesn't contain /../
            std::string normalized = filename;
            size_t pos = 0;
            while ((pos = normalized.find("/../", pos)) != std::string::npos) {
                throw vegas::InvalidInputException("Path traversal attempt detected: " + filename);
            }
            // Also check for leading ../
            if (normalized.find("../") == 0) {
                throw vegas::InvalidInputException("Path traversal attempt detected: " + filename);
            }
            // Check for trailing /..
            if (normalized.size() >= 3 && normalized.substr(normalized.size() - 3) == "/..") {
                throw vegas::InvalidInputException("Path traversal attempt detected: " + filename);
            }
        } else {
            // File exists, check it's within the project directory (parent of build directory)
            char cwd[PATH_MAX];
            if (getcwd(cwd, sizeof(cwd)) == nullptr) {
                throw vegas::FileIOException("Cannot get current working directory");
            }
            
            // Normalize both paths
            std::string resolved_str(resolved);
            std::string cwd_str(cwd);
            
            // Remove trailing slash if present
            if (!cwd_str.empty() && cwd_str.back() == '/') {
                cwd_str.pop_back();
            }
            
            // Find parent directory of cwd (project root)
            size_t last_slash = cwd_str.find_last_of('/');
            if (last_slash == std::string::npos) {
                // cwd is root directory '/', unlikely
                throw vegas::FileIOException("Unexpected current directory structure");
            }
            std::string parent_dir = cwd_str.substr(0, last_slash);
            
            // Ensure parent_dir ends with '/' for prefix check
            if (!parent_dir.empty() && parent_dir.back() != '/') {
                parent_dir += '/';
            }
            
            // Ensure resolved path starts with parent_dir
            if (resolved_str.find(parent_dir) != 0) {
                throw vegas::InvalidInputException("File path escapes project directory: " + filename);
            }
        }
        
        // Try to open the file (for reading)
        std::ifstream infile(filename);
        if (!infile.good()) {
            throw vegas::FileIOException("Cannot open file: " + filename + " (file doesn't exist or permission denied)");
        }
        
        // Check if file is readable by trying to read first character
        char test;
        infile.get(test);
        if (infile.fail() && !infile.eof()) {
            throw vegas::FileIOException("File is not readable: " + filename);
        }
    }

    void HEADER()
    {
        rlutil::saveDefaultColor();
        std::cout << std::endl;
        rlutil::setColor(rlutil::LIGHTCYAN);
        std::cout << " ________________________________________________________"
                  << std::endl;
        std::cout << "|                                                        |"
                  << std::endl;
        std::cout << "| ";
        rlutil::setColor(rlutil::LIGHTGREEN);
        std::cout << "******************************************************";
        rlutil::setColor(rlutil::LIGHTCYAN);
        std::cout << " |" << std::endl;
        std::cout << "| ";
        rlutil::setColor(rlutil::LIGHTGREEN);
        std::cout << "**                                                  **";
        rlutil::setColor(rlutil::LIGHTCYAN);
        std::cout << " |" << std::endl;
        std::cout << "| ";
        rlutil::setColor(rlutil::LIGHTGREEN);
        std::cout << "**                   ";
        rlutil::setColor(rlutil::LIGHTRED);
        std::cout << "VEGAS";
        rlutil::setColor(rlutil::LIGHTGREEN);
        std::cout << "                          **";
        rlutil::setColor(rlutil::LIGHTCYAN);
        std::cout << " |" << std::endl;
        std::cout << "| ";
        rlutil::setColor(rlutil::LIGHTGREEN);
        std::cout << "**                                                  **";
        rlutil::setColor(rlutil::LIGHTCYAN);
        std::cout <<" |" << std::endl;
        std::cout << "| ";
        rlutil::setColor(rlutil::LIGHTGREEN);
        std::cout << "******************************************************";
        rlutil::setColor(rlutil::LIGHTCYAN);
        std::cout << " |" << std::endl;
        std::cout << "|_______________________________"
                  << "_________________________|"
                  << std::endl;
        std::cout << std::endl;
        rlutil::resetColor();
    }

    void PRINT_VALUES(System& system_,
                      const std::string& sample,
                      const Index& mcs,
                      const std::string& out,
                      const Real& kb,
                      const Index& seed,
                      const std::string& initialstate,
                      const std::vector<std::string>& anisotropyfiles)
    {
        // The initial values are printed for the user visualization.
        std::cout << "\t\tSample file = \n\t\t\t" << sample << std::endl;
        std::cout << "\t\tNum MCS = \n\t\t\t" << mcs << std::endl;
        std::cout << "\t\tOut file = \n\t\t\t" << out << std::endl;
        if (initialstate != "")
        {
            std::cout << "\t\tInitial state file = \n\t\t\t" << initialstate << std::endl;
        }
        if (anisotropyfiles.size() != 0)
        {
            for (auto&& anisotropyfile : anisotropyfiles)
                std::cout << "\t\tAnisotropy file = \n\t\t\t" << anisotropyfile << std::endl;
        }

        std::cout << std::endl;
        std::cout << "\t\tNum Ions = \n\t\t\t" << system_.getLattice().getAtoms().size() << std::endl;

        // The amount of ions are printed for type ion.
        for (auto&& type : system_.getLattice().getMapTypeIndexes())
        {
            std::cout << "\t\tNum " << type.first << " Ions  = \n\t\t\t" << system_.getLattice().getSizesByIndex().at(type.second) << std::endl;
        }

        std::cout << "\t\tkb = \n\t\t\t" << kb << std::endl;
        std::cout << "\t\tseed = \n\t\t\t" << system_.getSeed() << std::endl;

        std::cout << std::endl;
        std::cout << std::endl;

    }

    System CREATE_SYSTEM(std::string jsonfile, bool print)
    {
        // Use the new ConfigParser to parse the configuration
        vegas::SimulationConfig config = vegas::ConfigParser::parse(jsonfile);
        
        // Use the SystemBuilder to construct the system
        vegas::SystemBuilder builder;
        builder.withConfig(config);
        
        System system_ = builder.build();
        
        // Print values if requested
        if (print) {
            PRINT_VALUES(system_, config.sampleFile, config.mcs, config.outputFile, 
                         config.kb, config.seed, config.initialStateFile, config.anisotropyFiles);
        }

        return system_;
    }
    
    System CREATE_SYSTEM(std::string jsonfile, bool print, const Json::Value& overrides)
    {
        (void)jsonfile; (void)print; (void)overrides; // suppress unused parameter warnings
        throw vegas::ConfigurationException("Configuration overrides are not yet implemented");
    }
};
