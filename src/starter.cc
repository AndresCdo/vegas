#include "../include/starter.h"
#include "../include/system.h"
#include "../include/config_parser.h"
#include "../include/system_builder.h"
#include <vector>
#include <string>
#include <stdexcept>

namespace STARTER {

    void CHECK(std::string sample,
               Index mcs,
               const std::vector<Real>& temps,
               const std::vector<Real>& fields)
    {
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
        
        // Basic check for path traversal attempts
        if (filename.find("..") != std::string::npos) {
            throw vegas::InvalidInputException("Filename contains path traversal attempt: " + filename);
        }
        
        // Check for absolute paths (optional security measure)
        if (!filename.empty() && filename[0] == '/') {
            throw vegas::InvalidInputException("Absolute paths are not allowed: " + filename);
        }
        
        // Try to open the file
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
