#include "../include/system.h"
#include "../include/rlutil.h"
#include "../include/starter.h"
#include "../include/exception.h"
#include "../include/config_parser.h"
#include "../include/system_builder.h"
#include <cxxopts.hpp>
#include "json/json.h"

#include <iostream>
#include <string>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <vector>
#include <optional>

// Many functions for checking, errors and read the information were
// defined in a namespace STARTER into the file starter.h
using namespace STARTER;

// Version information
constexpr const char* VEGAS_VERSION = "2.3.0";
constexpr const char* VEGAS_BUILD_DATE = __DATE__;
constexpr const char* VEGAS_BUILD_TIME = __TIME__;

// Forward declarations
void print_version();
void print_header(bool use_color = true);
void run_simulation(const std::string& config_file, bool verbose = true, bool use_color = true,
    const std::optional<Index>& mcsOverride = std::nullopt,
    const std::optional<Index>& seedOverride = std::nullopt,
    const std::optional<Real>& kbOverride = std::nullopt,
    const std::optional<std::string>& outputOverride = std::nullopt,
    const std::optional<std::string>& sampleOverride = std::nullopt,
    const std::optional<std::string>& initialStateOverride = std::nullopt,
    const std::optional<std::string>& resumeFile = std::nullopt);
void analyze_output(const std::string& h5_file, bool verbose = true, bool use_color = true);
void validate_config(const std::string& config_file, bool verbose = true, bool use_color = true);
void print_info(bool verbose = true);
void report_overrides(const cxxopts::ParseResult& result, bool verbose);

// Helper function implementations
void report_overrides(const cxxopts::ParseResult& result, bool verbose) {
    if (!verbose) return;
    
    static const char* override_options[] = {
        "mcs", "seed", "kb", "output", "sample", "initialstate"
    };
    std::vector<std::string> present_overrides;
    
    for (const char* opt : override_options) {
        if (result.count(opt)) {
            present_overrides.push_back(opt);
        }
    }
    
    if (!present_overrides.empty()) {
        std::cout << "Applying configuration overrides: ";
        for (size_t i = 0; i < present_overrides.size(); ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << "--" << present_overrides[i];
        }
        std::cout << "\n";
    }
}

// Function implementations
void print_version() {
    std::cout << "VEGAS (VEctor General Atomistic Simulator) version " << VEGAS_VERSION << "\n"
              << "Build date: " << VEGAS_BUILD_DATE << " " << VEGAS_BUILD_TIME << "\n"
              << "License: MIT\n"
              << "Website: https://github.com/jdalzatec/vegas\n";
}

void print_header(bool use_color) {
    if (!use_color) {
        std::cout << "\n"
                  << " ________________________________________________________\n"
                  << "|                                                        |\n"
                  << "| ****************************************************** |\n"
                  << "| **                                                  ** |\n"
                  << "| **                   VEGAS                          ** |\n"
                  << "| **                                                  ** |\n"
                  << "| ****************************************************** |\n"
                  << "|________________________________________________________|\n\n";
        return;
    }
    
    std::cout << '\a';
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
    std::cout << " |" << std::endl;
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

void run_simulation(const std::string& config_file, bool verbose, bool use_color,
    const std::optional<Index>& mcsOverride,
    const std::optional<Index>& seedOverride,
    const std::optional<Real>& kbOverride,
    const std::optional<std::string>& outputOverride,
    const std::optional<std::string>& sampleOverride,
    const std::optional<std::string>& initialStateOverride,
    const std::optional<std::string>& resumeFile)
{
    if (verbose && use_color) {
        print_header(use_color);
    } else if (verbose) {
        std::cout << "VEGAS simulation started\n";
    }
    
    vegas::SimulationConfig config = vegas::ConfigParser::parse(config_file);
    config.applyOverrides(mcsOverride, seedOverride, kbOverride, outputOverride, sampleOverride, initialStateOverride);
    
    vegas::SystemBuilder builder;
    builder.withConfig(config);
    System system_ = builder.build();
    
    // Handle resume from checkpoint
    if (resumeFile.has_value()) {
        Index tempIndex = 0;
        bool loaded = system_.loadCheckpoint(resumeFile.value(), tempIndex);
        if (!loaded) {
            throw vegas::SimulationException("Failed to load checkpoint from: " + resumeFile.value());
        }
        if (verbose) {
            std::cout << "Resumed from checkpoint at temperature index: " << tempIndex << "\n";
        }
    }
    
    system_.cycle();
    
    if (verbose) {
        if (use_color) {
            rlutil::setColor(rlutil::LIGHTGREEN);
            std::cout << "Successful completion !!!" << std::endl;
            rlutil::resetColor();
            std::cout << '\a';
        } else {
            std::cout << "Successful completion !!!" << std::endl;
        }
    }
}

void analyze_output(const std::string& h5_file, bool verbose, bool use_color) {
    std::ifstream h5_check(h5_file);
    if (!h5_check.good()) {
        throw vegas::FileIOException("HDF5 file not found: " + h5_file);
    }
    h5_check.close();
    
    std::string analyzer = "analyzers/vegas-analyzer-lite.py";
    std::ifstream analyzer_check(analyzer);
    if (!analyzer_check.good()) {
        if (use_color) rlutil::setColor(rlutil::LIGHTRED);
        std::cerr << "Analyzer script not found: " << analyzer << "\n";
        std::cerr << "Please ensure the analyzers/ directory exists and contains the Python scripts.\n";
        if (use_color) rlutil::resetColor();
        throw vegas::FileIOException("Analyzer script not found: " + analyzer);
    }
    analyzer_check.close();
    
    std::string python = "python3";
    std::ifstream venv_check(".venv/bin/python3");
    if (venv_check.good()) {
        python = ".venv/bin/python3";
    }
    venv_check.close();
    
    if (verbose) {
        std::cout << "Running analysis on: " << h5_file << "\n";
        std::cout << "Using analyzer: " << analyzer << "\n\n";
    }
    
    std::string cmd = python + " " + analyzer + " " + h5_file;
    int result = std::system(cmd.c_str());
    
    if (result != 0) {
        if (use_color) rlutil::setColor(rlutil::LIGHTRED);
        std::cerr << "Analyzer failed with exit code: " << result << "\n";
        if (use_color) rlutil::resetColor();
        throw vegas::SimulationException("Analyzer failed with exit code: " + std::to_string(result));
    }
    
    if (verbose) {
        std::cout << "\nAnalysis complete.\n";
    }
}

void validate_config(const std::string& config_file, bool verbose, bool use_color) {
    if (verbose) {
        std::cout << "Validating configuration: " << config_file << "\n\n";
    }
    
    try {
        vegas::SimulationConfig config = vegas::ConfigParser::parse(config_file);
        
        if (verbose) {
            std::cout << "Configuration valid!\n\n";
            std::cout << "  Sample file: " << config.sampleFile << "\n";
            std::cout << "  Output file: " << config.outputFile << "\n";
            std::cout << "  MCS: " << config.mcs << "\n";
            std::cout << "  Seed: " << config.seed << "\n";
            std::cout << "  Boltzmann constant: " << config.kb << "\n";
            std::cout << "  Temperature points: " << config.temperatures.size() << "\n";
            if (!config.temperatures.empty()) {
                std::cout << "    Range: " << config.temperatures.front() 
                          << " to " << config.temperatures.back() << "\n";
            }
            std::cout << "  Field points: " << config.fields.size() << "\n";
            if (!config.fields.empty()) {
                std::cout << "    Range: " << config.fields.front() 
                          << " to " << config.fields.back() << "\n";
            }
            if (config.hasInitialState) {
                std::cout << "  Initial state: " << config.initialStateFile << "\n";
            }
            if (config.hasAnisotropy) {
                std::cout << "  Anisotropy files: " << config.anisotropyFiles.size() << "\n";
                for (const auto& file : config.anisotropyFiles) {
                    std::cout << "    - " << file << "\n";
                }
            }
            std::cout << "\nValidation complete.\n";
        } else {
            std::cout << "Configuration valid: " << config_file << "\n";
        }
    } catch (const vegas::VEGASException& e) {
        if (use_color) rlutil::setColor(rlutil::LIGHTRED);
        std::cerr << "Validation failed: " << e.what() << "\n";
        if (use_color) rlutil::resetColor();
        throw;
    }
}

void print_info(bool verbose) {
    print_version();
    if (verbose) {
        std::cout << "\nAvailable spin models: Heisenberg, Ising, QuantumIsing, Adaptive, Cone, HN\n";
        std::cout << "Output format: HDF5 with magnetization, energy, temperature, field datasets\n";
        std::cout << "Python analyzers available in analyzers/ directory\n";
    }
}

int main(int argc, char const *argv[]) {
    // Variables accessible in catch blocks
    bool use_color = true;
    bool verbose = true;
    
    try {
        cxxopts::Options options("vegas", "VEGAS - VEctor General Atomistic Simulator");
        
        // Global options
        options.add_options()
            ("h,help", "Show help")
            ("v,version", "Show version information")
            ("V,verbose", "Verbose output (default)")
            ("q,quiet", "Quiet mode (minimal output)")
            ("no-color", "Disable colored output")
            ("c,config", "Configuration file (for backward compatibility)", cxxopts::value<std::string>())
            ("r,resume", "Resume from checkpoint", cxxopts::value<std::string>())
            // Simulation overrides (only used with 'run' command)
            ("mcs", "Override MCS (Monte Carlo steps)", cxxopts::value<Index>())
            ("seed", "Override random seed", cxxopts::value<Index>())
            ("kb", "Override Boltzmann constant", cxxopts::value<Real>())
            ("output", "Override output file path", cxxopts::value<std::string>())
            ("sample", "Override sample file path", cxxopts::value<std::string>())
            ("initialstate", "Override initial state file", cxxopts::value<std::string>());

        // Subcommands (positional arguments)
        options.add_options()
            ("command", "Subcommand: run, analyze, validate, info", cxxopts::value<std::string>())
            ("input", "Input file (config.json or output.h5)", cxxopts::value<std::string>());
        
        options.parse_positional({"command", "input"});
        options.positional_help("[command] [input]");
        options.show_positional_help();
        
        // Parse command line
        auto result = options.parse(argc, argv);
        
        // Handle global options
        use_color = !result.count("no-color");
        verbose = true;
        if (result.count("quiet")) verbose = false;
        if (result.count("verbose")) verbose = true;
        
        // Show version
        if (result.count("version")) {
            print_version();
            return 0;
        }
        
        // Show help
        if (result.count("help") || (argc == 1 && !result.count("config"))) {
            std::cout << options.help() << "\n\n";
            std::cout << "Examples:\n"
                      << "  ./vegas run config.json          Run simulation with config.json\n"
                      << "  ./vegas analyze output.h5        Analyze simulation results\n"
                      << "  ./vegas validate config.json     Validate configuration file\n"
                      << "  ./vegas info                     Show version and information\n"
                      << "  ./vegas config.json              Backward compatibility (runs simulation)\n"
                      << "  ./vegas --version                Show version\n"
                      << "  ./vegas --quiet config.json      Run quietly\n"
                      << "  ./vegas --no-color run config.json Disable colors\n";
            return 0;
        }
        
        // Backward compatibility: if a config file is provided without subcommand
        if (result.count("config")) {
            std::string config_file = result["config"].as<std::string>();
            if (verbose) {
                std::cout << "Note: Using backward compatibility mode. Consider using 'vegas run " << config_file << "'\n";
            }
            
            report_overrides(result, verbose);
            
            run_simulation(config_file, verbose, use_color,
                result.count("mcs") ? std::optional<Index>(result["mcs"].as<Index>()) : std::nullopt,
                result.count("seed") ? std::optional<Index>(result["seed"].as<Index>()) : std::nullopt,
                result.count("kb") ? std::optional<Real>(result["kb"].as<Real>()) : std::nullopt,
                result.count("output") ? std::optional<std::string>(result["output"].as<std::string>()) : std::nullopt,
                result.count("sample") ? std::optional<std::string>(result["sample"].as<std::string>()) : std::nullopt,
                result.count("initialstate") ? std::optional<std::string>(result["initialstate"].as<std::string>()) : std::nullopt
            );
            return 0;
        }
        
        // Handle subcommands
        if (!result.count("command")) {
            std::cerr << "Error: No command specified. Use --help for usage.\n";
            return EXIT_FAILURE;
        }
        
        std::string command = result["command"].as<std::string>();
        
        if (command == "run") {
            if (!result.count("input")) {
                std::cerr << "Error: run command requires a configuration file.\n";
                return EXIT_FAILURE;
            }
            std::string config_file = result["input"].as<std::string>();
            
            report_overrides(result, verbose);
            
            run_simulation(config_file, verbose, use_color,
                result.count("mcs") ? std::optional<Index>(result["mcs"].as<Index>()) : std::nullopt,
                result.count("seed") ? std::optional<Index>(result["seed"].as<Index>()) : std::nullopt,
                result.count("kb") ? std::optional<Real>(result["kb"].as<Real>()) : std::nullopt,
                result.count("output") ? std::optional<std::string>(result["output"].as<std::string>()) : std::nullopt,
                result.count("sample") ? std::optional<std::string>(result["sample"].as<std::string>()) : std::nullopt,
                result.count("initialstate") ? std::optional<std::string>(result["initialstate"].as<std::string>()) : std::nullopt
            );
        } else if (command == "analyze") {
            if (!result.count("input")) {
                std::cerr << "Error: analyze command requires an HDF5 output file.\n";
                return EXIT_FAILURE;
            }
            std::string h5_file = result["input"].as<std::string>();
            analyze_output(h5_file, verbose, use_color);
        } else if (command == "validate") {
            if (!result.count("input")) {
                std::cerr << "Error: validate command requires a configuration file.\n";
                return EXIT_FAILURE;
            }
            std::string config_file = result["input"].as<std::string>();
            validate_config(config_file, verbose, use_color);
        } else if (command == "info") {
            print_info(verbose);
        } else {
            // Check if the command is actually a .json file (backward compatibility)
            if (command.size() > 5 && command.substr(command.size() - 5) == ".json") {
                if (verbose) {
                    std::cout << "Note: Using backward compatibility mode. Consider using 'vegas run " << command << "'\n";
                }
                
                report_overrides(result, verbose);
                
                run_simulation(command, verbose, use_color,
                    result.count("mcs") ? std::optional<Index>(result["mcs"].as<Index>()) : std::nullopt,
                    result.count("seed") ? std::optional<Index>(result["seed"].as<Index>()) : std::nullopt,
                    result.count("kb") ? std::optional<Real>(result["kb"].as<Real>()) : std::nullopt,
                    result.count("output") ? std::optional<std::string>(result["output"].as<std::string>()) : std::nullopt,
                    result.count("sample") ? std::optional<std::string>(result["sample"].as<std::string>()) : std::nullopt,
                    result.count("initialstate") ? std::optional<std::string>(result["initialstate"].as<std::string>()) : std::nullopt
                );
            } else {
                std::cerr << "Error: Unknown command '" << command << "'. Use --help for usage.\n";
                return EXIT_FAILURE;
            }
        }
        
        return 0;
    } catch (const cxxopts::exceptions::exception& e) {
        std::cerr << "Error parsing command line: " << e.what() << "\n";
        return EXIT_FAILURE;
    } catch (const vegas::VEGASException& e) {
        if (use_color) rlutil::setColor(rlutil::LIGHTRED);
        std::cerr << "Simulation error: " << e.what() << "\n";
        if (use_color) rlutil::resetColor();
        return EXIT_FAILURE;
    } catch (const std::exception& e) {
        if (use_color) rlutil::setColor(rlutil::LIGHTRED);
        std::cerr << "Unexpected error: " << e.what() << "\n";
        if (use_color) rlutil::resetColor();
        return EXIT_FAILURE;
    }
}
