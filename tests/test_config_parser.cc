// ConfigParser unit tests for VEGAS simulation package

#include "../include/config_parser.h"
#include "../include/exception.h"
#include "../include/params.h"
#include <iostream>
#include <cassert>
#include <fstream>

void test_parse_basic_config()
{
    std::cout << "Testing ConfigParser basic config... ";
    
    vegas::SimulationConfig config = vegas::ConfigParser::parse("test_system/test_config.json");
    
    assert(config.sampleFile == "ferromagnetic_chain.txt");
    assert(config.mcs == 100);
    assert(config.kb == 1.0);
    assert(config.seed == 12345);
    assert(config.temperatures.size() == 3);
    assert(config.fields.size() == 3);
    assert(config.hasInitialState == true);
    assert(config.outputFile == "test_results.h5");
    
    std::cout << "PASS\n";
}

void test_parse_temperature_range()
{
    std::cout << "Testing ConfigParser temperature range... ";
    
    vegas::SimulationConfig config = vegas::ConfigParser::parse("test_system/range_config.json");
    
    assert(config.temperatures.size() == 5);
    assert(fp_equal(config.temperatures[0], 0.5));
    assert(fp_equal(config.temperatures[4], 2.0));
    
    std::cout << "PASS\n";
}

void test_parse_field_range()
{
    std::cout << "Testing ConfigParser field range... ";
    
    vegas::SimulationConfig config = vegas::ConfigParser::parse("test_system/range_config.json");
    
    assert(config.fields.size() == 5);
    assert(fp_equal(config.fields[0], 0.0));
    assert(fp_equal(config.fields[4], 1.0));
    
    std::cout << "PASS\n";
}

void test_parse_single_values()
{
    std::cout << "Testing ConfigParser single values... ";
    
    vegas::SimulationConfig config = vegas::ConfigParser::parse("test_system/valid_config.json");
    
    assert(config.temperatures.size() == 1);
    assert(fp_equal(config.temperatures[0], 1.0));
    assert(config.fields.size() == 1);
    assert(fp_equal(config.fields[0], 0.0));
    
    std::cout << "PASS\n";
}

void test_parse_arrays()
{
    std::cout << "Testing ConfigParser array values... ";
    
    vegas::SimulationConfig config = vegas::ConfigParser::parse("test_system/array_config.json");
    
    assert(config.temperatures.size() == 4);
    assert(fp_equal(config.temperatures[0], 0.5));
    assert(fp_equal(config.temperatures[3], 2.0));
    
    assert(config.fields.size() == 4);
    assert(fp_equal(config.fields[0], 0.0));
    assert(fp_equal(config.fields[3], 1.5));
    
    std::cout << "PASS\n";
}

void test_validate_missing_file()
{
    std::cout << "Testing ConfigParser missing file validation... ";
    
    bool caught = false;
    try {
        vegas::ConfigParser::parse("test_system/nonexistent_config.json");
    } catch (const vegas::FileIOException& e) {
        caught = true;
    } catch (...) {
        caught = true;
    }
    assert(caught);
    
    std::cout << "PASS\n";
}

void test_validate_invalid_mcs()
{
    std::cout << "Testing ConfigParser invalid MCS validation... ";
    
    bool caught = false;
    try {
        vegas::ConfigParser::parse("test_system/invalid_mcs_config.json");
    } catch (const vegas::InvalidInputException& e) {
        caught = true;
    } catch (...) {
        caught = true;
    }
    assert(caught);
    
    std::cout << "PASS\n";
}

void test_default_values()
{
    std::cout << "Testing ConfigParser default values... ";
    
    vegas::SimulationConfig config = vegas::ConfigParser::parse("test_system/valid_config.json");
    
    assert(config.mcs == 100);
    assert(config.seed == 42);
    assert(fp_equal(config.kb, 1.0));
    assert(config.hasAnisotropy == false);
    
    std::cout << "PASS\n";
}

int main()
{
    std::cout << "Running ConfigParser tests for vegas...\n\n";
    
    try {
        test_parse_basic_config();
        test_parse_temperature_range();
        test_parse_field_range();
        test_parse_single_values();
        test_parse_arrays();
        test_validate_missing_file();
        test_validate_invalid_mcs();
        test_default_values();
        
        std::cout << "\nAll ConfigParser tests passed!\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\nTest failed with exception: " << e.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "\nTest failed with unknown exception\n";
        return 1;
    }
}
