// SystemBuilder unit tests for VEGAS simulation package

#include "../include/system_builder.h"
#include "../include/config_parser.h"
#include "../include/exception.h"
#include "../include/params.h"
#include <iostream>
#include <cassert>

void test_builder_fluent_interface()
{
    std::cout << "Testing SystemBuilder fluent interface... ";
    
    vegas::SystemBuilder builder;
    
    vegas::SystemBuilder& ref1 = builder.setSampleFile("test.txt");
    vegas::SystemBuilder& ref2 = builder.setMcs(100);
    vegas::SystemBuilder& ref3 = builder.setSeed(42);
    vegas::SystemBuilder& ref4 = builder.setKb(1.0);
    vegas::SystemBuilder& ref5 = builder.setTemperatures({1.0});
    vegas::SystemBuilder& ref6 = builder.setFields({0.0});
    
    assert(&ref1 == &builder);
    assert(&ref2 == &builder);
    assert(&ref3 == &builder);
    assert(&ref4 == &builder);
    assert(&ref5 == &builder);
    assert(&ref6 == &builder);
    
    std::cout << "PASS\n";
}

void test_builder_set_sample_file()
{
    std::cout << "Testing SystemBuilder set sample file... ";
    
    vegas::SystemBuilder builder;
    builder.setSampleFile("test_system/ferromagnetic_chain.txt");
    
    std::cout << "PASS\n";
}

void test_builder_set_temperatures()
{
    std::cout << "Testing SystemBuilder set temperatures... ";
    
    vegas::SystemBuilder builder;
    std::vector<Real> temps = {0.5, 1.0, 1.5, 2.0};
    builder.setTemperatures(temps);
    
    std::cout << "PASS\n";
}

void test_builder_set_fields()
{
    std::cout << "Testing SystemBuilder set fields... ";
    
    vegas::SystemBuilder builder;
    std::vector<Real> fields = {0.0, 0.5, 1.0};
    builder.setFields(fields);
    
    std::cout << "PASS\n";
}

void test_builder_with_config()
{
    std::cout << "Testing SystemBuilder with config... ";
    
    vegas::SimulationConfig config;
    config.sampleFile = "test_system/ferromagnetic_chain.txt";
    config.mcs = 100;
    config.seed = 42;
    config.kb = 1.0;
    config.temperatures = {1.0};
    config.fields = {0.0};
    config.outputFile = "test_output.h5";
    config.hasInitialState = false;
    config.hasAnisotropy = false;
    
    vegas::SystemBuilder builder;
    builder.withConfig(config);
    
    std::cout << "PASS\n";
}

void test_builder_build_minimal()
{
    std::cout << "Testing SystemBuilder build minimal... ";
    
    vegas::SystemBuilder builder;
    builder.setSampleFile("test_system/ferromagnetic_chain.txt")
           .setTemperatures({1.0})
           .setFields({0.0})
           .setMcs(10)
           .setSeed(42)
           .setOutputFile("test_builder_output.h5");
    
    System system = builder.build();
    
    assert(system.getLattice().getAtoms().size() == 2);
    assert(system.getSeed() == 42);
    
    std::remove("test_builder_output.h5");
    
    std::cout << "PASS\n";
}

void test_builder_missing_sample()
{
    std::cout << "Testing SystemBuilder missing sample file... ";
    
    vegas::SystemBuilder builder;
    builder.setTemperatures({1.0})
           .setFields({0.0});
    
    bool caught = false;
    try {
        builder.build();
    } catch (const vegas::InvalidInputException& e) {
        caught = true;
    } catch (...) {
        caught = true;
    }
    assert(caught);
    
    std::cout << "PASS\n";
}

void test_builder_missing_temps_fields()
{
    std::cout << "Testing SystemBuilder missing temps/fields... ";
    
    vegas::SystemBuilder builder;
    builder.setSampleFile("test_system/ferromagnetic_chain.txt");
    
    bool caught = false;
    try {
        builder.build();
    } catch (const vegas::InvalidInputException& e) {
        caught = true;
    } catch (...) {
        caught = true;
    }
    assert(caught);
    
    std::cout << "PASS\n";
}

int main()
{
    std::cout << "Running SystemBuilder tests for vegas...\n\n";
    
    try {
        test_builder_fluent_interface();
        test_builder_set_sample_file();
        test_builder_set_temperatures();
        test_builder_set_fields();
        test_builder_with_config();
        test_builder_build_minimal();
        test_builder_missing_sample();
        test_builder_missing_temps_fields();
        
        std::cout << "\nAll SystemBuilder tests passed!\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\nTest failed with exception: " << e.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "\nTest failed with unknown exception\n";
        return 1;
    }
}
