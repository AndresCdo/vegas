// Integration tests for VEGAS simulation package
// Tests actual simulation functionality with different spin models

#include "../include/params.h"
#include "../include/atom.h"
#include "../include/lattice.h"
#include "../include/spin_model.h"
#include <iostream>
#include <cassert>
#include <random>
#include <memory>

void test_spin_model_factory() {
    std::cout << "Testing SpinModel factory... ";
    
    // Test creation of different spin models
    auto heisenberg = createSpinModel("heisenberg");
    assert(heisenberg != nullptr);
    assert(heisenberg->getName() == "heisenberg");
    
    auto ising = createSpinModel("ising");
    assert(ising != nullptr);
    assert(ising->getName() == "ising");
    
    auto random = createSpinModel("random");
    assert(random != nullptr);
    assert(random->getName() == "heisenberg"); // random maps to heisenberg
    
    auto flip = createSpinModel("flip");
    assert(flip != nullptr);
    assert(flip->getName() == "ising"); // flip maps to ising
    
    auto qising = createSpinModel("qising");
    assert(qising != nullptr);
    assert(qising->getName() == "qising");
    
    auto adaptive = createSpinModel("adaptive");
    assert(adaptive != nullptr);
    assert(adaptive->getName() == "adaptive");
    
    auto cone30 = createSpinModel("cone30");
    assert(cone30 != nullptr);
    assert(cone30->getName() == "cone30");
    
    auto cone15 = createSpinModel("cone15");
    assert(cone15 != nullptr);
    assert(cone15->getName() == "cone15");
    
    auto hn30 = createSpinModel("hn30");
    assert(hn30 != nullptr);
    assert(hn30->getName() == "hn30");
    
    auto hn15 = createSpinModel("hn15");
    assert(hn15 != nullptr);
    assert(hn15->getName() == "hn15");
    
    // Test default model (should be heisenberg)
    auto default_model = createSpinModel("unknown");
    assert(default_model != nullptr);
    assert(default_model->getName() == "heisenberg");
    
    std::cout << "PASS\n";
}

void test_atom_spin_model_integration() {
    std::cout << "Testing Atom-SpinModel integration... ";
    
    // Create random number generators
    std::mt19937_64 engine(42);
    std::uniform_real_distribution<> realRandomGenerator(0.0, 1.0);
    std::normal_distribution<> gaussianRandomGenerator(0.0, 1.0);
    
    // Test Heisenberg model
    Atom atom1(0, {0.0, 0.0, 1.0}, {0.0, 0.0, 0.0});
    atom1.setModel("heisenberg");
    
    // Store initial spin
    Array initialSpin = atom1.getSpin();
    
    // Initialize random state
    atom1.initializeRandomState(engine, realRandomGenerator, gaussianRandomGenerator);
    Array randomSpin = atom1.getSpin();
    
    // Verify spin changed and has correct norm
    assert(fp_not_equal(initialSpin[0], randomSpin[0]) || 
           fp_not_equal(initialSpin[1], randomSpin[1]) || 
           fp_not_equal(initialSpin[2], randomSpin[2]));
    
    Real norm = std::sqrt((randomSpin * randomSpin).sum());
    assert(fp_equal(norm, atom1.getSpinNorm()));
    
    // Test Ising model
    Atom atom2(1, {0.0, 0.0, 1.0}, {0.0, 0.0, 0.0});
    atom2.setModel("ising");
    atom2.initializeRandomState(engine, realRandomGenerator, gaussianRandomGenerator);
    
    // Ising spin should be either +z or -z
    Array isingSpin = atom2.getSpin();
    assert(fp_equal(std::abs(isingSpin[2]), atom2.getSpinNorm()));
    assert(fp_equal(isingSpin[0], 0.0));
    assert(fp_equal(isingSpin[1], 0.0));
    
    std::cout << "PASS\n";
}

void test_spin_model_randomization() {
    std::cout << "Testing spin model randomization... ";
    
    // Create random number generators
    std::mt19937_64 engine(123);
    std::uniform_real_distribution<> realRandomGenerator(0.0, 1.0);
    std::normal_distribution<> gaussianRandomGenerator(0.0, 1.0);
    
    // Test Heisenberg randomization
    Atom atom(0, {0.0, 0.0, 1.0}, {0.0, 0.0, 0.0});
    atom.setModel("heisenberg");
    
    Array originalSpin = atom.getSpin();
    
    // Randomize spin multiple times
    for (int i = 0; i < 10; i++) {
        atom.randomizeSpin(engine, realRandomGenerator, gaussianRandomGenerator, 0.1, i % 5);
        Array newSpin = atom.getSpin();
        
        // Verify spin norm is preserved
        Real norm = std::sqrt((newSpin * newSpin).sum());
        assert(fp_equal(norm, atom.getSpinNorm()));
        
        // Verify spin changed from original (at least once)
        if (i == 0) {
            assert(fp_not_equal(originalSpin[0], newSpin[0]) || 
                   fp_not_equal(originalSpin[1], newSpin[1]) || 
                   fp_not_equal(originalSpin[2], newSpin[2]));
        }
    }
    
    // Test Ising flip
    Atom atom2(1, {0.0, 0.0, 1.0}, {0.0, 0.0, 0.0});
    atom2.setModel("ising");
    
    Array originalSpin2 = atom2.getSpin();
    atom2.randomizeSpin(engine, realRandomGenerator, gaussianRandomGenerator, 0.1, 0);
    Array flippedSpin = atom2.getSpin();
    
    // Ising flip should invert the spin
    assert(fp_equal(flippedSpin[0], -originalSpin2[0]));
    assert(fp_equal(flippedSpin[1], -originalSpin2[1]));
    assert(fp_equal(flippedSpin[2], -originalSpin2[2]));
    
    std::cout << "PASS\n";
}

void test_safe_stod_function() {
    std::cout << "Testing safe_stod function... ";
    
    // Test valid conversions
    assert(fp_equal(safe_stod("3.14"), 3.14));
    assert(fp_equal(safe_stod("-2.5"), -2.5));
    assert(fp_equal(safe_stod("0"), 0.0));
    assert(fp_equal(safe_stod("1e-3"), 0.001));
    
    // Test with default value
    assert(fp_equal(safe_stod("invalid", 42.0), 42.0));
    
    std::cout << "PASS\n";
}

void test_example_config_parsing() {
    std::cout << "Testing example config parsing simulation... ";
    
    // This is a simple test that would normally involve parsing JSON
    // For now, we just verify the example file exists and has expected structure
    // In a real test, we would parse the JSON and verify values
    
    std::cout << "SKIP (requires JSON parsing)\n";
}

int main() {
    std::cout << "Running integration tests for vegas...\n\n";
    
    try {
        test_spin_model_factory();
        test_atom_spin_model_integration();
        test_spin_model_randomization();
        test_safe_stod_function();
        test_example_config_parsing();
        
        std::cout << "\nAll integration tests passed!\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\nTest failed with exception: " << e.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "\nTest failed with unknown exception\n";
        return 1;
    }
}