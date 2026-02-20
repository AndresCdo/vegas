// Simple test file to verify basic functionality
// This doesn't require Google Test installation

#include "../include/params.h"
#include "../include/atom.h"
#include "../include/spin_model.h"
#include "../include/lattice.h"
#include <iostream>
#include <cassert>

void test_floating_point_comparison() {
    std::cout << "Testing floating point comparison... ";
    assert(fp_equal(1.0, 1.0));
    assert(fp_equal(1.0, 1.0 + EPSILON/2.0));
    assert(!fp_equal(1.0, 1.0 + EPSILON*2.0));
    assert(fp_not_equal(1.0, 1.0 + EPSILON*2.0));
    std::cout << "PASS\n";
}

void test_constants() {
    std::cout << "Testing constants... ";
    assert(DEFAULT_MCS == 5000);
    assert(fp_equal(DEFAULT_TEMP_START, 0.001));
    assert(fp_equal(DEFAULT_TEMP_FINAL, 10.0));
    assert(fp_equal(DEFAULT_TEMP_DELTA, 0.1));
    assert(fp_equal(MAX_SIGMA, 60.0));
    assert(fp_equal(MIN_SIGMA, 1e-10));
    assert(fp_equal(SIGMA_ADJUSTMENT_FACTOR, 0.5));
    assert(fp_equal(MIN_REJECTION_RATE, 1e-10));
    assert(NUM_SPIN_MODELS == 5);
    assert(fp_equal(SPIN_NORM_ROUNDING_PRECISION, 10000.0));
    std::cout << "PASS\n";
}

void test_array_operations() {
    std::cout << "Testing array operations... ";
    Array a = {1.0, 2.0, 3.0};
    Array b = {4.0, 5.0, 6.0};
    
    Array sum = a + b;
    assert(fp_equal(sum[0], 5.0));
    assert(fp_equal(sum[1], 7.0));
    assert(fp_equal(sum[2], 9.0));
    
    Array product = a * b;
    assert(fp_equal(product[0], 4.0));
    assert(fp_equal(product[1], 10.0));
    assert(fp_equal(product[2], 18.0));
    
    Real dot_product = (a * b).sum();
    assert(fp_equal(dot_product, 32.0));
    std::cout << "PASS\n";
}

void test_zero_constant() {
    std::cout << "Testing ZERO constant... ";
    assert(fp_equal(ZERO[0], 0.0));
    assert(fp_equal(ZERO[1], 0.0));
    assert(fp_equal(ZERO[2], 0.0));
    std::cout << "PASS\n";
}

void test_atom_basic() {
    std::cout << "Testing Atom basic operations... ";
    Array spin = {0.0, 0.0, 1.0};
    Array position = {0.0, 0.0, 0.0};
    Atom atom(0, spin, position);
    
    assert(atom.getIndex() == 0);
    assert(fp_equal(atom.getSpin()[0], 0.0));
    assert(fp_equal(atom.getSpin()[1], 0.0));
    // Spin may be initialized to -spinNorm (first possible projection)
    assert(fp_equal(std::abs(atom.getSpin()[2]), 1.0));
    assert(fp_equal(atom.getSpinNorm(), 1.0));
    assert(fp_equal(atom.getPosition()[0], 0.0));
    // Test setting external field
    Array field = {1.0, 0.0, 0.0};
    atom.setExternalField(field);
    assert(fp_equal(atom.getExternalField()[0], 1.0));
    // Test setting type
    atom.setType("Fe");
    assert(atom.getType() == "Fe");
    // Test setting model (Heisenberg)
    atom.setModel("heisenberg");
    assert(atom.getModel() == "heisenberg");
    
    std::cout << "PASS\n";
}

void test_lattice_loading() {
    std::cout << "Testing Lattice loading... ";
    // Use the test lattice file from test_system directory
    Lattice lattice("test_system/ferromagnetic_chain.txt");
    
    auto& atoms = lattice.getAtoms();
    assert(atoms.size() == 2);
    
    // Check atom indices and types
    assert(atoms[0].getIndex() == 0);
    assert(atoms[1].getIndex() == 1);
    assert(atoms[0].getType() == "Fe");
    assert(atoms[1].getType() == "Fe");
    
    // Check spin norms (should be 1.0)
    assert(fp_equal(atoms[0].getSpinNorm(), 1.0));
    assert(fp_equal(atoms[1].getSpinNorm(), 1.0));
    
    // Check positions
    assert(fp_equal(atoms[0].getPosition()[0], 0.0));
    assert(fp_equal(atoms[1].getPosition()[0], 1.0));
    
    // Check neighbor counts (only atom 0 has a neighbor)
    assert(atoms[0].getNeighborIndexes().size() == 1);
    assert(atoms[1].getNeighborIndexes().size() == 0);
    
    // Check exchanges
    assert(atoms[0].getExchanges().size() == 1);
    assert(fp_equal(atoms[0].getExchanges()[0], -1.0));
    assert(atoms[1].getExchanges().size() == 0);
    
    // Check type mappings
    const auto& typeIndexes = lattice.getMapTypeIndexes();
    assert(typeIndexes.size() == 1);
    assert(typeIndexes.at("Fe") == 0);
    
    const auto& indexTypes = lattice.getMapIndexTypes();
    assert(indexTypes.at(0) == "Fe");
    
    const auto& sizes = lattice.getSizesByIndex();
    assert(sizes.size() == 1);
    assert(sizes[0] == 2);
    
    std::cout << "PASS\n";
}

int main() {
    std::cout << "Running simple tests for vegas...\n";
    
    try {
        test_floating_point_comparison();
        test_constants();
        test_array_operations();
        test_zero_constant();
        test_atom_basic();
        test_lattice_loading();
        
        std::cout << "\nAll tests passed!\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\nTest failed with exception: " << e.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "\nTest failed with unknown exception\n";
        return 1;
    }
}