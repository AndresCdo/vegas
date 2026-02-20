// System unit tests for VEGAS simulation package

#include "../include/system.h"
#include "../include/params.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdio>

void test_system_constructor()
{
    std::cout << "Testing System constructor... ";
    
    System system("test_system/ferromagnetic_chain.txt", {1.0}, {0.0}, 10, 42, "test_sys_output.h5", 1.0);
    
    assert(system.getSeed() == 42);
    assert(system.getLattice().getAtoms().size() == 2);
    
    std::remove("test_sys_output.h5");
    
    std::cout << "PASS\n";
}

void test_system_get_lattice()
{
    std::cout << "Testing System get lattice... ";
    
    System system("test_system/ferromagnetic_chain.txt", {1.0}, {0.0}, 10, 42, "test_lattice_output.h5", 1.0);
    
    Lattice& lattice = system.getLattice();
    assert(lattice.getAtoms().size() == 2);
    assert(lattice.getMapTypeIndexes().size() == 1);
    
    std::remove("test_lattice_output.h5");
    
    std::cout << "PASS\n";
}

void test_system_randomize_spins()
{
    std::cout << "Testing System randomize spins... ";
    
    System system("test_system/ferromagnetic_chain.txt", {1.0}, {0.0}, 10, 42, "test_random_output.h5", 1.0);
    system.randomizeSpins();
    
    auto& atoms = system.getLattice().getAtoms();
    for (auto& atom : atoms) {
        Real norm = std::sqrt((atom.getSpin() * atom.getSpin()).sum());
        assert(fp_equal(norm, atom.getSpinNorm()));
    }
    
    std::remove("test_random_output.h5");
    
    std::cout << "PASS\n";
}

void test_system_local_energy()
{
    std::cout << "Testing System local energy... ";
    
    System system("test_system/ferromagnetic_chain.txt", {1.0}, {0.0}, 10, 42, "test_energy_output.h5", 1.0);
    system.randomizeSpins();
    
    Real energy = system.localEnergy(0, 0.0);
    
    assert(std::isfinite(energy));
    
    std::remove("test_energy_output.h5");
    
    std::cout << "PASS\n";
}

void test_system_total_energy()
{
    std::cout << "Testing System total energy... ";
    
    System system("test_system/ferromagnetic_chain.txt", {1.0}, {0.0}, 10, 42, "test_total_energy_output.h5", 1.0);
    system.randomizeSpins();
    
    Real energy = system.totalEnergy(0.0);
    
    assert(std::isfinite(energy));
    
    std::remove("test_total_energy_output.h5");
    
    std::cout << "PASS\n";
}

void test_system_monte_carlo_step()
{
    std::cout << "Testing System Monte Carlo step... ";
    
    System system("test_system/ferromagnetic_chain.txt", {1.0}, {0.0}, 10, 42, "test_mc_output.h5", 1.0);
    system.randomizeSpins();
    
    Array initialSpin = system.getLattice().getAtoms()[0].getSpin();
    
    system.monteCarloStep(1.0, 0.0);
    
    std::remove("test_mc_output.h5");
    
    std::cout << "PASS\n";
}

void test_system_cycle_minimal()
{
    std::cout << "Testing System minimal cycle... ";
    
    System system("test_system/ferromagnetic_chain.txt", {1.0}, {0.0}, 10, 42, "test_cycle_output.h5", 1.0);
    system.randomizeSpins();
    
    system.cycle();
    
    std::remove("test_cycle_output.h5");
    
    std::cout << "PASS\n";
}

void test_system_set_state()
{
    std::cout << "Testing System set state... ";
    
    System system("test_system/ferromagnetic_chain.txt", {1.0}, {0.0}, 10, 42, "test_state_output.h5", 1.0);
    system.setState("test_system/initial_spins.txt");
    
    auto& atoms = system.getLattice().getAtoms();
    assert(fp_equal(atoms[0].getSpin()[2], 1.0));
    assert(fp_equal(atoms[1].getSpin()[2], 1.0));
    
    std::remove("test_state_output.h5");
    
    std::cout << "PASS\n";
}

int main()
{
    std::cout << "Running System tests for vegas...\n\n";
    
    try {
        test_system_constructor();
        test_system_get_lattice();
        test_system_randomize_spins();
        test_system_local_energy();
        test_system_total_energy();
        test_system_monte_carlo_step();
        test_system_cycle_minimal();
        test_system_set_state();
        
        std::cout << "\nAll System tests passed!\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\nTest failed with exception: " << e.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "\nTest failed with unknown exception\n";
        return 1;
    }
}
