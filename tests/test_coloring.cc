// Graph coloring unit tests for checkerboard decomposition
// Tests Welsh-Powell greedy coloring on different lattice topologies

#include "../include/params.h"
#include "../include/lattice.h"
#include <iostream>
#include <cassert>

void test_1d_chain_coloring() {
    std::cout << "Testing 1D chain coloring (bipartite, 2 colors)... ";
    
    // Load 1D chain lattice with bidirectional bonds (100 atoms)
    Lattice lattice("../benchmarks/lattices/ising_chain_100.txt");
    
    const auto& colors = lattice.getAtomColors();
    Index numColors = lattice.getNumColors();
    
    assert(colors.size() == 100);
    assert(numColors == 2);  // 1D chain is bipartite
    
    // Verify coloring
    assert(lattice.verifyColoring());
    
    // Check that neighboring atoms have different colors
    // For a chain, each atom should have opposite color to its neighbors
    for (Index i = 0; i < 100; ++i) {
        const auto& neighbors = lattice.getNeighborIndexes()[i];
        for (Index neighbor : neighbors) {
            assert(colors[i] != colors[neighbor]);
        }
    }
    
    std::cout << "PASS\n";
}

void test_2d_square_coloring() {
    std::cout << "Testing 2D square lattice coloring (bipartite, 2 colors)... ";
    
    // Load 2D square lattice (20x20 from benchmarks)
    Lattice lattice("../benchmarks/lattices/ising_square_20x20_pbc.txt");
    
    const auto& colors = lattice.getAtomColors();
    Index numColors = lattice.getNumColors();
    
    assert(colors.size() == 400);  // 20x20 = 400 atoms
    assert(numColors == 2);  // Square lattice is bipartite
    
    // Verify coloring
    assert(lattice.verifyColoring());
    
    // Quick sanity: count atoms per color
    Index colorCount[3] = {0, 0, 0};
    for (Index color : colors) {
        assert(color >= 0 && color < 3);
        colorCount[color]++;
    }
    
    // Square lattice should have roughly equal number of atoms per color
    // (exactly equal for even LxL)
    assert(colorCount[0] == 200);
    assert(colorCount[1] == 200);
    assert(colorCount[2] == 0);
    
    std::cout << "PASS\n";
}

void test_triangular_coloring() {
    std::cout << "Testing triangular lattice coloring (non-bipartite, 3 colors)... ";
    
    // Load triangular lattice (3x3 with PBC)
    Lattice lattice("../test_system/triangular_3x3_pbc.txt");
    
    const auto& colors = lattice.getAtomColors();
    Index numColors = lattice.getNumColors();
    
    assert(colors.size() == 9);
    assert(numColors == 3);  // Triangular lattice requires 3 colors
    
    // Verify coloring
    assert(lattice.verifyColoring());
    
    // Count atoms per color
    Index colorCount[4] = {0, 0, 0, 0};
    for (Index color : colors) {
        assert(color >= 0 && color < 4);
        colorCount[color]++;
    }
    
    // Triangular 3x3 should have 3 atoms of each color (balanced)
    assert(colorCount[0] == 3);
    assert(colorCount[1] == 3);
    assert(colorCount[2] == 3);
    assert(colorCount[3] == 0);
    
    std::cout << "PASS\n";
}

void test_coloring_correctness() {
    std::cout << "Testing coloring correctness (no adjacent same color)... ";
    
    // Test on all three lattices
    {
        Lattice lattice1("../benchmarks/lattices/ising_chain_100.txt");
        assert(lattice1.verifyColoring());
    }
    {
        Lattice lattice2("../benchmarks/lattices/ising_square_20x20_pbc.txt");
        assert(lattice2.verifyColoring());
    }
    {
        Lattice lattice3("../test_system/triangular_3x3_pbc.txt");
        assert(lattice3.verifyColoring());
    }
    
    std::cout << "PASS\n";
}

void test_3x3_square_pbc_coloring() {
    std::cout << "Testing 3x3 square PBC lattice coloring (non-bipartite, 3 colors)... ";
    
    // Load 3x3 square lattice with PBC (odd-length cycle)
    Lattice lattice("../test_system/square_3x3_pbc.txt");
    
    const auto& colors = lattice.getAtomColors();
    Index numColors = lattice.getNumColors();
    
    assert(colors.size() == 9);
    assert(numColors == 3);  // 3x3 square with PBC has odd cycle, requires 3 colors
    
    // Verify coloring
    assert(lattice.verifyColoring());
    
    // Count atoms per color
    Index colorCount[4] = {0, 0, 0, 0};
    for (Index color : colors) {
        assert(color >= 0 && color < 4);
        colorCount[color]++;
    }
    
    // Should have 3 atoms per color (balanced)
    assert(colorCount[0] == 3);
    assert(colorCount[1] == 3);
    assert(colorCount[2] == 3);
    assert(colorCount[3] == 0);
    
    std::cout << "PASS\n";
}

void test_color_statistics_output() {
    std::cout << "Testing color statistics output... ";
    
    // This test just ensures the lattice constructor computes coloring
    // and doesn't crash. The actual output to std::cout is not captured.
    Lattice lattice("../benchmarks/lattices/ising_chain_100.txt");
    
    // Verify we have valid colors
    assert(lattice.getNumColors() > 0);
    assert(lattice.getAtomColors().size() == 100);
    
    std::cout << "PASS\n";
}

int main() {
    std::cout << "Running graph coloring tests for checkerboard decomposition...\n";
    
    try {
        test_1d_chain_coloring();
        test_2d_square_coloring();
        test_triangular_coloring();
        test_3x3_square_pbc_coloring();
        test_coloring_correctness();
        test_color_statistics_output();
        
        std::cout << "\nAll graph coloring tests passed!\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\nTest failed with exception: " << e.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "\nTest failed with unknown exception\n";
        return 1;
    }
}