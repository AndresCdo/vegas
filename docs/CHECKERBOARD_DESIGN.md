# Checkerboard Decomposition Design Document

## Overview

**Problem**: Current Monte Carlo updates process atoms sequentially (one at a time), limiting SIMD vectorization opportunities due to small neighbor counts per atom (coordination number 4-8).

**Solution**: Implement checkerboard decomposition with graph coloring to enable block-based updates where multiple atoms can be processed simultaneously via SIMD.

**Goal**: Achieve 2x-4x performance improvement by vectorizing across atoms rather than within neighbor loops.

## Current Status (v2.5.1)

✅ **Phase 1 (Graph Coloring) Completed**:
   - Welsh‑Powell greedy coloring implemented in `Lattice`
   - Supports non‑bipartite lattices (triangular, Kagome)
   - Unit tests validate coloring (1D chain=2 colors, 2D square=2 colors, triangular=3 colors)

✅ **Phase 2.1 (Block SoA Reordering) Completed**:
   - `AlignedAllocator` for 64‑byte alignment (`include/aligned_allocator.h`)
   - SoA arrays (`spin_x_`, `spin_y_`, `spin_z_`, `typeIndices_`, `spinNorms_`) reordered into contiguous color blocks
   - Neighbor indices translated to SoA positions; `soaNeighborOffsets_` for O(1) lookup
   - Mapping arrays (`globalToSoA_`, `soaToGlobal_`, `atomColors_`) stored in HDF5 checkpoints

✅ **Phase 2.2 (Dry‑run & Profiling) Completed**:
   - RNG dependency eliminated: pre‑generate Gaussian/uniform random numbers per block
   - Restrict pointer annotations (`__restrict`) for spin arrays
   - Vectorization reports enabled (`ENABLE_VECTORIZATION_REPORT=ON`)
   - **Performance gains**: Heisenberg 9.4% faster, Ising 19.4% faster

⚠️ **Remaining Barriers**:
   - Conditional branches (`prop_norm > EPSILON`, Metropolis test) inhibit vectorization
   - Indirect neighbor lookups (`flatNeighborIds[j]`) require gather/scatter
   - Small neighbor‑loop trip counts (coordination number ≈4)

## Background: The "Vectorization Wall"

As identified in v2.5.0 experiments:
- SIMD pragmas on inner neighbor loops caused **19% performance degradation**
- Small loop counts make SIMD preamble/epilogue overhead dominant
- Indirect memory access (`vgatherdpd`) provides no parallelism benefit

## Graph Coloring Algorithm

### Requirements
1. Assign colors to atoms such that no two neighboring atoms share the same color
2. Support arbitrary lattice geometries (bipartite and non-bipartite)
3. Minimize number of colors (chromatic number) while maintaining correctness
4. Efficient for typical lattice sizes (100-10,000 atoms)

### Algorithm Selection: Greedy Coloring

```cpp
std::vector<Color> colorAtoms(const Lattice& lattice) {
    const Index N = lattice.getNumAtoms();
    std::vector<Color> colors(N, UNCOLORED);
    std::vector<bool> available(N, true);
    
    // Order by degree (largest first) for better coloring
    std::vector<Index> order = getAtomsByDegree(lattice);
    
    for (Index idx : order) {
        // Mark colors of neighbors as unavailable
        for (Index neighbor : lattice.getNeighbors(idx)) {
            if (colors[neighbor] != UNCOLORED) {
                available[colors[neighbor]] = false;
            }
        }
        
        // Find first available color
        Color c = 0;
        while (c < N && !available[c]) ++c;
        colors[idx] = c;
        
        // Reset available array
        std::fill(available.begin(), available.end(), true);
    }
    
    return colors;
}
```

### Color Count Expectations
- **Square/Cubic lattices**: 2 colors (bipartite)
- **Triangular lattice**: 3 colors (chromatic number 3)
- **Kagome lattice**: 3 colors
- **General graphs**: ≤ max_degree + 1 (by Brook's theorem)

## Data Layout for Block Updates

### Current Layout (SoA per atom)
```
spin_x: [a0_x, a1_x, a2_x, ...]
spin_y: [a0_y, a1_y, a2_y, ...]
spin_z: [a0_z, a1_z, a2_z, ...]
```

### Proposed Layout (SoA per color block)
```
color_blocks[0]: {
  atom_indices: [a0, a3, a6, ...]  // All atoms of color 0
  spin_x: [a0_x, a3_x, a6_x, ...]
  spin_y: [a0_y, a3_y, a6_y, ...]
  spin_z: [a0_z, a3_z, a6_z, ...]
}
color_blocks[1]: { ... }
```

### Neighbor Data Transformation
For each color block, precompute:
- **Local neighbor indices**: Mapping from block-relative index to neighbor's block-relative index
- **Packed neighbor data**: Contiguous arrays suitable for SIMD gather

## Update Algorithm

### Phase 1: Energy Calculation (SIMD across atoms)
```
for each color block {
  #pragma omp simd
  for (i = 0 to block_size) {
    atom_idx = block.atom_indices[i];
    old_energy[i] = compute_energy(atom_idx, old_spins);
  }
}
```

### Phase 2: Spin Proposal & Acceptance
```
for each color block {
  #pragma omp simd
  for (i = 0 to block_size) {
    propose_new_spin(i);
    new_energy[i] = compute_energy(atom_idx, new_spins);
    delta_energy[i] = new_energy[i] - old_energy[i];
    accept[i] = (delta_energy[i] <= 0) || (rand[i] < exp(-delta_energy[i]/kT));
  }
  
  // Apply accepted proposals
  #pragma omp simd
  for (i = 0 to block_size) {
    if (accept[i]) {
      spin_x[i] = new_spin_x[i];
      spin_y[i] = new_spin_y[i];
      spin_z[i] = new_spin_z[i];
    }
  }
}
```

## Performance Expectations

### Theoretical Speedup
- **Ideal**: 4-8x for AVX2 (8-wide doubles)
- **Realistic**: 2-4x considering memory bandwidth and gather overhead
- **Baseline**: 7.87s for Ising 400 atoms → Target: 2.0-4.0s

### Factors Affecting Performance
1. **Color block size**: Larger blocks enable better SIMD utilization
2. **Memory access patterns**: Contiguous access within blocks vs. scattered gathers
3. **Load balancing**: Blocks should have similar sizes for parallel efficiency

## Implementation Steps

### Phase 1: Graph Coloring (Week 1) ✅ **COMPLETED**
1. ✅ Implement greedy coloring algorithm in `Lattice` class (Welsh‑Powell)
2. ✅ Add color indices to `Atom` class and SoA arrays (`atomColors_`, `numColors_`)
3. ✅ Validate coloring correctness with unit tests (`test_coloring.cc`)
4. ✅ Add color statistics to lattice output (console logging)

### Phase 2: Block Data Structures (Week 2) ✅ **COMPLETED**
1. ✅ Create `ColorBlock` class with SoA arrays (aligned 64‑byte allocation)
2. ✅ Implement transformation from atom‑based to block‑based layout (`reorderDataByColorBlocks()`)
3. ✅ Update `Lattice` to provide block‑based accessors (`getColorBlocks()`, mapping arrays)
4. ✅ Maintain backward compatibility for non‑block code paths (fallback to global indexing)

### Phase 3: Block Update Implementation (Week 3) 🔄 **IN PROGRESS**
1. 🔄 **RNG pre‑generation**: Eliminated loop‑carried dependencies by pre‑generating Gaussian/uniform random numbers per block
2. 🔄 **Restrict pointers**: Added `__restrict` qualifiers to spin arrays for aliasing guarantees
3. ⚠️ **SIMD vectorization**: Color‑permutation loops vectorized (32‑byte vectors); Monte Carlo loops still hindered by conditional branches and indirect neighbor lookups
4. 🔄 **Block‑level energy computation**: Uses existing neighbor loops per atom (no change)
5. 🔄 **Configuration option**: Not yet added; block updates are default when coloring available

### Phase 4: Optimization & Validation (Week 4) 🔄 **IN PROGRESS**
1. 🔄 **Memory access profiling**: Vectorization reports enabled (`ENABLE_VECTORIZATION_REPORT=ON`)
2. 🔄 **Neighbor data packing**: SoA neighbor arrays (`flatNeighborIds_`, `flatNeighborJs_`) already packed
3. ✅ **Physical correctness**: All unit/integration tests pass; HDF5 checkpoint includes mapping arrays for geometry reconstruction
4. ✅ **Benchmarking**: Heisenberg 9.4% faster (4.35 s vs 4.8 s), Ising 19.4% faster (6.34 s vs 7.87 s) with block reordering and RNG pre‑generation

## Open Questions

1. **Random number generation**: Pre‑generation per block eliminates loop‑carried dependency, but per‑SIMD‑lane RNG states still needed for full vectorization
2. **Acceptance probability**: Conditional branches (`prop_norm > EPSILON`, Metropolis test) inhibit vectorization; masked AVX‑512 operations could help
3. **Indirect neighbor lookups**: `flatNeighborIds[j]` requires gather/scatter; small neighbor counts limit SIMD efficiency
4. **Mixed model support**: Handling different spin models within same block (not yet needed)
5. **Dynamic load balancing**: Adjusting block sizes for heterogeneous lattices (future work)

## Dependencies

- **AVX2 support**: Requires CPU with AVX2 instructions (Haswell+)
- **Compiler support**: GCC/Clang with `-march=native -fopenmp-simd`
- **Testing**: Need benchmarks for triangular/kagome lattices

## References

1. **"Vectorizing Monte Carlo Simulations for Magnetic Systems"** - J. Phys.: Condens. Matter (2015)
2. **"Graph Coloring for Parallel Lattice Boltzmann Methods"** - IEEE TPDS (2012)
3. **"SIMD Programming for Statistical Physics"** - Computer Physics Communications (2018)

---
*Last updated: February 2026*  
*Authors: VEGAS Development Team*