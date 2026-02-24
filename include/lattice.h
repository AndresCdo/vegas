#ifndef LATTICE
#define LATTICE

#include <vector>
#include <fstream>
#include <map>

#include "atom.h"
#include "color_block.h"
#include "aligned_allocator.h"
#include "params.h"

// Flat neighbor structure for SoA hot path
struct NeighborInteraction {
    Index id;
    Real J;
};

// Aligned vector types for SIMD
using AlignedRealVector = std::vector<Real, AlignedAllocator<Real, SIMD_ALIGNMENT>>;
using AlignedIndexVector = std::vector<Index, AlignedAllocator<Index, SIMD_ALIGNMENT>>;

class Lattice
{
public:
    Lattice(std::string fileName);
    ~Lattice();
    
    // Move constructor and assignment
    Lattice(Lattice&& other) noexcept;
    Lattice& operator=(Lattice&& other) noexcept;
    
    // Delete copy constructor and assignment
    Lattice(const Lattice&) = delete;
    Lattice& operator=(const Lattice&) = delete;

    std::vector<Atom>& getAtoms();

    const std::map<std::string, Index>& getMapTypeIndexes() const;
    const std::map<Index, std::string>& getMapIndexTypes() const;
    const std::vector<Index>& getSizesByIndex() const;

    // SoA accessors for hot path
    Index getNumAtoms() const { return atoms_.size(); }
    
    // Direct pointer access for hot path - returns raw pointers to SoA arrays
    // Interleaved layout access (experimental)
    const Real* getSpinXYZ() const { return spin_xyz_.data(); }
    Real* getSpinXYZ() { return spin_xyz_.data(); }
    const Real* getOldSpinXYZ() const { return old_xyz_.data(); }
    Real* getOldSpinXYZ() { return old_xyz_.data(); }
    
    // Legacy separate array access (deprecated, will be removed)
    const Real* getSpinX() const { return spin_x_.data(); }
    const Real* getSpinY() const { return spin_y_.data(); }
    const Real* getSpinZ() const { return spin_z_.data(); }
    Real* getSpinX() { return spin_x_.data(); }
    Real* getSpinY() { return spin_y_.data(); }
    Real* getSpinZ() { return spin_z_.data(); }
    
    const Real* getOldSpinX() const { return old_x_.data(); }
    const Real* getOldSpinY() const { return old_y_.data(); }
    const Real* getOldSpinZ() const { return old_z_.data(); }
    Real* getOldSpinX() { return old_x_.data(); }
    Real* getOldSpinY() { return old_y_.data(); }
    Real* getOldSpinZ() { return old_z_.data(); }
    
    // Get neighbor data (unchanged - still AoS for now)
    const std::vector<std::vector<Index>>& getNeighborIndexes() const { return neighborIndexes_; }
    const std::vector<std::vector<Real>>& getExchanges() const { return exchanges_; }
    
    // Flat neighbor accessors for SoA hot path (Jagged Array with AoS)
    const std::vector<Index>& getNeighborOffsets() const { return neighborOffsets_; }
    const std::vector<Index>& getSoANeighborOffsets() const { return soaNeighborOffsets_; }
    const std::vector<NeighborInteraction>& getFlatNeighbors() const { return flatNeighbors_; }
    
    // SoA neighbor accessors for SIMD vectorization
    const Index* getFlatNeighborIds() const { return flatNeighborIds_.data(); }
    const Real* getFlatNeighborJs() const { return flatNeighborJs_.data(); }
    
    // Get type indices for all atoms (for SoA hot path)
    const std::vector<Index>& getTypeIndices() const { return typeIndices_; }
    
    // Get spin norms (for SoA hot path)
    const std::vector<Real>& getSpinNorms() const { return spinNorms_; }
    
    // Sync methods
    void syncAtomsToSoA();   // After JSON parsing: Atom → SoA
    void syncSoAToAtoms();   // Before HDF5: SoA → Atom (lazy)
    void syncSeparateToInterleaved(); // Copy separate arrays → interleaved arrays
    void syncInterleavedToSeparate(); // Copy interleaved arrays → separate arrays
    
    // Graph coloring for checkerboard decomposition
    const std::vector<Index>& getAtomColors() const { return atomColors_; }
    Index getNumColors() const { return numColors_; }
    bool verifyColoring() const;
    
    // Color block access for checkerboard decomposition
    const std::vector<ColorBlock>& getColorBlocks() const { return colorBlocks_; }
    const std::vector<Index>& getGlobalToSoA() const { return globalToSoA_; }
    const std::vector<Index>& getSoAToGlobal() const { return soaToGlobal_; }

private:
    std::vector<Atom> atoms_;
    std::map<std::string, Index> mapTypeIndexes_;
    std::map<Index, std::string> mapIndexTypes_;
    std::vector<Index> sizesByIndex_;

    // SoA arrays for spins - Source of Truth for hot path
    // Interleaved layout: [x0,y0,z0, x1,y1,z1, ...] for cache locality (experimental)
    AlignedRealVector spin_xyz_;   // size = 3 * N
    AlignedRealVector old_xyz_;    // size = 3 * N, for rollback
    // Legacy separate arrays (kept for compatibility during transition)
    AlignedRealVector spin_x_, spin_y_, spin_z_;
    AlignedRealVector old_x_, old_y_, old_z_;
    
    // Type indices and spin norms for SoA hot path
    std::vector<Index> typeIndices_;
    std::vector<Real> spinNorms_;
    
    // Neighbor data (kept as-is for now)
    std::vector<std::vector<Index>> neighborIndexes_;
    std::vector<std::vector<Real>> exchanges_;
    
    // Flat neighbor data for SoA hot path (Jagged Array with AoS)
    std::vector<Index> neighborOffsets_;      // [N+1] indexed by global atom index (original ordering)
    std::vector<Index> soaNeighborOffsets_;   // [N+1] indexed by SoA position (color-block ordering)
    std::vector<NeighborInteraction> flatNeighbors_;  // [total_neighbors]
    
    // SoA neighbor arrays for SIMD vectorization
    AlignedIndexVector flatNeighborIds_;      // [total_neighbors] stores SoA indices of neighbors
    AlignedRealVector flatNeighborJs_;        // [total_neighbors]
    
    // Graph coloring for checkerboard decomposition
    std::vector<Index> atomColors_;           // [num_atoms]
    Index numColors_;                         // Number of colors used
    
    // Block reordering for checkerboard decomposition
    std::vector<Index> globalToSoA_;           // [num_atoms] maps global atom index to SoA position
    std::vector<Index> soaToGlobal_;           // [num_atoms] maps SoA position to global atom index
    std::vector<ColorBlock> colorBlocks_;      // One block per color
    
    // Private helper methods
    void computeGraphColoring();
    void reorderDataByColorBlocks();
    void verifyAlignment() const;

};

#endif