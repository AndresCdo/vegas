#ifndef LATTICE
#define LATTICE

#include <vector>
#include <fstream>
#include <map>

#include "atom.h"

// Flat neighbor structure for SoA hot path
struct NeighborInteraction {
    Index id;
    Real J;
};

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
    const std::vector<NeighborInteraction>& getFlatNeighbors() const { return flatNeighbors_; }
    
    // Get type indices for all atoms (for SoA hot path)
    const std::vector<Index>& getTypeIndices() const { return typeIndices_; }
    
    // Get spin norms (for SoA hot path)
    const std::vector<Real>& getSpinNorms() const { return spinNorms_; }
    
    // Sync methods
    void syncAtomsToSoA();   // After JSON parsing: Atom → SoA
    void syncSoAToAtoms();   // Before HDF5: SoA → Atom (lazy)

private:
    std::vector<Atom> atoms_;
    std::map<std::string, Index> mapTypeIndexes_;
    std::map<Index, std::string> mapIndexTypes_;
    std::vector<Index> sizesByIndex_;

    // SoA arrays for spins - Source of Truth for hot path
    std::vector<Real> spin_x_, spin_y_, spin_z_;
    std::vector<Real> old_x_, old_y_, old_z_;
    
    // Type indices and spin norms for SoA hot path
    std::vector<Index> typeIndices_;
    std::vector<Real> spinNorms_;
    
    // Neighbor data (kept as-is for now)
    std::vector<std::vector<Index>> neighborIndexes_;
    std::vector<std::vector<Real>> exchanges_;
    
    // Flat neighbor data for SoA hot path (Jagged Array with AoS)
    std::vector<Index> neighborOffsets_;  // [N+1]
    std::vector<NeighborInteraction> flatNeighbors_;  // [total_neighbors]

};

#endif