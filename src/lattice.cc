#include "../include/lattice.h"
#include "../include/exception.h"
#include "../include/starter.h"

#include <algorithm>
#include <iostream>
#include <cstdint>


Lattice::Lattice(std::string fileName)
{
    STARTER::CHECKFILE(fileName);
    std::ifstream file(fileName);
    if (!file.is_open()) {
        throw vegas::FileIOException("Cannot open lattice file: " + fileName);
    }

    Index num_ions = 0;
    Index num_interactions = 0;
    Index num_types = 0;

    if (!(file >> num_ions >> num_interactions >> num_types)) {
        throw vegas::InvalidInputException("Invalid lattice file format (header): " + fileName);
    }

    if (num_ions == 0) {
        throw vegas::InvalidInputException("Lattice file has zero atoms: " + fileName);
    }
    if (num_types == 0) {
        throw vegas::InvalidInputException("Lattice file has zero atom types: " + fileName);
    }

    for (Index i = 0; i < num_types; ++i)
    {
        std::string type;
        file >> type;
        this -> mapTypeIndexes_[type] = i;
        this -> mapIndexTypes_[i] = type;
    }

    this -> sizesByIndex_ = std::vector<Index>(num_types);

    this -> atoms_ = std::vector<Atom>(num_ions);
    
    // Allocate SoA arrays
    this -> spin_x_.resize(num_ions);
    this -> spin_y_.resize(num_ions);
    this -> spin_z_.resize(num_ions);
    this -> old_x_.resize(num_ions);
    this -> old_y_.resize(num_ions);
    this -> old_z_.resize(num_ions);
    this -> typeIndices_.resize(num_ions);
    this -> spinNorms_.resize(num_ions);
    
    // Allocate interleaved SoA arrays (experimental)
    this -> spin_xyz_.resize(3 * num_ions);
    this -> old_xyz_.resize(3 * num_ions);
    
    // Reserve neighbor data
    this -> neighborIndexes_.resize(num_ions);
    this -> exchanges_.resize(num_ions);
    
    Real px;
    Real py;
    Real pz;
    Real spinNorm;
    Real hx;
    Real hy;
    Real hz;
    std::string type;
    std::string model;

    for (Index index = 0; index < num_ions; ++index)
    {
        if (!(file >> index >> px >> py >> pz >> spinNorm >> hx >> hy >> hz >> type >> model)) {
            throw vegas::InvalidInputException("Invalid lattice file format at atom " + std::to_string(index) + ": " + fileName);
        }

        if (spinNorm <= 0.0) {
            throw vegas::InvalidInputException("Invalid spin norm at atom " + std::to_string(index) + " in: " + fileName);
        }

        auto typeIt = this -> mapTypeIndexes_.find(type);
        if (typeIt == this -> mapTypeIndexes_.end()) {
            throw vegas::InvalidInputException("Unknown atom type '" + type + "' at atom " + std::to_string(index) + " in: " + fileName);
        }

        Array position({px, py, pz});

        Array spin({0.0, 0.0, spinNorm}); // ALWAYS THE INITIAL SPIN WILL BE IN THE Z-DIRECTION

        std::transform(model.begin(), model.end(), model.begin(), tolower);

        Atom atom(index, spin, position);
        atom.setType(type);
        atom.setExternalField({hx, hy, hz});
        atom.setModel(model);
        atom.setTypeIndex(this -> mapTypeIndexes_.at(type));

        this -> atoms_.at(index) = std::move(atom);


        this -> sizesByIndex_.at(this -> mapTypeIndexes_.at(type)) += 1;
        
        // Initialize SoA arrays
        spin_x_[index] = 0.0;
        spin_y_[index] = 0.0;
        spin_z_[index] = spinNorm;
        typeIndices_[index] = this -> mapTypeIndexes_.at(type);
        spinNorms_[index] = spinNorm;
    }

    Index index;
    Index nbh;
    Real exchange;
    for (Index _ = 0; _ < num_interactions; ++_)
    {
        if (!(file >> index >> nbh >> exchange)) {
            throw vegas::InvalidInputException("Invalid lattice file format at interaction " + std::to_string(_) + ": " + fileName);
        }
        if (index >= num_ions) {
            throw vegas::InvalidInputException("Atom index " + std::to_string(index) + " out of bounds in: " + fileName);
        }
        if (nbh >= num_ions) {
            throw vegas::InvalidInputException("Neighbor index " + std::to_string(nbh) + " out of bounds in: " + fileName);
        }
        this -> atoms_.at(index).addNeighborIndex(nbh);
        this -> atoms_.at(index).addExchange(exchange);
        
        // Also store in SoA-friendly format
        this -> neighborIndexes_[index].push_back(nbh);
        this -> exchanges_[index].push_back(exchange);
    }

    if (file.fail() && !file.eof()) {
        throw vegas::InvalidInputException("Error reading lattice file: " + fileName);
    }
    
    // Build flat neighbor arrays (Jagged Array with AoS)
    // Count total neighbors first
    Index totalNeighbors = 0;
    for (Index i = 0; i < num_ions; ++i) {
        totalNeighbors += neighborIndexes_[i].size();
    }
    
    // Allocate flat arrays
    neighborOffsets_.resize(num_ions + 1);
    flatNeighbors_.reserve(totalNeighbors);
    flatNeighborIds_.reserve(totalNeighbors);
    flatNeighborJs_.reserve(totalNeighbors);
    
    // Fill flat arrays
    Index currentOffset = 0;
    for (Index i = 0; i < num_ions; ++i) {
        neighborOffsets_[i] = currentOffset;
        for (size_t j = 0; j < neighborIndexes_[i].size(); ++j) {
            Index nid = neighborIndexes_[i][j];
            Real Jval = exchanges_[i][j];
            flatNeighbors_.push_back({nid, Jval});
            flatNeighborIds_.push_back(nid);
            flatNeighborJs_.push_back(Jval);
            currentOffset++;
        }
    }
    neighborOffsets_[num_ions] = currentOffset;  // Sentinel
    
    // Compute graph coloring for checkerboard decomposition
    computeGraphColoring();
    reorderDataByColorBlocks();
    verifyAlignment();
}

Lattice::~Lattice()
{

}

// Move constructor
Lattice::Lattice(Lattice&& other) noexcept
    : atoms_(std::move(other.atoms_)),
      mapTypeIndexes_(std::move(other.mapTypeIndexes_)),
      mapIndexTypes_(std::move(other.mapIndexTypes_)),
      sizesByIndex_(std::move(other.sizesByIndex_)),
       spin_x_(std::move(other.spin_x_)),
       spin_y_(std::move(other.spin_y_)),
       spin_z_(std::move(other.spin_z_)),
       old_x_(std::move(other.old_x_)),
       old_y_(std::move(other.old_y_)),
       old_z_(std::move(other.old_z_)),
       spin_xyz_(std::move(other.spin_xyz_)),
       old_xyz_(std::move(other.old_xyz_)),
      neighborIndexes_(std::move(other.neighborIndexes_)),
      exchanges_(std::move(other.exchanges_)),
      typeIndices_(std::move(other.typeIndices_)),
      spinNorms_(std::move(other.spinNorms_)),
       neighborOffsets_(std::move(other.neighborOffsets_)),
       soaNeighborOffsets_(std::move(other.soaNeighborOffsets_)),
       flatNeighbors_(std::move(other.flatNeighbors_)),
       flatNeighborIds_(std::move(other.flatNeighborIds_)),
       flatNeighborJs_(std::move(other.flatNeighborJs_)),
       atomColors_(std::move(other.atomColors_)),
       numColors_(other.numColors_),
       globalToSoA_(std::move(other.globalToSoA_)),
       soaToGlobal_(std::move(other.soaToGlobal_)),
       colorBlocks_(std::move(other.colorBlocks_))
{
}

// Move assignment operator
Lattice& Lattice::operator=(Lattice&& other) noexcept
{
    if (this != &other) {
        atoms_ = std::move(other.atoms_);
        mapTypeIndexes_ = std::move(other.mapTypeIndexes_);
        mapIndexTypes_ = std::move(other.mapIndexTypes_);
        sizesByIndex_ = std::move(other.sizesByIndex_);
        spin_x_ = std::move(other.spin_x_);
        spin_y_ = std::move(other.spin_y_);
        spin_z_ = std::move(other.spin_z_);
        old_x_ = std::move(other.old_x_);
        old_y_ = std::move(other.old_y_);
        old_z_ = std::move(other.old_z_);
        spin_xyz_ = std::move(other.spin_xyz_);
        old_xyz_ = std::move(other.old_xyz_);
        neighborIndexes_ = std::move(other.neighborIndexes_);
        exchanges_ = std::move(other.exchanges_);
        typeIndices_ = std::move(other.typeIndices_);
        spinNorms_ = std::move(other.spinNorms_);
        neighborOffsets_ = std::move(other.neighborOffsets_);
        soaNeighborOffsets_ = std::move(other.soaNeighborOffsets_);
        flatNeighbors_ = std::move(other.flatNeighbors_);
        flatNeighborIds_ = std::move(other.flatNeighborIds_);
        flatNeighborJs_ = std::move(other.flatNeighborJs_);
        atomColors_ = std::move(other.atomColors_);
        numColors_ = other.numColors_;
        globalToSoA_ = std::move(other.globalToSoA_);
        soaToGlobal_ = std::move(other.soaToGlobal_);
        colorBlocks_ = std::move(other.colorBlocks_);
    }
    return *this;
}

std::vector<Atom>& Lattice::getAtoms()
{
    return this -> atoms_;
}

const std::map<std::string, Index>& Lattice::getMapTypeIndexes() const
{
    return this -> mapTypeIndexes_;
}

const std::map<Index, std::string>& Lattice::getMapIndexTypes() const
{
    return this -> mapIndexTypes_;
}

const std::vector<Index>& Lattice::getSizesByIndex() const
{
    return this -> sizesByIndex_;
}

// Sync from Atoms to SoA (call after JSON parsing)
void Lattice::syncAtomsToSoA()
{
    Index N = atoms_.size();
    for (Index globalIdx = 0; globalIdx < N; ++globalIdx) {
        const Array& spin = atoms_[globalIdx].getSpin();
        Index soaIdx = globalToSoA_[globalIdx];
        spin_x_[soaIdx] = spin[0];
        spin_y_[soaIdx] = spin[1];
        spin_z_[soaIdx] = spin[2];
        old_x_[soaIdx] = spin[0];
        old_y_[soaIdx] = spin[1];
        old_z_[soaIdx] = spin[2];
        // Interleaved arrays
        spin_xyz_[3*soaIdx] = spin[0];
        spin_xyz_[3*soaIdx + 1] = spin[1];
        spin_xyz_[3*soaIdx + 2] = spin[2];
        old_xyz_[3*soaIdx] = spin[0];
        old_xyz_[3*soaIdx + 1] = spin[1];
        old_xyz_[3*soaIdx + 2] = spin[2];
    }
}

// Sync from SoA to Atoms (call before HDF5 output)
void Lattice::syncSoAToAtoms()
{
    Index N = atoms_.size();
    for (Index globalIdx = 0; globalIdx < N; ++globalIdx) {
        Index soaIdx = globalToSoA_[globalIdx];
        Array spin({spin_xyz_[3*soaIdx], spin_xyz_[3*soaIdx + 1], spin_xyz_[3*soaIdx + 2]});
        atoms_[globalIdx].setSpin(spin);
        Array oldSpin({old_xyz_[3*soaIdx], old_xyz_[3*soaIdx + 1], old_xyz_[3*soaIdx + 2]});
        atoms_[globalIdx].setOldSpin(oldSpin);
    }
}

// Copy separate arrays → interleaved arrays
void Lattice::syncSeparateToInterleaved()
{
    Index N = atoms_.size();
    for (Index soaIdx = 0; soaIdx < N; ++soaIdx) {
        spin_xyz_[3*soaIdx] = spin_x_[soaIdx];
        spin_xyz_[3*soaIdx + 1] = spin_y_[soaIdx];
        spin_xyz_[3*soaIdx + 2] = spin_z_[soaIdx];
        old_xyz_[3*soaIdx] = old_x_[soaIdx];
        old_xyz_[3*soaIdx + 1] = old_y_[soaIdx];
        old_xyz_[3*soaIdx + 2] = old_z_[soaIdx];
    }
}

// Copy interleaved arrays → separate arrays
void Lattice::syncInterleavedToSeparate()
{
    Index N = atoms_.size();
    for (Index soaIdx = 0; soaIdx < N; ++soaIdx) {
        spin_x_[soaIdx] = spin_xyz_[3*soaIdx];
        spin_y_[soaIdx] = spin_xyz_[3*soaIdx + 1];
        spin_z_[soaIdx] = spin_xyz_[3*soaIdx + 2];
        old_x_[soaIdx] = old_xyz_[3*soaIdx];
        old_y_[soaIdx] = old_xyz_[3*soaIdx + 1];
        old_z_[soaIdx] = old_xyz_[3*soaIdx + 2];
    }
}

// Graph coloring for checkerboard decomposition (Welsh-Powell greedy algorithm)
void Lattice::computeGraphColoring()
{
    Index N = atoms_.size();
    atomColors_.resize(N);
    std::fill(atomColors_.begin(), atomColors_.end(), N);  // Use N as uncolored sentinel
    
    // Compute degrees and create order (indices sorted by degree descending)
    std::vector<Index> order(N);
    std::vector<Index> degrees(N);
    for (Index i = 0; i < N; ++i) {
        order[i] = i;
        degrees[i] = neighborIndexes_[i].size();
    }
    
    // Sort indices by degree descending (Welsh-Powell ordering)
    std::stable_sort(order.begin(), order.end(),
        [&degrees](Index a, Index b) { return degrees[a] > degrees[b]; });
    
    // Boolean array for available colors (size N+1 to avoid resizing)
    std::vector<bool> available(N, true);
    
    // Greedy coloring
    for (Index idx : order) {
        // Mark colors of already-colored neighbors as unavailable
        for (Index nbh : neighborIndexes_[idx]) {
            Index neighborColor = atomColors_[nbh];
            if (neighborColor < N) {  // Neighbor is colored
                available[neighborColor] = false;
            }
        }
        
        // Find first available color
        Index color = 0;
        while (color < N && !available[color]) ++color;
        
        // Assign color to atom
        atomColors_[idx] = color;
        
        // Reset available array for next atom
        std::fill(available.begin(), available.end(), true);
    }
    
    // Compute number of colors used
    numColors_ = 0;
    for (Index color : atomColors_) {
        if (color + 1 > numColors_) {
            numColors_ = color + 1;
        }
    }
    
    // Output coloring statistics (debug)
    std::cout << "[Lattice] Graph coloring complete: Found " << numColors_ 
              << " colors for " << N << " atoms." << std::endl;
}

void Lattice::reorderDataByColorBlocks()
{
    Index N = atoms_.size();
    globalToSoA_.resize(N);
    soaToGlobal_.resize(N);
    colorBlocks_.clear();
    colorBlocks_.resize(numColors_);
    
    // Initialize each block with its color index
    for (Index c = 0; c < numColors_; ++c) {
        colorBlocks_[c] = ColorBlock(c);
    }
    
    // Assign SoA positions in color order
    Index soaPos = 0;
    for (Index color = 0; color < numColors_; ++color) {
        for (Index globalIdx = 0; globalIdx < N; ++globalIdx) {
            if (atomColors_[globalIdx] == color) {
                globalToSoA_[globalIdx] = soaPos;
                soaToGlobal_[soaPos] = globalIdx;
                colorBlocks_[color].addAtom(globalIdx);
                soaPos++;
            }
        }
    }
    if (soaPos != N) {
        throw vegas::SimulationException("Color block reordering mismatch: " + 
                                         std::to_string(soaPos) + " != " + std::to_string(N));
    }
    
    // DEBUG: Output block sizes
    std::cout << "[Lattice] Color block sizes:";
    for (Index c = 0; c < numColors_; ++c) {
        std::cout << " " << colorBlocks_[c].getSize();
    }
    std::cout << std::endl;
    
    // --- Physically reorder SoA data into contiguous color blocks ---
    
    // 1. Reorder spin arrays
    AlignedRealVector new_spin_x(N), new_spin_y(N), new_spin_z(N);
    AlignedRealVector new_old_x(N), new_old_y(N), new_old_z(N);
    AlignedRealVector new_spin_xyz(3*N), new_old_xyz(3*N);
    std::vector<Index> new_typeIndices(N);
    std::vector<Real> new_spinNorms(N);
    
    for (Index soaIdx = 0; soaIdx < N; ++soaIdx) {
        Index globalIdx = soaToGlobal_[soaIdx];
        new_spin_x[soaIdx] = spin_x_[globalIdx];
        new_spin_y[soaIdx] = spin_y_[globalIdx];
        new_spin_z[soaIdx] = spin_z_[globalIdx];
        new_old_x[soaIdx] = old_x_[globalIdx];
        new_old_y[soaIdx] = old_y_[globalIdx];
        new_old_z[soaIdx] = old_z_[globalIdx];
        // Interleaved arrays
        new_spin_xyz[3*soaIdx] = spin_xyz_[3*globalIdx];
        new_spin_xyz[3*soaIdx + 1] = spin_xyz_[3*globalIdx + 1];
        new_spin_xyz[3*soaIdx + 2] = spin_xyz_[3*globalIdx + 2];
        new_old_xyz[3*soaIdx] = old_xyz_[3*globalIdx];
        new_old_xyz[3*soaIdx + 1] = old_xyz_[3*globalIdx + 1];
        new_old_xyz[3*soaIdx + 2] = old_xyz_[3*globalIdx + 2];
        new_typeIndices[soaIdx] = typeIndices_[globalIdx];
        new_spinNorms[soaIdx] = spinNorms_[globalIdx];
    }
    
    // Swap with original arrays
    spin_x_.swap(new_spin_x);
    spin_y_.swap(new_spin_y);
    spin_z_.swap(new_spin_z);
    old_x_.swap(new_old_x);
    old_y_.swap(new_old_y);
    old_z_.swap(new_old_z);
    spin_xyz_.swap(new_spin_xyz);
    old_xyz_.swap(new_old_xyz);
    typeIndices_.swap(new_typeIndices);
    spinNorms_.swap(new_spinNorms);
    
    // 2. Update neighbor indices to SoA indexing
    for (size_t i = 0; i < flatNeighborIds_.size(); ++i) {
        Index globalNeighborId = flatNeighborIds_[i];
        flatNeighborIds_[i] = globalToSoA_[globalNeighborId];
    }
    
    // Also update flatNeighbors_ AoS structure
    for (auto& neighbor : flatNeighbors_) {
        neighbor.id = globalToSoA_[neighbor.id];
    }
    
    // 3. Build soaNeighborOffsets_ (prefix sum of neighbor counts in SoA order)
    soaNeighborOffsets_.resize(N + 1);
    soaNeighborOffsets_[0] = 0;
    for (Index soaIdx = 0; soaIdx < N; ++soaIdx) {
        Index globalIdx = soaToGlobal_[soaIdx];
        Index neighStart = neighborOffsets_[globalIdx];
        Index neighEnd = neighborOffsets_[globalIdx + 1];
        Index neighborCount = neighEnd - neighStart;
        soaNeighborOffsets_[soaIdx + 1] = soaNeighborOffsets_[soaIdx] + neighborCount;
    }
    
    // Verify total neighbor count matches
    if (soaNeighborOffsets_[N] != flatNeighborIds_.size()) {
        throw vegas::SimulationException("Neighbor count mismatch after SoA reordering");
    }
    
    std::cout << "[Lattice] SoA data reordered into " << numColors_ << " color blocks." << std::endl;
}

// Verify that no two neighboring atoms share the same color
bool Lattice::verifyColoring() const
{
    Index N = atoms_.size();
    for (Index i = 0; i < N; ++i) {
        Index color_i = atomColors_[i];
        if (color_i >= N) {
            std::cerr << "[Lattice] Atom " << i << " is uncolored!" << std::endl;
            return false;
        }
        for (Index nbh : neighborIndexes_[i]) {
            if (atomColors_[nbh] == color_i) {
                std::cerr << "[Lattice] Coloring violation: Atoms " << i 
                          << " and " << nbh << " share color " << color_i << std::endl;
                return false;
            }
        }
    }
    return true;
}

void Lattice::verifyAlignment() const {
#ifndef NDEBUG
    auto check = [](const void* ptr, const char* name) {
        uintptr_t addr = reinterpret_cast<uintptr_t>(ptr);
        if (addr % SIMD_ALIGNMENT != 0) {
            std::cerr << "[Lattice] WARNING: " << name << " is not aligned to " 
                      << SIMD_ALIGNMENT << " bytes (address: " << addr << ")" << std::endl;
        }
    };
    
    check(spin_x_.data(), "spin_x");
    check(spin_y_.data(), "spin_y");
    check(spin_z_.data(), "spin_z");
    check(old_x_.data(), "old_x");
    check(old_y_.data(), "old_y");
    check(old_z_.data(), "old_z");
    check(spin_xyz_.data(), "spin_xyz");
    check(old_xyz_.data(), "old_xyz");
    check(flatNeighborIds_.data(), "flatNeighborIds");
    check(flatNeighborJs_.data(), "flatNeighborJs");
#endif
}


