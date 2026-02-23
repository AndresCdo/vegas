#include "../include/lattice.h"
#include "../include/exception.h"

#include <algorithm>


Lattice::Lattice(std::string fileName)
{
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
    
    // Fill flat arrays
    Index currentOffset = 0;
    for (Index i = 0; i < num_ions; ++i) {
        neighborOffsets_[i] = currentOffset;
        for (size_t j = 0; j < neighborIndexes_[i].size(); ++j) {
            flatNeighbors_.push_back({neighborIndexes_[i][j], exchanges_[i][j]});
            currentOffset++;
        }
    }
    neighborOffsets_[num_ions] = currentOffset;  // Sentinel

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
      neighborIndexes_(std::move(other.neighborIndexes_)),
      exchanges_(std::move(other.exchanges_)),
      typeIndices_(std::move(other.typeIndices_)),
      spinNorms_(std::move(other.spinNorms_)),
      neighborOffsets_(std::move(other.neighborOffsets_)),
      flatNeighbors_(std::move(other.flatNeighbors_))
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
        neighborIndexes_ = std::move(other.neighborIndexes_);
        exchanges_ = std::move(other.exchanges_);
        typeIndices_ = std::move(other.typeIndices_);
        spinNorms_ = std::move(other.spinNorms_);
        neighborOffsets_ = std::move(other.neighborOffsets_);
        flatNeighbors_ = std::move(other.flatNeighbors_);
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
    for (Index i = 0; i < N; ++i) {
        const Array& spin = atoms_[i].getSpin();
        spin_x_[i] = spin[0];
        spin_y_[i] = spin[1];
        spin_z_[i] = spin[2];
        old_x_[i] = spin[0];
        old_y_[i] = spin[1];
        old_z_[i] = spin[2];
    }
}

// Sync from SoA to Atoms (call before HDF5 output)
void Lattice::syncSoAToAtoms()
{
    Index N = atoms_.size();
    for (Index i = 0; i < N; ++i) {
        Array spin({spin_x_[i], spin_y_[i], spin_z_[i]});
        atoms_[i].setSpin(spin);
        Array oldSpin({old_x_[i], old_y_[i], old_z_[i]});
        atoms_[i].setOldSpin(oldSpin);
    }
}
