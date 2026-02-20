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
    }

    if (file.fail() && !file.eof()) {
        throw vegas::InvalidInputException("Error reading lattice file: " + fileName);
    }

}

Lattice::~Lattice()
{

}

// Move constructor
Lattice::Lattice(Lattice&& other) noexcept
    : atoms_(std::move(other.atoms_)),
      mapTypeIndexes_(std::move(other.mapTypeIndexes_)),
      mapIndexTypes_(std::move(other.mapIndexTypes_)),
      sizesByIndex_(std::move(other.sizesByIndex_))
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
