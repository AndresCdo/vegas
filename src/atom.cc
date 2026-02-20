#include "../include/atom.h"

Real dot(const Array& A, const Array& B)
{
    return (A * B).sum();
}

Array cross(const Array& A, const Array& B)
{
    return Array{A[1]*B[2] - A[2]*B[1], A[2]*B[0] - A[0]*B[2], A[0]*B[1] - A[1]*B[0]};
}

Real norm(const Array& A)
{
    return std::sqrt((A*A).sum());
}


Atom::Atom() : Atom(0, ZERO, ZERO)
{

}

Atom::Atom(Index index, Array spin, Array position)
{

    this -> index_ = index;
    this -> spin_ = spin;
    this -> neighborIndexes_ = std::vector<Index>();
    this -> exchanges_ = std::vector<Real>();
    this -> position_ = position;
    this -> spinNorm_ = sqrt((this -> spin_ * this -> spin_).sum());
    this -> type_ = "nothing";

    this -> projections_ = std::vector<double>(0);
    this -> possibleProjections_ = std::vector<double>(0);
    for (double p = - this -> spinNorm_; p <= this -> spinNorm_; p += 1.0)
    {
        this -> projections_.push_back(p);
        this -> possibleProjections_.push_back(p);
    }

    this -> Sproj_ = 0;
    this -> setSpin({0.0, 0.0, this -> getPossibleProjections()[0]}); // ALWAYS THE INITIAL SPIN WILL BE IN THE Z-DIRECTION
    this -> removePossibleProjection(0);

    // Initialize with default spin model (Heisenberg)
    this -> model_ = "heisenberg";
    this -> spinModel_ = createSpinModel(this -> model_);

    // Initialize anisotropy members with defaults
    this -> anisotropyUnit_ = Array{0.0, 0.0, 1.0};  // Default unit vector along z
    this -> typeAnisotropy_ = "";
    this -> kan_ = 0.0;

    this -> oldSpin_ = this -> spin_;
}

Atom::~Atom()
{

}

// Move constructor
Atom::Atom(Atom&& other) noexcept
    : index_(std::move(other.index_))
    , spin_(std::move(other.spin_))
    , neighborIndexes_(std::move(other.neighborIndexes_))
    , exchanges_(std::move(other.exchanges_))
    , position_(std::move(other.position_))
    , spinNorm_(std::move(other.spinNorm_))
    , oldSpin_(std::move(other.oldSpin_))
    , externalField_(std::move(other.externalField_))
    , type_(std::move(other.type_))
    , typeIndex_(std::move(other.typeIndex_))
    , projections_(std::move(other.projections_))
    , possibleProjections_(std::move(other.possibleProjections_))
    , model_(std::move(other.model_))
    , spinModel_(std::move(other.spinModel_))
    , Sproj_(std::move(other.Sproj_))
    , anisotropyTerms_(std::move(other.anisotropyTerms_))
    , anisotropyUnit_(std::move(other.anisotropyUnit_))
    , typeAnisotropy_(std::move(other.typeAnisotropy_))
    , kan_(std::move(other.kan_))
{
}

// Move assignment operator
Atom& Atom::operator=(Atom&& other) noexcept
{
    if (this != &other) {
        index_ = std::move(other.index_);
        spin_ = std::move(other.spin_);
        neighborIndexes_ = std::move(other.neighborIndexes_);
        exchanges_ = std::move(other.exchanges_);
        position_ = std::move(other.position_);
        spinNorm_ = std::move(other.spinNorm_);
        oldSpin_ = std::move(other.oldSpin_);
        externalField_ = std::move(other.externalField_);
        type_ = std::move(other.type_);
        typeIndex_ = std::move(other.typeIndex_);
        projections_ = std::move(other.projections_);
        possibleProjections_ = std::move(other.possibleProjections_);
        model_ = std::move(other.model_);
        spinModel_ = std::move(other.spinModel_);
        Sproj_ = std::move(other.Sproj_);
        anisotropyTerms_ = std::move(other.anisotropyTerms_);
        anisotropyUnit_ = std::move(other.anisotropyUnit_);
        typeAnisotropy_ = std::move(other.typeAnisotropy_);
        kan_ = std::move(other.kan_);
    }
    return *this;
}

const Array& Atom::getPosition() const
{
    return this -> position_;
}

const Index& Atom::getIndex() const
{
    return this -> index_;
}

const std::vector<Index>& Atom::getNeighborIndexes() const
{
    return this -> neighborIndexes_;
}

const Real& Atom::getSpinNorm() const
{
    return this -> spinNorm_;
}

const Array& Atom::getSpin() const
{
    return this -> spin_;
}

const std::vector<Real>& Atom::getExchanges() const
{
    return this -> exchanges_;
}

const std::string& Atom::getType() const
{
    return this -> type_;
}

const std::vector<double>& Atom::getProjections() const
{
    return this -> projections_;
}

const std::vector<double>& Atom::getPossibleProjections() const
{
    return this -> possibleProjections_;
}

const Array& Atom::getExternalField() const
{
    return this -> externalField_;
}

void Atom::setPosition(const Array& position)
{
    this -> position_ = position;
}

void Atom::setNeighborIndexes(const std::vector<Index>& indices)
{
    this -> neighborIndexes_ = indices;
}

void Atom::setSpin(const Array& spin)
{
    this -> spin_ = spin;
    // this -> spinNorm_ = sqrt((this -> spin_ * this -> spin_).sum());
}

void Atom::setExchanges(const std::vector<Real>& exchanges)
{
    this -> exchanges_ = exchanges;
}

void Atom::setType(const std::string& type)
{
    this -> type_ = type;
}

void Atom::addNeighborIndex(Index neighborIndex)
{
    this -> neighborIndexes_.push_back(neighborIndex);
}

void Atom::addExchange(Real exchange)
{
    this -> exchanges_.push_back(exchange);
}

void Atom::setExternalField(const Array& externalField)
{
    this -> externalField_ = externalField;
}

void Atom::changeProjection(Index i, Real value)
{
    this -> possibleProjections_.at(i) = value;
}

void Atom::removePossibleProjection(Index i)
{
    this -> possibleProjections_.erase(this -> possibleProjections_.begin() + i);
}


const std::string& Atom::getModel() const
{
    return this -> model_;
}

void Atom::setModel(const std::string& model)
{
    this -> model_ = model;
    this -> spinModel_ = createSpinModel(model);
}

Real Atom::getExchangeEnergy(const std::vector<Atom>& atoms) const
{
    Real energy = 0.0;
    Index nbh_c = 0;
    for (auto&& nbhIndex : this -> neighborIndexes_)
        energy -= this -> exchanges_.at(nbh_c++) * (this -> spin_ * atoms.at(nbhIndex).getSpin()).sum();
    return energy;
}

Real Atom::getZeemanEnergy(const Real& H) const
{
    return - H * (this -> getSpin() * this -> getExternalField()).sum();
}

Real Atom::getAnisotropyEnergy(const Atom& atom) const
{
    Real anisotropyEnergy = 0.0;
    for (auto&& func : this -> anisotropyTerms_)
        anisotropyEnergy += func(atom);
    return anisotropyEnergy;
}


void Atom::setOldSpin(const Array& oldSpin)
{
    this -> oldSpin_ = oldSpin;
}

const Index& Atom::getSproj() const
{
    return this -> Sproj_;
}

void Atom::setSproj(const Index& Sproj)
{
    this -> Sproj_ = Sproj;
}

const Array& Atom::getOldSpin() const
{
    return this -> oldSpin_;
}

void Atom::revertSpin()
{
    this -> changeProjection(this -> Sproj_, this -> spin_[2]);
    this -> spin_ = this -> oldSpin_;
}

void Atom::addAnisotropyTerm(const std::function<Real(const Atom&)>& func)
{
    this -> anisotropyTerms_.push_back(func);
}

const Index& Atom::getTypeIndex() const
{
    return this -> typeIndex_;
}

void Atom::setTypeIndex(const Index& typeIndex)
{
    this -> typeIndex_ = typeIndex;
}

// Spin model method implementations
void Atom::initializeRandomState(
    std::mt19937_64& engine,
    std::uniform_real_distribution<>& realRandomGenerator,
    std::normal_distribution<>& gaussianRandomGenerator)
{
    if (this -> spinModel_) {
        this -> spinModel_->initializeRandomState(engine, realRandomGenerator, gaussianRandomGenerator, *this);
    }
}

void Atom::randomizeSpin(
    std::mt19937_64& engine,
    std::uniform_real_distribution<>& realRandomGenerator,
    std::normal_distribution<>& gaussianRandomGenerator,
    Real sigma,
    Index num)
{
    if (this -> spinModel_) {
        this -> spinModel_->randomizeSpin(engine, realRandomGenerator, gaussianRandomGenerator, sigma, *this, num);
    }
}

const Array& Atom::getAnisotropyUnit() const
{
    return this -> anisotropyUnit_;
}

const std::string& Atom::getTypeAnisotropy() const
{
    return this -> typeAnisotropy_;
}

const Real& Atom::getKan() const
{
    return this -> kan_;
}

void Atom::setAnisotropyUnit(const Array& unit)
{
    this -> anisotropyUnit_ = unit;
}

void Atom::setTypeAnisotropy(const std::string& type)
{
    this -> typeAnisotropy_ = type;
}

void Atom::setKan(Real kan)
{
    this -> kan_ = kan;
}
