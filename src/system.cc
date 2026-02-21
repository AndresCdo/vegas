#include "../include/system.h"
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "../include/rlutil.h"
#include "../include/error.h"
#include <functional>
#include <stdexcept>

// Message to exit and launch an error.


// Safe conversion from string to double with error handling
Real safe_stod(const std::string& str, const std::string& context)
{
    try {
        size_t pos = 0;
        double value = std::stod(str, &pos);
        
        // Check if entire string was consumed
        if (pos != str.length()) {
            throw std::invalid_argument("Extra characters after number");
        }
        
        return value;
    } catch (const std::invalid_argument& e) {
        throw vegas::NumericConversionException("Invalid number format in " + context + ": '" + str + "'");
    } catch (const std::out_of_range& e) {
        throw vegas::NumericConversionException("Number out of range in " + context + ": '" + str + "'");
    }
    
    // Should never reach here
    return 0.0;
}

// Safe conversion from string to double with default value
Real safe_stod(const std::string& str, Real default_value)
{
    try {
        size_t pos = 0;
        double value = std::stod(str, &pos);
        
        // Check if entire string was consumed
        if (pos != str.length()) {
            return default_value;
        }
        
        return value;
    } catch (const std::invalid_argument&) {
        return default_value;
    } catch (const std::out_of_range&) {
        return default_value;
    }
}

std::vector<std::string> split(const std::string& s, char delim)
{
    std::stringstream ss(s);
    std::string item;
    std::vector<std::string> tokens;
    while (std::getline(ss, item, delim)) {
        tokens.push_back(item);
    }
    return tokens;
}

const std::string ETA_seconds(const Index& seconds)
{
    return "ETR: " + std::to_string(int(seconds / 3600)) + ":"
                   + std::to_string(int((seconds % 3600) / 60)) + ":"
                   + (((seconds % 60)< 10)? "0":"") + std::to_string(seconds % 60);
}

System::System(std::string fileName,
               std::vector<Real> temps,
               std::vector<Real> fields,
               Index mcs,
               Index seed,
               std::string outName,
               Real kb) : lattice_(fileName)
{
    this -> mcs_ = mcs;
    this -> thermalizationFraction_ = DEFAULT_THERMALIZATION_FRACTION;
    this -> thermalizationSteps_ = static_cast<Index>(mcs * this -> thermalizationFraction_);
    this -> measurementSteps_ = mcs - this -> thermalizationSteps_;
    this -> kb_ = kb;
    this -> temps_ = temps;
    this -> fields_ = fields;
    this -> intRandomGenerator_ = std::uniform_int_distribution<>(0, this -> lattice_.getAtoms().size() - 1);
    this -> gaussianRandomGenerator_ = std::normal_distribution<>(0.0, 1.0);
    this -> outName_ = outName;

    this -> num_types_ = this -> lattice_.getMapTypeIndexes().size();

    this -> engine_.seed(seed);
    this -> seed_ = seed;

    this -> sigma_ = std::vector<Real>(this -> num_types_);
    this -> counterRejections_ = std::vector<Index>(this -> num_types_);
    this -> magnetizationByTypeIndex_ = std::vector<Array>(this -> num_types_ + 1);
    for (Index i = 0; i < this -> sigma_.size(); ++i)
    {
        this -> sigma_.at(i) = MAX_SIGMA;
        this -> counterRejections_.at(i) = 0;
        this -> magnetizationByTypeIndex_.at(i) = ZERO;
    }
    this -> magnetizationByTypeIndex_.at(this -> num_types_) = ZERO; // for total magnetization

    std::remove(outName.c_str());


    this -> reporter_ = Reporter(this -> outName_,
                                 this -> magnetizationByTypeIndex_,
                                 this -> lattice_,
                                 this -> temps_,
                                 this -> fields_,
                                 this -> measurementSteps_,
                                 this -> seed_,
                                 this -> kb_);
}

System::~System()
{

}

// Move constructor
System::System(System&& other) noexcept
    : lattice_(std::move(other.lattice_)),
      mcs_(std::move(other.mcs_)),
      thermalizationSteps_(std::move(other.thermalizationSteps_)),
      measurementSteps_(std::move(other.measurementSteps_)),
      thermalizationFraction_(std::move(other.thermalizationFraction_)),
      kb_(std::move(other.kb_)),
      seed_(std::move(other.seed_)),
      temps_(std::move(other.temps_)),
      fields_(std::move(other.fields_)),
      outName_(std::move(other.outName_)),
      magnetizationByTypeIndex_(std::move(other.magnetizationByTypeIndex_)),
      engine_(std::move(other.engine_)),
      realRandomGenerator_(std::move(other.realRandomGenerator_)),
      intRandomGenerator_(std::move(other.intRandomGenerator_)),
      gaussianRandomGenerator_(std::move(other.gaussianRandomGenerator_)),
      reporter_(std::move(other.reporter_)),
      sigma_(std::move(other.sigma_)),
      counterRejections_(std::move(other.counterRejections_)),
      num_types_(std::move(other.num_types_))
{
    // Update random generator distributions with new bounds if needed
    intRandomGenerator_ = std::uniform_int_distribution<>(0, lattice_.getAtoms().size() - 1);
}

// Move assignment operator
System& System::operator=(System&& other) noexcept
{
    if (this != &other) {
        lattice_ = std::move(other.lattice_);
        mcs_ = std::move(other.mcs_);
        thermalizationSteps_ = std::move(other.thermalizationSteps_);
        measurementSteps_ = std::move(other.measurementSteps_);
        thermalizationFraction_ = std::move(other.thermalizationFraction_);
        kb_ = std::move(other.kb_);
        seed_ = std::move(other.seed_);
        temps_ = std::move(other.temps_);
        fields_ = std::move(other.fields_);
        outName_ = std::move(other.outName_);
        magnetizationByTypeIndex_ = std::move(other.magnetizationByTypeIndex_);
        engine_ = std::move(other.engine_);
        realRandomGenerator_ = std::move(other.realRandomGenerator_);
        intRandomGenerator_ = std::move(other.intRandomGenerator_);
        gaussianRandomGenerator_ = std::move(other.gaussianRandomGenerator_);
        reporter_ = std::move(other.reporter_);
        sigma_ = std::move(other.sigma_);
        counterRejections_ = std::move(other.counterRejections_);
        num_types_ = std::move(other.num_types_);
        
        // Update random generator distributions with new bounds
        intRandomGenerator_ = std::uniform_int_distribution<>(0, lattice_.getAtoms().size() - 1);
    }
    return *this;
}

void System::ComputeMagnetization()
{
    for (auto& val : this -> magnetizationByTypeIndex_)
        val = 0.0;

    for (auto& atom : this -> lattice_.getAtoms())
    {
        this -> magnetizationByTypeIndex_.at(atom.getTypeIndex()) += atom.getSpin();
        this -> magnetizationByTypeIndex_.at(this -> num_types_) += atom.getSpin();
    }
}

Real System::localEnergy(Index index, Real H)
{
    const Atom& atom = this -> lattice_.getAtoms().at(index);
    return this -> localEnergy(atom, H);
}

Real System::localEnergy(const Atom& atom, Real H)
{
    Real energy = 0.0;
    energy += atom.getExchangeEnergy(this -> lattice_.getAtoms());
    energy += atom.getAnisotropyEnergy(atom);
    energy += atom.getZeemanEnergy(H);
    return energy;
}

Real System::totalEnergy(Real H)
{
    Real exchange_energy = 0.0;
    Real other_energy = 0.0;
    for (auto& atom : this -> lattice_.getAtoms())
    {
        exchange_energy += atom.getExchangeEnergy(this -> lattice_.getAtoms());
        other_energy += atom.getAnisotropyEnergy(atom);
        other_energy += atom.getZeemanEnergy(H);
    }
    // Apply 0.5 correction because bonds are stored bidirectionally
    // (each bond counted twice: once for each direction)
    return 0.5 * exchange_energy + other_energy;
}

void System::randomizeSpins()
{
    for (auto& atom : this -> lattice_.getAtoms())
        atom.initializeRandomState(this -> engine_,
            this -> realRandomGenerator_,
            this -> gaussianRandomGenerator_);
}


void System::monteCarloStep(Real T, Real H)
{
    Index num = Index(this -> realRandomGenerator_(this -> engine_) * NUM_SPIN_MODELS);
    for (Index _ = 0; _ < this -> lattice_.getAtoms().size(); ++_)
    {
        Index randIndex = this -> intRandomGenerator_(this -> engine_);
        Atom& atom = this -> lattice_.getAtoms().at(randIndex);
        Real oldEnergy = this -> localEnergy(atom, H);
        atom.randomizeSpin(this -> engine_,
            this -> realRandomGenerator_,
            this -> gaussianRandomGenerator_,
            this -> sigma_.at(atom.getTypeIndex()), num);
        Real newEnergy = this -> localEnergy(atom, H);
        Real deltaEnergy = newEnergy - oldEnergy;

        if (deltaEnergy > 0 && this -> realRandomGenerator_(this -> engine_) > std::exp(- deltaEnergy / (this -> kb_ * T)))
        {
            atom.revertSpin();
            this -> counterRejections_.at(atom.getTypeIndex()) += 1;
        }
    }
}

void System::cycle()
{

    std::vector< std::vector<Real> > histMag_x(this -> num_types_ + 1);
    std::vector< std::vector<Real> > histMag_y(this -> num_types_ + 1);
    std::vector< std::vector<Real> > histMag_z(this -> num_types_ + 1);

    for (Index i = 0; i <= this -> num_types_; ++i)
    {
        histMag_x.at(i) = std::vector<Real>(0);
        histMag_y.at(i) = std::vector<Real>(0);
        histMag_z.at(i) = std::vector<Real>(0);
    }

    Index initial_time = 0;
    Index final_time = 0;
    Real av_time_per_step = 0.0;

    for (Index index = 0; index < this -> temps_.size(); ++index)
    {
        initial_time = time(NULL);

        Real T = this -> temps_.at(index);
        Real H = this -> fields_.at(index);
        std::vector<Real> enes;

        for (Index i = 0; i <= this -> num_types_; ++i)
        {
            histMag_x.at(i).clear();
            histMag_y.at(i).clear();
            histMag_z.at(i).clear();
        }

        // Reset sigma at the start of each temperature/field point
        this -> resetSigma();

        // Phase 1: Thermalization (adapt sigma, no measurements)
        for (Index _ = 0; _ < this -> thermalizationSteps_; ++_)
        {
            this -> monteCarloStep(T, H);
            this -> adaptSigma();
        }

        // Phase 2: Measurement (fixed sigma, store observables)
        for (Index _ = 0; _ < this -> measurementSteps_; ++_)
        {
            this -> monteCarloStep(T, H);
            enes.push_back(this -> totalEnergy(H));
            this -> ComputeMagnetization();

            Index i = 0;
            for (auto& val : this -> counterRejections_)
            {
                // Store magnetization for each type
                histMag_x.at(i).push_back(this -> magnetizationByTypeIndex_.at(i)[0]);
                histMag_y.at(i).push_back(this -> magnetizationByTypeIndex_.at(i)[1]);
                histMag_z.at(i).push_back(this -> magnetizationByTypeIndex_.at(i)[2]);
                
                // Reset rejection counters (but don't adapt sigma during measurement!)
                this -> counterRejections_.at(i) = 0;
                i++;
            }

            // Store total magnetization
            histMag_x.at(this -> num_types_).push_back(this -> magnetizationByTypeIndex_.at(this -> num_types_)[0]);
            histMag_y.at(this -> num_types_).push_back(this -> magnetizationByTypeIndex_.at(this -> num_types_)[1]);
            histMag_z.at(this -> num_types_).push_back(this -> magnetizationByTypeIndex_.at(this -> num_types_)[2]);

        }
        this -> reporter_.partial_report(enes, histMag_x, histMag_y, histMag_z, this -> lattice_, index);


        final_time = time(NULL);
        av_time_per_step = (av_time_per_step*index + final_time - initial_time) / (index + 1);

        rlutil::saveDefaultColor();
        rlutil::setColor(rlutil::YELLOW);

        std::cout << ETA_seconds(av_time_per_step * (this -> temps_.size() - index));
        rlutil::setColor(rlutil::LIGHTBLUE);
        std::cout << std::setprecision(5) << std::fixed;
        std::cout << "\t("
                  << 100.0 * (index + 1) / (this -> temps_.size())
                  << "%)";
        rlutil::resetColor();
        // std::cout << "\t==>\tT = " << T << "; H = " << H << std::endl;
        std::cout << "\t==>\tT = " << T << "; H = " << H << " ";
        Index i = 0;
        for (auto& element : this -> lattice_.getMapTypeIndexes())
        {
            std::cout << element.first << " " << this -> sigma_.at(i) << " ";
            i++;
        }
        std::cout << std::endl;

    }

    this -> reporter_.close();

}

void System::resetSigma()
{
    for (Index i = 0; i < this -> sigma_.size(); ++i)
    {
        this -> sigma_.at(i) = MAX_SIGMA;
        this -> counterRejections_.at(i) = 0;
    }
}

void System::adaptSigma()
{
    Real rejection;
    Real sigma_temp;
    
    Index i = 0;
    for (auto& val : this -> counterRejections_)
    {
        rejection = val / Real(this -> lattice_.getSizesByIndex().at(i));
        // Avoid division by zero or very small rejection rates
        if (rejection < MIN_REJECTION_RATE)
        {
            sigma_temp = MAX_SIGMA; // Use maximum sigma when rejection rate is negligible
        }
        else
        {
            sigma_temp = this -> sigma_.at(i) * (SIGMA_ADJUSTMENT_FACTOR / rejection);
            sigma_temp = std::max(MIN_SIGMA, std::min(sigma_temp, MAX_SIGMA));
        }
        this -> sigma_.at(i) = sigma_temp;
        this -> counterRejections_.at(i) = 0;
        i++;
    }
}

Lattice& System::getLattice()
{
    return this -> lattice_;
}

Index System::getSeed() const
{
    return this -> seed_;
}

void System::setState(std::string fileState)
{
    std::ifstream file(fileState);

    Array spin;
    std::string sx;
    std::string sy;
    std::string sz;
    for (auto& atom : this -> lattice_.getAtoms())
    {
        file >> sx >> sy >> sz;
        Array spin({safe_stod(sx, "spin x-component"), 
                     safe_stod(sy, "spin y-component"), 
                     safe_stod(sz, "spin z-component")});
        Real norm = std::round(std::sqrt((spin*spin).sum()) * SPIN_NORM_ROUNDING_PRECISION) / SPIN_NORM_ROUNDING_PRECISION;
        if (fp_not_equal(norm, atom.getSpinNorm()))
        {
            std::cout << norm << " " << atom.getSpinNorm() << " "
                      << fp_equal(norm, atom.getSpinNorm()) << " " << spin
                      << std::endl;
            throw vegas::SimulationException("The spin norm of the site " + std::to_string(atom.getIndex()) + " does not match with the initial state given !!!");
        }

        atom.setSpin(spin);
    }
}

void System::setAnisotropies(std::vector<std::string> anisotropyfiles)
{
    for (auto& fileName : anisotropyfiles)
    {
        std::ifstream file(fileName);
        for (Index i = 0; i < this -> lattice_.getAtoms().size(); ++i)
        {
            std::string line;
            std::getline(file, line);
            std::vector<std::string> sep = split(line, ' ');

            if (sep.size() == 4) // add an uniaxial term
            {
                Real ax = safe_stod(sep[0], "anisotropy x-component");
                Real ay = safe_stod(sep[1], "anisotropy y-component");
                Real az = safe_stod(sep[2], "anisotropy z-component");
                Real kan = safe_stod(sep[3], "anisotropy constant");

                // Normalize anisotropy unit vector
                Real norm = std::sqrt(ax*ax + ay*ay + az*az);
                if (norm > 0.0) {
                    ax /= norm;
                    ay /= norm;
                    az /= norm;
                }

                // Store anisotropy unit vector and constant in atom
                this -> lattice_.getAtoms().at(i).setAnisotropyUnit(Array{ax, ay, az});
                this -> lattice_.getAtoms().at(i).setKan(kan);
                this -> lattice_.getAtoms().at(i).setTypeAnisotropy("uniaxial");

                std::function<Real(const Atom&)> func = [kan, ax, ay, az](const Atom& atom){
                   return - kan * (ax * atom.getSpin()[0] + ay * atom.getSpin()[1] + az * atom.getSpin()[2]) * (ax * atom.getSpin()[0] + ay * atom.getSpin()[1] + az * atom.getSpin()[2]);
                };
                this -> lattice_.getAtoms().at(i).addAnisotropyTerm(func);
            }
            else if (sep.size() == 7)
            {
                Real Ax = safe_stod(sep[0], "anisotropy A x-component");
                Real Ay = safe_stod(sep[1], "anisotropy A y-component");
                Real Az = safe_stod(sep[2], "anisotropy A z-component");
                Real Bx = safe_stod(sep[3], "anisotropy B x-component");
                Real By = safe_stod(sep[4], "anisotropy B y-component");
                Real Bz = safe_stod(sep[5], "anisotropy B z-component");

                Real Cx = Ay*Bz - Az*By;
                Real Cy = Az*Bx - Ax*Bz;
                Real Cz = Ax*By - Ay*Bx;

                Array A = {Ax, Ay, Az};
                Array B = {Bx, By, Bz};
                Array C = {Cx, Cy, Cz};

                Real kan = safe_stod(sep[6], "anisotropy constant");

                // Store anisotropy unit vector (use A direction) and constant in atom
                // Normalize A for storage
                Real normA = std::sqrt(Ax*Ax + Ay*Ay + Az*Az);
                if (normA > 0.0) {
                    this -> lattice_.getAtoms().at(i).setAnisotropyUnit(Array{Ax/normA, Ay/normA, Az/normA});
                }
                this -> lattice_.getAtoms().at(i).setKan(kan);
                this -> lattice_.getAtoms().at(i).setTypeAnisotropy("biaxial");

                std::function<Real(const Atom&)> func = [kan, A, B, C](const Atom& atom){
                    return - kan * ((atom.getSpin() * A).sum()*(atom.getSpin() * A).sum()*(atom.getSpin() * B).sum()*(atom.getSpin() * B).sum()
                    + (atom.getSpin() * A).sum()*(atom.getSpin() * A).sum()*(atom.getSpin() * C).sum()*(atom.getSpin() * C).sum()
                    + (atom.getSpin() * B).sum()*(atom.getSpin() * B).sum()*(atom.getSpin() * C).sum()*(atom.getSpin() * C).sum());
                };

                this -> lattice_.getAtoms().at(i).addAnisotropyTerm(func);

            }
            else
            {
                throw vegas::InvalidInputException("The anisotropy file with name " + fileName + " does not have the correct format !!!");
            }
        }
    }
}
