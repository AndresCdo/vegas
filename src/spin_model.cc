#include "../include/spin_model.h"
#include "../include/atom.h"
#include "../include/params.h"

#include <cmath>
#include <random>

// HeisenbergSpinModel implementation
void HeisenbergSpinModel::initializeRandomState(
    std::mt19937_64& engine,
    std::uniform_real_distribution<>& realRandomGenerator,
    std::normal_distribution<>& gaussianRandomGenerator,
    Atom& atom) const
{
    atom.setOldSpin(atom.getSpin());
    Array gamma({gaussianRandomGenerator(engine), gaussianRandomGenerator(engine), gaussianRandomGenerator(engine)});
    Array unitArray = gamma / std::sqrt((gamma * gamma).sum());
    atom.setSpin(atom.getSpinNorm() * unitArray);
}

void HeisenbergSpinModel::randomizeSpin(
    std::mt19937_64& engine,
    std::uniform_real_distribution<>& realRandomGenerator,
    std::normal_distribution<>& gaussianRandomGenerator,
    Real sigma,
    Atom& atom,
    Index num) const
{
    atom.setOldSpin(atom.getSpin());
    Array gamma({gaussianRandomGenerator(engine), gaussianRandomGenerator(engine), gaussianRandomGenerator(engine)});
    Array unitArray = gamma / std::sqrt((gamma * gamma).sum());
    atom.setSpin(atom.getSpinNorm() * unitArray);
}

// IsingSpinModel implementation
void IsingSpinModel::initializeRandomState(
    std::mt19937_64& engine,
    std::uniform_real_distribution<>& realRandomGenerator,
    std::normal_distribution<>& gaussianRandomGenerator,
    Atom& atom) const
{
    atom.setOldSpin(atom.getSpin());
    if (realRandomGenerator(engine) < 0.5)
        atom.setSpin({0.0, 0.0, -atom.getSpinNorm()});
    else
        atom.setSpin({0.0, 0.0, atom.getSpinNorm()});
}

void IsingSpinModel::randomizeSpin(
    std::mt19937_64& engine,
    std::uniform_real_distribution<>& realRandomGenerator,
    std::normal_distribution<>& gaussianRandomGenerator,
    Real sigma,
    Atom& atom,
    Index num) const
{
    atom.setOldSpin(atom.getSpin());
    atom.setSpin(-atom.getSpin());
}

// Quantum Ising model (qising)
class QuantumIsingSpinModel : public SpinModel {
public:
    void initializeRandomState(
        std::mt19937_64& engine,
        std::uniform_real_distribution<>& realRandomGenerator,
        std::normal_distribution<>& gaussianRandomGenerator,
        Atom& atom) const override
    {
        atom.setOldSpin(atom.getSpin());
        atom.setSproj(int(realRandomGenerator(engine) * atom.getPossibleProjections().size()));
        atom.setSpin({0.0, 0.0, atom.getPossibleProjections()[atom.getSproj()]});
        atom.changeProjection(atom.getSproj(), atom.getOldSpin()[2]);
    }
    
    void randomizeSpin(
        std::mt19937_64& engine,
        std::uniform_real_distribution<>& realRandomGenerator,
        std::normal_distribution<>& gaussianRandomGenerator,
        Real sigma,
        Atom& atom,
        Index num) const override
    {
        atom.setOldSpin(atom.getSpin());
        atom.setSproj(int(realRandomGenerator(engine) * atom.getPossibleProjections().size()));
        atom.setSpin({0.0, 0.0, atom.getPossibleProjections()[atom.getSproj()]});
        atom.changeProjection(atom.getSproj(), atom.getOldSpin()[2]);
    }
    
    std::string getName() const override { return "qising"; }
};

// Adaptive model (combines random and flip)
class AdaptiveSpinModel : public SpinModel {
public:
    void initializeRandomState(
        std::mt19937_64& engine,
        std::uniform_real_distribution<>& realRandomGenerator,
        std::normal_distribution<>& gaussianRandomGenerator,
        Atom& atom) const override
    {
        atom.setOldSpin(atom.getSpin());
        Array gamma({gaussianRandomGenerator(engine), gaussianRandomGenerator(engine), gaussianRandomGenerator(engine)});
        Array unitArray = gamma / std::sqrt((gamma * gamma).sum());
        atom.setSpin(atom.getSpinNorm() * unitArray);
    }
    
    void randomizeSpin(
        std::mt19937_64& engine,
        std::uniform_real_distribution<>& realRandomGenerator,
        std::normal_distribution<>& gaussianRandomGenerator,
        Real sigma,
        Atom& atom,
        Index num) const override
    {
        atom.setOldSpin(atom.getSpin());
        
        // Adaptive logic: combine random and flip based on sigma
        if (realRandomGenerator(engine) < sigma) {
            // Random move
            Array gamma({gaussianRandomGenerator(engine), gaussianRandomGenerator(engine), gaussianRandomGenerator(engine)});
            Array unitArray = gamma / std::sqrt((gamma * gamma).sum());
            atom.setSpin(atom.getSpinNorm() * unitArray);
        } else {
            // Flip move
            atom.setSpin(-atom.getSpin());
        }
    }
    
    std::string getName() const override { return "adaptive"; }
};

// Cone models (cone30, cone15)
class ConeSpinModel : public SpinModel {
private:
    Real coneAngle_;
    
public:
    ConeSpinModel(Real coneAngle) : coneAngle_(coneAngle) {}
    
    void initializeRandomState(
        std::mt19937_64& engine,
        std::uniform_real_distribution<>& realRandomGenerator,
        std::normal_distribution<>& gaussianRandomGenerator,
        Atom& atom) const override
    {
        atom.setOldSpin(atom.getSpin());
        Array gamma({gaussianRandomGenerator(engine), gaussianRandomGenerator(engine), gaussianRandomGenerator(engine)});
        Array unitArray = gamma / std::sqrt((gamma * gamma).sum());
        atom.setSpin(atom.getSpinNorm() * unitArray);
    }
    
    void randomizeSpin(
        std::mt19937_64& engine,
        std::uniform_real_distribution<>& realRandomGenerator,
        std::normal_distribution<>& gaussianRandomGenerator,
        Real sigma,
        Atom& atom,
        Index num) const override
    {
        atom.setOldSpin(atom.getSpin());
        
        // Cone move: random rotation within a cone
        Real theta = coneAngle_ * realRandomGenerator(engine);
        Real phi = 2.0 * M_PI * realRandomGenerator(engine);
        
        Array newSpin = atom.getSpinNorm() * Array{
            std::sin(theta) * std::cos(phi),
            std::sin(theta) * std::sin(phi),
            std::cos(theta)
        };
        
        atom.setSpin(newSpin);
    }
    
    std::string getName() const override { 
        if (fp_equal(coneAngle_, 30.0 * M_PI / 180.0)) {
            return "cone30";
        } else {
            return "cone15";
        }
    }
};

// HN models (hn30, hn15) - Hybrid Nematic models
class HNSpinModel : public SpinModel {
private:
    Real coneAngle_;
    
public:
    HNSpinModel(Real coneAngle) : coneAngle_(coneAngle) {}
    
    void initializeRandomState(
        std::mt19937_64& engine,
        std::uniform_real_distribution<>& realRandomGenerator,
        std::normal_distribution<>& gaussianRandomGenerator,
        Atom& atom) const override
    {
        atom.setOldSpin(atom.getSpin());
        Array gamma({gaussianRandomGenerator(engine), gaussianRandomGenerator(engine), gaussianRandomGenerator(engine)});
        Array unitArray = gamma / std::sqrt((gamma * gamma).sum());
        atom.setSpin(atom.getSpinNorm() * unitArray);
    }
    
    void randomizeSpin(
        std::mt19937_64& engine,
        std::uniform_real_distribution<>& realRandomGenerator,
        std::normal_distribution<>& gaussianRandomGenerator,
        Real sigma,
        Atom& atom,
        Index num) const override
    {
        atom.setOldSpin(atom.getSpin());
        
        // HN move: combination of cone and flip
        if (realRandomGenerator(engine) < 0.5) {
            // Cone move
            Real theta = coneAngle_ * realRandomGenerator(engine);
            Real phi = 2.0 * M_PI * realRandomGenerator(engine);
            
            Array newSpin = atom.getSpinNorm() * Array{
                std::sin(theta) * std::cos(phi),
                std::sin(theta) * std::sin(phi),
                std::cos(theta)
            };
            
            atom.setSpin(newSpin);
        } else {
            // Flip move
            atom.setSpin(-atom.getSpin());
        }
    }
    
    std::string getName() const override { 
        if (fp_equal(coneAngle_, 30.0 * M_PI / 180.0)) {
            return "hn30";
        } else {
            return "hn15";
        }
    }
};

// Factory function implementation
std::unique_ptr<SpinModel> createSpinModel(const std::string& modelName) {
    if (modelName == "heisenberg" || modelName == "random") {
        return std::make_unique<HeisenbergSpinModel>();
    } else if (modelName == "ising" || modelName == "flip") {
        return std::make_unique<IsingSpinModel>();
    } else if (modelName == "qising") {
        return std::make_unique<QuantumIsingSpinModel>();
    } else if (modelName == "adaptive") {
        return std::make_unique<AdaptiveSpinModel>();
    } else if (modelName == "cone30") {
        return std::make_unique<ConeSpinModel>(30.0 * M_PI / 180.0);
    } else if (modelName == "cone15") {
        return std::make_unique<ConeSpinModel>(15.0 * M_PI / 180.0);
    } else if (modelName == "hn30") {
        return std::make_unique<HNSpinModel>(30.0 * M_PI / 180.0);
    } else if (modelName == "hn15") {
        return std::make_unique<HNSpinModel>(15.0 * M_PI / 180.0);
    } else {
        // Default to Heisenberg model
        return std::make_unique<HeisenbergSpinModel>();
    }
}