#include "../include/spin_model.h"
#include "../include/atom.h"
#include "../include/params.h"
#include "../include/random.h"

#include <cmath>
#include <random>

namespace {
    Array cross(const Array& A, const Array& B)
    {
        return Array{A[1]*B[2] - A[2]*B[1], A[2]*B[0] - A[0]*B[2], A[0]*B[1] - A[1]*B[0]};
    }
    
    Real dot(const Array& A, const Array& B)
    {
        return (A * B).sum();
    }
    
    Real norm(const Array& A)
    {
        return std::sqrt((A * A).sum());
    }
    
    Array rotateVectorToAlignWith(const Array& v, const Array& targetDirection)
    {
        Array zAxis = {0.0, 0.0, 1.0};
        Real targetNorm = norm(targetDirection);
        
        if (targetNorm < EPSILON) {
            return v;
        }
        
        Array normalizedTarget = targetDirection / targetNorm;
        Real cosAngle = dot(zAxis, normalizedTarget);
        
        if (std::abs(cosAngle - 1.0) < EPSILON) {
            return v;
        }
        
        if (std::abs(cosAngle + 1.0) < EPSILON) {
            return -v;
        }
        
        Array rotationAxis = cross(zAxis, normalizedTarget);
        Real axisNorm = norm(rotationAxis);
        
        if (axisNorm < EPSILON) {
            return v;
        }
        
        rotationAxis = rotationAxis / axisNorm;
        Real sinAngle = std::sqrt(1.0 - cosAngle * cosAngle);
        
        Array rotated = v * cosAngle + 
                        cross(rotationAxis, v) * sinAngle + 
                        rotationAxis * dot(rotationAxis, v) * (1.0 - cosAngle);
        
        return rotated;
    }
}

// HeisenbergSpinModel implementation
void HeisenbergSpinModel::initializeRandomState(
    vegas::Xoshiro256StarStar& engine,
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
    vegas::Xoshiro256StarStar& engine,
    std::uniform_real_distribution<>& realRandomGenerator,
    std::normal_distribution<>& gaussianRandomGenerator,
    Real sigma,
    Atom& atom,
    Index num) const
{
    atom.setOldSpin(atom.getSpin());
    Array currentSpin = atom.getSpin();
    Real spinNorm = atom.getSpinNorm();
    
    Array perturbation({
        sigma * gaussianRandomGenerator(engine),
        sigma * gaussianRandomGenerator(engine),
        sigma * gaussianRandomGenerator(engine)
    });
    
    Array newSpin = currentSpin + perturbation;
    Real newNorm = norm(newSpin);
    
    if (newNorm > EPSILON) {
        atom.setSpin(spinNorm * newSpin / newNorm);
    } else {
        Array gamma({gaussianRandomGenerator(engine), gaussianRandomGenerator(engine), gaussianRandomGenerator(engine)});
        Array unitArray = gamma / std::sqrt((gamma * gamma).sum());
        atom.setSpin(spinNorm * unitArray);
    }
}

// IsingSpinModel implementation
void IsingSpinModel::initializeRandomState(
    vegas::Xoshiro256StarStar& engine,
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
    vegas::Xoshiro256StarStar& engine,
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
        vegas::Xoshiro256StarStar& engine,
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
        vegas::Xoshiro256StarStar& engine,
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
        vegas::Xoshiro256StarStar& engine,
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
        vegas::Xoshiro256StarStar& engine,
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
        vegas::Xoshiro256StarStar& engine,
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
        vegas::Xoshiro256StarStar& engine,
        std::uniform_real_distribution<>& realRandomGenerator,
        std::normal_distribution<>& gaussianRandomGenerator,
        Real sigma,
        Atom& atom,
        Index num) const override
    {
        atom.setOldSpin(atom.getSpin());
        Array currentSpin = atom.getSpin();
        Real spinNorm = atom.getSpinNorm();
        
        Real theta = coneAngle_ * realRandomGenerator(engine);
        Real phi = 2.0 * M_PI * realRandomGenerator(engine);
        
        Array localPerturbation = Array{
            std::sin(theta) * std::cos(phi),
            std::sin(theta) * std::sin(phi),
            std::cos(theta)
        };
        
        Array rotatedPerturbation = rotateVectorToAlignWith(localPerturbation, currentSpin);
        
        atom.setSpin(spinNorm * rotatedPerturbation);
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
        vegas::Xoshiro256StarStar& engine,
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
        vegas::Xoshiro256StarStar& engine,
        std::uniform_real_distribution<>& realRandomGenerator,
        std::normal_distribution<>& gaussianRandomGenerator,
        Real sigma,
        Atom& atom,
        Index num) const override
    {
        atom.setOldSpin(atom.getSpin());
        Array currentSpin = atom.getSpin();
        Real spinNorm = atom.getSpinNorm();
        
        if (realRandomGenerator(engine) < 0.5) {
            Real theta = coneAngle_ * realRandomGenerator(engine);
            Real phi = 2.0 * M_PI * realRandomGenerator(engine);
            
            Array localPerturbation = Array{
                std::sin(theta) * std::cos(phi),
                std::sin(theta) * std::sin(phi),
                std::cos(theta)
            };
            
            Array rotatedPerturbation = rotateVectorToAlignWith(localPerturbation, currentSpin);
            atom.setSpin(spinNorm * rotatedPerturbation);
        } else {
            atom.setSpin(-currentSpin);
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