#ifndef SPIN_MODEL_H
#define SPIN_MODEL_H

#include "params.h"
#include "random.h"
#include <functional>
#include <memory>
#include <random>
#include <string>

// Forward declaration
class Atom;

// Base class for spin models
class SpinModel {
public:
    virtual ~SpinModel() = default;
    
    // Initialize random spin state
    virtual void initializeRandomState(
        vegas::Xoshiro256StarStar& engine,
        std::uniform_real_distribution<>& realRandomGenerator,
        std::normal_distribution<>& gaussianRandomGenerator,
        Atom& atom) const = 0;
    
    // Randomize spin
    virtual void randomizeSpin(
        vegas::Xoshiro256StarStar& engine,
        std::uniform_real_distribution<>& realRandomGenerator,
        std::normal_distribution<>& gaussianRandomGenerator,
        Real sigma,
        Atom& atom,
        Index num) const = 0;
    
    // Get model name
    virtual std::string getName() const = 0;
};

// Concrete spin models
class HeisenbergSpinModel : public SpinModel {
public:
    void initializeRandomState(
        vegas::Xoshiro256StarStar& engine,
        std::uniform_real_distribution<>& realRandomGenerator,
        std::normal_distribution<>& gaussianRandomGenerator,
        Atom& atom) const override;
    
    void randomizeSpin(
        vegas::Xoshiro256StarStar& engine,
        std::uniform_real_distribution<>& realRandomGenerator,
        std::normal_distribution<>& gaussianRandomGenerator,
        Real sigma,
        Atom& atom,
        Index num) const override;
    
    std::string getName() const override { return "heisenberg"; }
};

class IsingSpinModel : public SpinModel {
public:
    void initializeRandomState(
        vegas::Xoshiro256StarStar& engine,
        std::uniform_real_distribution<>& realRandomGenerator,
        std::normal_distribution<>& gaussianRandomGenerator,
        Atom& atom) const override;
    
    void randomizeSpin(
        vegas::Xoshiro256StarStar& engine,
        std::uniform_real_distribution<>& realRandomGenerator,
        std::normal_distribution<>& gaussianRandomGenerator,
        Real sigma,
        Atom& atom,
        Index num) const override;
    
    std::string getName() const override { return "ising"; }
};

// Factory function to create spin models
std::unique_ptr<SpinModel> createSpinModel(const std::string& modelName);

#endif // SPIN_MODEL_H