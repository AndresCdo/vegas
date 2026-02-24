#ifndef SYSTEM
#define SYSTEM

#include <random>
#include <cmath>
#include "lattice.h"
#include "reporter.h"
#include "random.h"
#include "simd_rng.h"
#include "spin_model_tags.h"

class System
{
public:
    System(std::string fileName,
           std::vector<Real> temps,
           std::vector<Real> fields,
           Index mcs,
           Index seed,
           std::string outName,
           Real kb);
    ~System();
    
    // Move constructor and assignment
    System(System&& other) noexcept;
    System& operator=(System&& other) noexcept;
    
    // Delete copy constructor and assignment
    System(const System&) = delete;
    System& operator=(const System&) = delete;

    void ComputeMagnetization();
    Real localEnergy(Index index, Real H);
    Real localEnergy(const Atom& atom, Real H);
    Real localEnergySoA(Index index, Real H);  // SoA-native version (Heisenberg)


    Real totalEnergy(Real H);
    
    void randomizeSpins();
    
    void monteCarloStep(Real T, Real H);

    // Template implementations for dispatch
    template<typename ModelTag>
    void monteCarloStepImpl(Real T, Real H);

    void monteCarloStep_ising(Real T, Real H);
    
    void cycle();
    
    void adaptSigma();
    
    void resetSigma();

    Lattice& getLattice();

    Index getSeed() const;

    const std::map<std::string, Array>& getMagnetizationType() const;

    void setState(std::string fileState);
    
    void writeCheckpoint(std::string filename, Index tempIndex);
    bool loadCheckpoint(std::string filename, Index& tempIndex);
    
    void setAnisotropies(std::vector<std::string> anisotropyfiles);

private:
    Lattice lattice_;
    Index mcs_;
    Index thermalizationSteps_;
    Index measurementSteps_;
    Real thermalizationFraction_;
    Real kb_;
    Index seed_;
    std::vector<Real> temps_;
    std::vector<Real> fields_;
    std::string outName_;

    std::vector<Array> magnetizationByTypeIndex_;

    vegas::Xoshiro256StarStar engine_;
    vegas::SimdRNG<SIMD_LANES> simdEngine_;
    std::uniform_real_distribution<> realRandomGenerator_;
    std::uniform_int_distribution<> intRandomGenerator_;
    std::normal_distribution<> gaussianRandomGenerator_;

    Reporter reporter_;

    std::vector<Real> sigma_;
    std::vector<Index> counterRejections_;

    Index num_types_;
    bool isIsing_;  // Model type for template dispatch
    Index resumeIndex_;  // Temperature index to resume from (0 = start fresh)
};

#endif
