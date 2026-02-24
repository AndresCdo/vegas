#include "../include/system.h"
#include "../include/random.h"
#include "../include/spin_model_tags.h"
#include "../include/starter.h"
#include <hdf5.h>
#if defined(__AVX2__) && defined(VEGAS_AVX2_EXPERIMENTAL)
#include <immintrin.h>
#endif
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
#include <filesystem>

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
               Real kb) : lattice_(fileName),
                           engine_(0), simdEngine_(0)
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
    
    // Detect if system is Ising or Heisenberg for template dispatch
    const auto& atoms = this -> lattice_.getAtoms();
    if (!atoms.empty()) {
        std::string model = atoms[0].getModel();
        this -> isIsing_ = (model == "ising");
    } else {
        this -> isIsing_ = false;
    }
    
    this -> resumeIndex_ = 0;  // Default: start from beginning

    std::random_device rd;
    uint64_t base_seed = (seed != 0) ? static_cast<uint64_t>(seed) : rd();
    this -> engine_ = vegas::Xoshiro256StarStar(base_seed);
    this -> simdEngine_ = vegas::SimdRNG<SIMD_LANES>(base_seed);
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

    // Validate output path to prevent path traversal
    if (outName.find("..") != std::string::npos) {
        throw vegas::InvalidInputException("Output filename contains path traversal attempt: " + outName);
    }
    std::remove(outName.c_str());

    // Initialize SoA from Atoms after lattice is loaded
    this -> lattice_.syncAtomsToSoA();

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
       simdEngine_(std::move(other.simdEngine_)),
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
        simdEngine_ = std::move(other.simdEngine_);
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

// SoA-native local energy calculation (Heisenberg model)
// Uses direct array access instead of Atom objects
Real System::localEnergySoA(Index i, Real H)
{
    // Get SoA data pointers
    const Real* spin_x = lattice_.getSpinX();
    const Real* spin_y = lattice_.getSpinY();
    const Real* spin_z = lattice_.getSpinZ();
    
    const auto& neighborIndexes = lattice_.getNeighborIndexes()[i];
    const auto& exchanges = lattice_.getExchanges()[i];
    
    // Current spin components
    Real sx = spin_x[i];
    Real sy = spin_y[i];
    Real sz = spin_z[i];
    
    // Exchange energy
    Real exchange_energy = 0.0;
    for (size_t j = 0; j < neighborIndexes.size(); ++j) {
        Index nbh = neighborIndexes[j];
        exchange_energy -= exchanges[j] * (
            sx * spin_x[nbh] + 
            sy * spin_y[nbh] + 
            sz * spin_z[nbh]
        );
    }
    
    // Zeeman energy (assuming external field in z-direction for simplicity)
    // For full support, would need external field from lattice
    Real zeeman_energy = 0.0;
    
    // Apply 0.5 correction for bidirectional bonds
    return 0.5 * exchange_energy + zeeman_energy;
}

void System::randomizeSpins()
{
    for (auto& atom : this -> lattice_.getAtoms())
        atom.initializeRandomState(this -> engine_,
            this -> realRandomGenerator_,
            this -> gaussianRandomGenerator_);
    // Sync updated spins to SoA
    this -> lattice_.syncAtomsToSoA();
}

// Template dispatch: Heisenberg implementation
template<>
void System::monteCarloStepImpl<HeisenbergTag>(Real T, Real H)
{
    Index N = this -> lattice_.getNumAtoms();
    Real* __restrict spin_xyz = lattice_.getSpinXYZ();
    
    const std::vector<Index>& typeIndices = lattice_.getTypeIndices();
    const std::vector<Real>& spinNorms = lattice_.getSpinNorms();
    const std::vector<Index>& soaNeighborOffsets = lattice_.getSoANeighborOffsets();
    const Index* flatNeighborIds = lattice_.getFlatNeighborIds();
    const Real* flatNeighborJs = lattice_.getFlatNeighborJs();
    const std::vector<ColorBlock>& colorBlocks = lattice_.getColorBlocks();
    Index numColors = lattice_.getNumColors();
    const std::vector<Index>& globalToSoA = lattice_.getGlobalToSoA();
    
    // Random permutation of colors to avoid bias (optional)
    std::vector<Index> colorOrder(numColors);
    for (Index c = 0; c < numColors; ++c) colorOrder[c] = c;
    std::shuffle(colorOrder.begin(), colorOrder.end(), engine_);
    
    for (Index colorIdx : colorOrder) {
        const ColorBlock& block = colorBlocks[colorIdx];
        Index blockSize = block.getSize();
        const std::vector<Index>& blockAtomIndices = block.getAtomIndices();  // global indices
        
        // Process atoms in groups of SIMD_LANES
        auto& simdEngine = this->simdEngine_;
        for (Index blockPos = 0; blockPos < blockSize; blockPos += SIMD_LANES) {
            Index laneCount = std::min(static_cast<Index>(SIMD_LANES), blockSize - blockPos);
            
            // Generate random numbers for each lane
            auto gaussian1 = simdEngine.gaussian(); // array<Real, SIMD_LANES>
            auto gaussian2 = simdEngine.gaussian();
            auto gaussian3 = simdEngine.gaussian();
            auto uniform = simdEngine.uniform();
            
            // Process each lane
            for (Index lane = 0; lane < laneCount; ++lane) {
                Index globalIdx = blockAtomIndices[blockPos + lane];
                Index soaIdx = globalToSoA[globalIdx];
                
                Index typeIdx = typeIndices[soaIdx];
                Real sigma = this -> sigma_.at(typeIdx);
                
                Real old_sx = spin_xyz[3*soaIdx];
                Real old_sy = spin_xyz[3*soaIdx + 1];
                Real old_sz = spin_xyz[3*soaIdx + 2];
                
                Index neighStart = soaNeighborOffsets[soaIdx];
                Index neighEnd = soaNeighborOffsets[soaIdx + 1];
                
                // Heisenberg: full 3D dot product
                Real oldEnergy = 0.0;
                for (Index j = neighStart; j < neighEnd; ++j) {
                    Index nbh = flatNeighborIds[j];
                    Real J = flatNeighborJs[j];
                    oldEnergy -= J * (
                        old_sx * spin_xyz[3*nbh] + 
                        old_sy * spin_xyz[3*nbh + 1] + 
                        old_sz * spin_xyz[3*nbh + 2]
                    );
                }
                oldEnergy *= 0.5;
                
                // Heisenberg: 3D Gaussian perturbation
                Real prop_sx = old_sx + sigma * gaussian1[lane];
                Real prop_sy = old_sy + sigma * gaussian2[lane];
                Real prop_sz = old_sz + sigma * gaussian3[lane];
                
                Real prop_norm = std::sqrt(prop_sx * prop_sx + prop_sy * prop_sy + prop_sz * prop_sz);
                if (prop_norm > EPSILON) {
                    Real snorm = spinNorms[soaIdx];
                    prop_sx = snorm * prop_sx / prop_norm;
                    prop_sy = snorm * prop_sy / prop_norm;
                    prop_sz = snorm * prop_sz / prop_norm;
                }
                
                spin_xyz[3*soaIdx] = prop_sx;
                spin_xyz[3*soaIdx + 1] = prop_sy;
                spin_xyz[3*soaIdx + 2] = prop_sz;
                
                // Heisenberg energy
                Real newEnergy = 0.0;
                for (Index j = neighStart; j < neighEnd; ++j) {
                    Index nbh = flatNeighborIds[j];
                    Real J = flatNeighborJs[j];
                    newEnergy -= J * (
                         prop_sx * spin_xyz[3*nbh] + 
                         prop_sy * spin_xyz[3*nbh + 1] + 
                         prop_sz * spin_xyz[3*nbh + 2]
                    );
                }
                newEnergy *= 0.5;
                
                Real deltaEnergy = newEnergy - oldEnergy;
                
                if (deltaEnergy > 0 && uniform[lane] > std::exp(- deltaEnergy / (this -> kb_ * T)))
                {
                    spin_xyz[3*soaIdx] = old_sx;
                    spin_xyz[3*soaIdx + 1] = old_sy;
                    spin_xyz[3*soaIdx + 2] = old_sz;
                    this -> counterRejections_.at(typeIdx) += 1;
                }
            }
        }
    }
}

// Template dispatch: Ising implementation
template<>
void System::monteCarloStepImpl<IsingTag>(Real T, Real H)
{
    Index N = this -> lattice_.getNumAtoms();
    Real* __restrict spin_z = lattice_.getSpinZ();
    
    const std::vector<Index>& typeIndices = lattice_.getTypeIndices();
    const std::vector<Real>& spinNorms = lattice_.getSpinNorms();
    const std::vector<Index>& soaNeighborOffsets = lattice_.getSoANeighborOffsets();
    const Index* flatNeighborIds = lattice_.getFlatNeighborIds();
    const Real* flatNeighborJs = lattice_.getFlatNeighborJs();
    const std::vector<ColorBlock>& colorBlocks = lattice_.getColorBlocks();
    Index numColors = lattice_.getNumColors();
    const std::vector<Index>& globalToSoA = lattice_.getGlobalToSoA();
    
    // Random permutation of colors to avoid bias (optional)
    std::vector<Index> colorOrder(numColors);
    for (Index c = 0; c < numColors; ++c) colorOrder[c] = c;
    std::shuffle(colorOrder.begin(), colorOrder.end(), engine_);
    
    for (Index colorIdx : colorOrder) {
        const ColorBlock& block = colorBlocks[colorIdx];
        Index blockSize = block.getSize();
        const std::vector<Index>& blockAtomIndices = block.getAtomIndices();  // global indices
        
        // Process atoms in groups of SIMD_LANES
        auto& simdEngine = this->simdEngine_;
        for (Index blockPos = 0; blockPos < blockSize; blockPos += SIMD_LANES) {
            Index laneCount = std::min(static_cast<Index>(SIMD_LANES), blockSize - blockPos);
            
            // Generate random numbers for each lane
            auto uniform = simdEngine.uniform();
            
#if defined(__AVX2__) && defined(VEGAS_AVX2_EXPERIMENTAL)
            // AVX2 experimental implementation for full vector (4 lanes)
            if (laneCount == SIMD_LANES) {
                // Load global indices for 4 atoms
                Index globalIdx0 = blockAtomIndices[blockPos];
                Index globalIdx1 = blockAtomIndices[blockPos + 1];
                Index globalIdx2 = blockAtomIndices[blockPos + 2];
                Index globalIdx3 = blockAtomIndices[blockPos + 3];
                
                // Convert to SoA indices
                Index soaIdx0 = globalToSoA[globalIdx0];
                Index soaIdx1 = globalToSoA[globalIdx1];
                Index soaIdx2 = globalToSoA[globalIdx2];
                Index soaIdx3 = globalToSoA[globalIdx3];
                
                // Load spin norms
                __m256d snorm_vec = _mm256_set_pd(
                    spinNorms[soaIdx3], spinNorms[soaIdx2], 
                    spinNorms[soaIdx1], spinNorms[soaIdx0]
                );
                
                // Load current spins
                __m256d old_sz_vec = _mm256_set_pd(
                    spin_z[soaIdx3], spin_z[soaIdx2],
                    spin_z[soaIdx1], spin_z[soaIdx0]
                );
                
                // Flip direction: prop_sz = -old_sz
                __m256d prop_sz_vec = _mm256_mul_pd(old_sz_vec, _mm256_set1_pd(-1.0));
                
                // Get neighbor offsets for the 4 atoms
                Index neighStart0 = soaNeighborOffsets[soaIdx0];
                Index neighStart1 = soaNeighborOffsets[soaIdx1];
                Index neighStart2 = soaNeighborOffsets[soaIdx2];
                Index neighStart3 = soaNeighborOffsets[soaIdx3];
                
                // Assume uniform coordination (check first atom's count)
                Index neighCount = soaNeighborOffsets[soaIdx0 + 1] - neighStart0;
                
                // Vectorized local field accumulation
                __m256d h_local_vec = _mm256_setzero_pd();
                
                // Process each neighbor position
                for (Index k = 0; k < neighCount; ++k) {
                    // Gather neighbor indices for the 4 atoms
                    __m256i neighbor_idx_vec = _mm256_set_epi64x(
                        flatNeighborIds[neighStart3 + k],
                        flatNeighborIds[neighStart2 + k],
                        flatNeighborIds[neighStart1 + k],
                        flatNeighborIds[neighStart0 + k]
                    );
                    
                    // Gather neighbor spins
                    __m256d neighbor_spin_vec = _mm256_i64gather_pd(
                        spin_z, neighbor_idx_vec, 8
                    );
                    
                    // Gather coupling constants J
                    __m256d J_vec = _mm256_set_pd(
                        flatNeighborJs[neighStart3 + k],
                        flatNeighborJs[neighStart2 + k],
                        flatNeighborJs[neighStart1 + k],
                        flatNeighborJs[neighStart0 + k]
                    );
                    
                    // Accumulate: h_local += J * neighbor_spin
                    h_local_vec = _mm256_fmadd_pd(J_vec, neighbor_spin_vec, h_local_vec);
                }
                
                // Compute deltaEnergy = 2 * old_sz * h_local
                __m256d deltaEnergy_vec = _mm256_mul_pd(
                    _mm256_mul_pd(_mm256_set1_pd(2.0), old_sz_vec),
                    h_local_vec
                );
                
                // Extract deltaEnergy to array for scalar exp
                alignas(32) Real deltaEnergy_arr[4];
                _mm256_store_pd(deltaEnergy_arr, deltaEnergy_vec);
                
                // Compute Boltzmann factors per lane (scalar)
                Real kT = this -> kb_ * T;
                Real boltzmann_arr[4];
                for (int lane = 0; lane < 4; ++lane) {
                    if (deltaEnergy_arr[lane] > 0.0) {
                        boltzmann_arr[lane] = std::exp(-deltaEnergy_arr[lane] / kT);
                    } else {
                        boltzmann_arr[lane] = 1.0;  // Always accept if energy decreases
                    }
                }
                
                // Load uniform random numbers and boltzmann factors
                __m256d uniform_vec = _mm256_set_pd(
                    uniform[3], uniform[2], uniform[1], uniform[0]
                );
                __m256d boltzmann_vec = _mm256_load_pd(boltzmann_arr);
                
                // Acceptance condition:
                // 1. deltaEnergy > 0 (energy increases)
                // 2. uniform > boltzmann (reject with probability 1 - exp(-ΔE/kT))
                __m256d zero = _mm256_setzero_pd();
                __m256d condition1 = _mm256_cmp_pd(deltaEnergy_vec, zero, _CMP_GT_OQ);
                __m256d condition2 = _mm256_cmp_pd(uniform_vec, boltzmann_vec, _CMP_GT_OQ);
                __m256d reject_mask = _mm256_and_pd(condition1, condition2);
                
                // Create mask for store: where reject_mask is true, we restore old_sz
                int mask = _mm256_movemask_pd(reject_mask);
                
                // Store proposed spins (will be overwritten for rejected lanes)
                _mm256_store_pd(&spin_z[soaIdx0], prop_sz_vec);
                
                // Restore old spins for rejected moves
                if (mask != 0) {
                    // Blend old and new based on mask
                    __m256d result_vec = _mm256_blendv_pd(prop_sz_vec, old_sz_vec, reject_mask);
                    _mm256_store_pd(&spin_z[soaIdx0], result_vec);
                    
                    // Update rejection counts per atom type
                    for (int lane = 0; lane < SIMD_LANES; ++lane) {
                        if (mask & (1 << lane)) {
                            Index soaIdx = (lane == 0) ? soaIdx0 :
                                          (lane == 1) ? soaIdx1 :
                                          (lane == 2) ? soaIdx2 : soaIdx3;
                            Index typeIdx = typeIndices[soaIdx];
                            this -> counterRejections_.at(typeIdx) += 1;
                        }
                    }
                }
                continue;  // Skip scalar loop
            }
#endif // defined(__AVX2__) && defined(VEGAS_AVX2_EXPERIMENTAL)
            
            // Scalar fallback for partial vectors or non-AVX2
            for (Index lane = 0; lane < laneCount; ++lane) {
                Index globalIdx = blockAtomIndices[blockPos + lane];
                Index soaIdx = globalToSoA[globalIdx];
                
                Index typeIdx = typeIndices[soaIdx];
                
                // Ising: flip spin direction
                Real old_sz = spin_z[soaIdx];
                Real prop_sz = -old_sz;  // Flip direction
                
                Index neighStart = soaNeighborOffsets[soaIdx];
                Index neighEnd = soaNeighborOffsets[soaIdx + 1];
                
                // Compute local field: h_local = Σ_{j neighbors} J * spin_z[nbh]
                Real h_local = 0.0;
                for (Index j = neighStart; j < neighEnd; ++j) {
                    Index nbh = flatNeighborIds[j];
                    Real J = flatNeighborJs[j];
                    h_local += J * spin_z[nbh];
                }
                
                // Energy change: ΔE = 2 * old_sz * h_local
                Real deltaEnergy = 2.0 * old_sz * h_local;
                
                // Temporarily store proposed spin
                spin_z[soaIdx] = prop_sz;
                
                if (deltaEnergy > 0 && uniform[lane] > std::exp(- deltaEnergy / (this -> kb_ * T)))
                {
                    spin_z[soaIdx] = old_sz;  // Reject: restore old value
                    this -> counterRejections_.at(typeIdx) += 1;
                }
            }
        }
    }
}

void System::monteCarloStep(Real T, Real H)
{
    // Dispatch to specialized implementation based on model type
    if (this -> isIsing_) {
        monteCarloStepImpl<IsingTag>(T, H);
    } else {
        monteCarloStepImpl<HeisenbergTag>(T, H);
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

    // Determine starting index for resume
    Index startIndex = 0;
    if (resumeIndex_ > 0) {
        if (resumeIndex_ >= this -> temps_.size()) {
            std::cout << "Checkpoint shows simulation already complete. Nothing to do.\n";
            return;
        }
        startIndex = resumeIndex_ + 1;
        std::cout << "Resuming from checkpoint: Starting at Temperature Index " 
                  << startIndex << " (T = " << this -> temps_.at(startIndex) << ")\n";
    }

    for (Index index = startIndex; index < this -> temps_.size(); ++index)
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
        // Sync SoA to Atoms before HDF5 output
        this -> lattice_.syncSoAToAtoms();
        this -> reporter_.partial_report(enes, histMag_x, histMag_y, histMag_z, this -> lattice_, index);
        
        // Write checkpoint after each temperature point
        this -> writeCheckpoint(this -> outName_, index);


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
        Index size = this->lattice_.getSizesByIndex().at(i);
        if (size == 0) {
            // No atoms of this type, keep sigma unchanged
            i++;
            continue;
        }
        rejection = val / Real(size);
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
    STARTER::CHECKFILE(fileState);
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
        STARTER::CHECKFILE(fileName);
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

void System::writeCheckpoint(std::string filename, Index tempIndex)
{
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file_id < 0) {
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }
    
    if (file_id < 0) {
        throw vegas::HDF5Exception("Cannot create checkpoint file: " + filename);
    }
    
    // Check if checkpoint group exists and delete it first
    H5E_auto_t old_func;
    void* old_client_data;
    H5Eget_auto2(H5E_DEFAULT, &old_func, &old_client_data);
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
    
    hid_t group_id = H5Gopen2(file_id, "/checkpoint", H5P_DEFAULT);
    H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);
    
    if (group_id >= 0) {
        H5Gclose(group_id);
        H5Ldelete(file_id, "/checkpoint", H5P_DEFAULT);
    }
    
    // Create checkpoint group
    group_id = H5Gcreate2(file_id, "/checkpoint", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group_id < 0) {
        H5Fclose(file_id);
        throw vegas::HDF5Exception("Cannot create checkpoint group");
    }
    
    Index N = this -> lattice_.getNumAtoms();
    // Ensure separate arrays are up-to-date (interleaved arrays are source of truth)
    this -> lattice_.syncInterleavedToSeparate();
    hsize_t dims[1] = {static_cast<hsize_t>(N)};
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    
    // Write spin_x
    double* spin_x = const_cast<double*>(lattice_.getSpinX());
    hid_t dataset_id = H5Dcreate2(group_id, "spin_x", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, spin_x);
    H5Dclose(dataset_id);
    
    // Write spin_y
    double* spin_y = const_cast<double*>(lattice_.getSpinY());
    dataset_id = H5Dcreate2(group_id, "spin_y", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, spin_y);
    H5Dclose(dataset_id);
    
    // Write spin_z
    double* spin_z = const_cast<double*>(lattice_.getSpinZ());
    dataset_id = H5Dcreate2(group_id, "spin_z", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, spin_z);
    H5Dclose(dataset_id);
    
    // Write mapping arrays for checkerboard decomposition
    const std::vector<Index>& globalToSoA = lattice_.getGlobalToSoA();
    const std::vector<Index>& soaToGlobal = lattice_.getSoAToGlobal();
    const std::vector<Index>& atomColors = lattice_.getAtomColors();
    
    // Write globalToSoA
    dataset_id = H5Dcreate2(group_id, "globalToSoA", H5T_NATIVE_UINT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        std::cerr << "[Warning] Failed to create globalToSoA dataset" << std::endl;
    } else {
        H5Dwrite(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, globalToSoA.data());
        H5Dclose(dataset_id);
    }
    
    // Write soaToGlobal
    dataset_id = H5Dcreate2(group_id, "soaToGlobal", H5T_NATIVE_UINT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, soaToGlobal.data());
    H5Dclose(dataset_id);
    
    // Write atomColors
    dataset_id = H5Dcreate2(group_id, "atomColors", H5T_NATIVE_UINT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, atomColors.data());
    H5Dclose(dataset_id);
    
    H5Sclose(dataspace_id);
    
    // Write RNG state
    auto rngState = this -> engine_.getState();
    hsize_t rng_dims[1] = {4};
    hid_t rng_space = H5Screate_simple(1, rng_dims, NULL);
    dataset_id = H5Dcreate2(group_id, "rng_state", H5T_NATIVE_UINT64, rng_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, rngState.data());
    H5Dclose(dataset_id);
    H5Sclose(rng_space);
    
    // Write SIMD RNG state
    auto simdRngState = this -> simdEngine_.getState();
    hsize_t simd_dims[2] = {SIMD_LANES, 4};
    hid_t simd_space = H5Screate_simple(2, simd_dims, NULL);
    dataset_id = H5Dcreate2(group_id, "simd_rng_state", H5T_NATIVE_UINT64, simd_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(simdRngState[0][0]));
    H5Dclose(dataset_id);
    H5Sclose(simd_space);
    
    // Write temperature index
    hsize_t ti_dims[1] = {1};
    hid_t ti_space = H5Screate_simple(1, ti_dims, NULL);
    dataset_id = H5Dcreate2(group_id, "temperature_index", H5T_NATIVE_INT, ti_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    int tempIdx = static_cast<int>(tempIndex);
    H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tempIdx);
    H5Dclose(dataset_id);
    H5Sclose(ti_space);
    
    // Write sigma values
    hsize_t sigma_dims[1] = {static_cast<hsize_t>(sigma_.size())};
    hid_t sigma_space = H5Screate_simple(1, sigma_dims, NULL);
    dataset_id = H5Dcreate2(group_id, "sigma", H5T_NATIVE_DOUBLE, sigma_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, sigma_.data());
    H5Dclose(dataset_id);
    H5Sclose(sigma_space);
    
    H5Gclose(group_id);
    H5Fclose(file_id);
}

bool System::loadCheckpoint(std::string filename, Index& tempIndex)
{
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        return false;
    }
    
    // Check if checkpoint group exists
    hid_t group_id = H5Gopen2(file_id, "/checkpoint", H5P_DEFAULT);
    if (group_id < 0) {
        H5Fclose(file_id);
        return false;
    }
    
    Index N = this -> lattice_.getNumAtoms();
    
    // Validate checkpoint dimensions (rank and size)
    {
        // Check spin_x dimensions
        hid_t check_dataset = H5Dopen2(group_id, "spin_x", H5P_DEFAULT);
        hid_t check_space = H5Dget_space(check_dataset);
        int rank = H5Sget_simple_extent_ndims(check_space);
        if (rank != 1) {
            H5Sclose(check_space);
            H5Dclose(check_dataset);
            H5Gclose(group_id);
            H5Fclose(file_id);
            throw vegas::HDF5Exception("Checkpoint spin_x is not a 1D array");
        }
        hsize_t dims[1];
        H5Sget_simple_extent_dims(check_space, dims, NULL);
        if (dims[0] != N) {
            H5Sclose(check_space);
            H5Dclose(check_dataset);
            H5Gclose(group_id);
            H5Fclose(file_id);
            throw vegas::SimulationException("Checkpoint size mismatch: checkpoint has " + 
                std::to_string(dims[0]) + " atoms but lattice has " + std::to_string(N));
        }
        H5Sclose(check_space);
        H5Dclose(check_dataset);
        
        // Check sigma dimensions
        check_dataset = H5Dopen2(group_id, "sigma", H5P_DEFAULT);
        check_space = H5Dget_space(check_dataset);
        rank = H5Sget_simple_extent_ndims(check_space);
        if (rank != 1) {
            H5Sclose(check_space);
            H5Dclose(check_dataset);
            H5Gclose(group_id);
            H5Fclose(file_id);
            throw vegas::HDF5Exception("Checkpoint sigma is not a 1D array");
        }
        H5Sget_simple_extent_dims(check_space, dims, NULL);
        if (dims[0] != sigma_.size()) {
            H5Sclose(check_space);
            H5Dclose(check_dataset);
            H5Gclose(group_id);
            H5Fclose(file_id);
            throw vegas::SimulationException("Checkpoint sigma size mismatch");
        }
        H5Sclose(check_space);
        H5Dclose(check_dataset);
    }
    
    // Read spin_x
    hid_t dataset_id = H5Dopen2(group_id, "spin_x", H5P_DEFAULT);
    hid_t dataspace_id = H5Dget_space(dataset_id);
    H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, lattice_.getSpinX());
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    
    // Read spin_y
    dataset_id = H5Dopen2(group_id, "spin_y", H5P_DEFAULT);
    dataspace_id = H5Dget_space(dataset_id);
    H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, lattice_.getSpinY());
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    
    // Read spin_z
    dataset_id = H5Dopen2(group_id, "spin_z", H5P_DEFAULT);
    dataspace_id = H5Dget_space(dataset_id);
    H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, lattice_.getSpinZ());
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    
    // Copy separate arrays to interleaved arrays (interleaved arrays are source of truth)
    this -> lattice_.syncSeparateToInterleaved();
    
    // Read mapping arrays if present (checkpoint version 2)
    const std::vector<Index>& latticeGlobalToSoA = lattice_.getGlobalToSoA();
    const std::vector<Index>& latticeSoAToGlobal = lattice_.getSoAToGlobal();
    const std::vector<Index>& latticeAtomColors = lattice_.getAtomColors();
    
    // Helper to check dataset existence
    auto dataset_exists = [&](const char* name) -> bool {
        return H5Lexists(group_id, name, H5P_DEFAULT) > 0;
    };
    
    if (dataset_exists("globalToSoA") && dataset_exists("soaToGlobal") && dataset_exists("atomColors")) {
        // Read mapping arrays
        std::vector<Index> checkpointGlobalToSoA(N);
        dataset_id = H5Dopen2(group_id, "globalToSoA", H5P_DEFAULT);
        H5Dread(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, checkpointGlobalToSoA.data());
        H5Dclose(dataset_id);
        
        std::vector<Index> checkpointSoAToGlobal(N);
        dataset_id = H5Dopen2(group_id, "soaToGlobal", H5P_DEFAULT);
        H5Dread(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, checkpointSoAToGlobal.data());
        H5Dclose(dataset_id);
        
        std::vector<Index> checkpointAtomColors(N);
        dataset_id = H5Dopen2(group_id, "atomColors", H5P_DEFAULT);
        H5Dread(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, checkpointAtomColors.data());
        H5Dclose(dataset_id);
        
        // Verify mapping consistency
        if (checkpointGlobalToSoA != latticeGlobalToSoA ||
            checkpointSoAToGlobal != latticeSoAToGlobal ||
            checkpointAtomColors != latticeAtomColors) {
            throw vegas::SimulationException(
                "Checkpoint mapping arrays differ from current lattice coloring. "
                "Cannot resume with different graph coloring.");
        }
        // Mapping matches, proceed
    } else {
        // Old checkpoint (no mapping arrays). Assume identity mapping.
        // Since we have already reordered data, the checkpoint spins are in original order.
        // We need to reorder spins from original order to SoA order.
        // This is a rare case; for simplicity, we reject.
        throw vegas::SimulationException(
            "Old checkpoint format detected (no mapping arrays). "
            "Cannot resume with checkerboard decomposition. Please start a new simulation.");
    }
    
    // Read RNG state
    std::array<uint64_t, 4> rngState;
    dataset_id = H5Dopen2(group_id, "rng_state", H5P_DEFAULT);
    H5Dread(dataset_id, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, rngState.data());
    this -> engine_.setState(rngState);
    H5Dclose(dataset_id);
    
    // Read SIMD RNG state (checkpoint version 3)
    if (dataset_exists("simd_rng_state")) {
        std::array<std::array<uint64_t, 4>, SIMD_LANES> simdRngState;
        dataset_id = H5Dopen2(group_id, "simd_rng_state", H5P_DEFAULT);
        H5Dread(dataset_id, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(simdRngState[0][0]));
        this -> simdEngine_.setState(simdRngState);
        H5Dclose(dataset_id);
    } else {
        // Old checkpoint: reconstruct SIMD RNG from scalar RNG state
        vegas::Xoshiro256StarStar temp = this -> engine_; // copy
        std::array<std::array<uint64_t, 4>, SIMD_LANES> laneStates;
        for (size_t i = 0; i < SIMD_LANES; ++i) {
            laneStates[i] = temp.getState();
            temp.jump();
        }
        this -> simdEngine_.setState(laneStates);
    }
    
    // Read temperature index
    int tempIdx = 0;
    dataset_id = H5Dopen2(group_id, "temperature_index", H5P_DEFAULT);
    H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tempIdx);
    tempIndex = static_cast<Index>(tempIdx);
    H5Dclose(dataset_id);
    
    // Read sigma values
    dataset_id = H5Dopen2(group_id, "sigma", H5P_DEFAULT);
    H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, sigma_.data());
    H5Dclose(dataset_id);
    
    H5Gclose(group_id);
    H5Fclose(file_id);
    
    return true;
}
