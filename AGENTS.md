# VEGAS Project - Agent Instructions

## Project Overview

VEGAS (VEctor General Atomistic Simulator) is a Monte Carlo simulation package for magnetic materials. The project implements atomistic simulations with multiple spin models, magnetic interactions, and HDF5 data output.

## Current State (February 2026)

### ✅ **COMPLETED IMPROVEMENTS**

#### **Critical Fixes:**
1. **Circular Dependency**: Fixed between `system.h` and `starter.h` using forward declaration
2. **Unsafe `atof()`**: Replaced all instances with `safe_stod()` with proper error handling
3. **Floating Point Comparisons**: Implemented `fp_equal()` and `fp_not_equal()` with epsilon
4. **Division by Zero**: Prevented in sigma adjustment (`system.cc:232-244`)
5. **Magic Numbers**: Defined as constants in `params.h**
6. **HDF5 Resource Leaks**: Fixed in `reporter.cc` (implemented RAII, ensured all HDF5 handles are properly closed)
7. **Duplicate EXIT() Definitions**: Fixed ODR violation by creating `error.h` with inline function

#### **Build System:**
1. **CMake Configuration**: Works with system packages (libjsoncpp-dev, libhdf5-dev)
2. **HDF5 Linking**: Fixed to include both C and C++ libraries
3. **JSON Library**: Multiple detection methods (pkg-config, find_package, manual)
4. **Testing Framework**: Basic tests integrated with CMake
5. **Python Environment**: Set up with uv virtual environment and requirements.txt
6. **Validation Pipeline**: Complete end-to-end validation script (`validate.sh`)

#### **Code Quality:**
1. **Move Semantics**: Implemented for `Lattice` and `System` classes
2. **Input Validation**: Enhanced `CHECKFILE()` with security checks
3. **Error Handling**: Improved throughout codebase
4. **Documentation**: Comprehensive README, FORMATS.md, example config
5. **SpinModel Abstraction**: Created base class and concrete implementations (Heisenberg, Ising, QuantumIsing, Adaptive, Cone, HN)
6. **Integration Tests**: Added comprehensive tests for SpinModel factory and Atom-SpinModel integration
7. **End-to-End Validation**: Complete test system with minimal example and verification scripts
8. **Unit Tests Expanded**: Added comprehensive unit tests for Atom and Lattice classes
9. **Exception Hierarchy**: Implemented comprehensive exception classes, replaced all EXIT() calls with specific exceptions
10. **CLI Enhancement**: Added professional command-line interface with subcommands (run, analyze, validate, info), verbosity control, colored output, configuration overrides, and backward compatibility

#### **2026-02-20: Physical Validation Suite (v2.4.0)**

**Critical Bugs Fixed:**
1. **HDF5 Double-Close Bug**: Reporter class lacked move semantics - temporary objects destroyed HDF5 handles prematurely, causing zero data output. Fixed with proper move constructor/assignment.
2. **Lattice Bond Storage Bug**: Bonds stored unidirectionally (i→j only), breaking detailed balance. Fixed by storing bonds bidirectionally (i→j AND j→i).

**Physical Validation Suite:**
- 6 comprehensive benchmarks comparing against exact solutions
- B1: 1D Ising energy vs exact (1.7% error)
- B2: 2D Ising Tc vs Onsager (deviation 0.138)
- B3: Ferromagnet ground state (exact match)
- B4: Ergodicity test (both branches converge)
- B5: Sigma freeze verification (code inspected)
- B6: Heisenberg isotropy (all moments within tolerance)

**New Files:**
- `benchmarks/` directory with configs, lattices, and results
- `generate_lattices.py` - Lattice file generator with bidirectional bonds
- `run_all_benchmarks.py` - Benchmark runner with automatic analysis
- `VALIDATION_REPORT.md` - Comprehensive results

**VEGAS Convention Documented:**
- J > 0 = ferromagnetic (standard physics convention)
- J < 0 = antiferromagnetic

#### **2026-02-23: I/O Scaling & RNG Infrastructure (v2.3.2)**

**I/O Scaling Profiling:**
- Profiled HDF5 parallel write scaling across 1, 4, 8, 16 processes
- Results: 1 process baseline (0.213s), 4 processes (0.223s, ~95% efficiency), 8 processes (0.271s, ~79% efficiency), 16 processes (0.469s, ~45% efficiency)
- Conclusion: Process-level parallelism valid; 4-8 processes optimal

**RNG Infrastructure:**
- Created `include/random.h` with SplitMix64 and Xoshiro256StarStar implementations
- SplitMix64 for high-quality seeding from std::random_device
- Xoshiro256** for fast, high-quality random numbers in simulation hot path
- **INTEGRATED**: Updated atom.h, atom.cc, spin_model.h, spin_model.cc, system.h, system.cc

**Chaos Testing:**
- Tested error paths with malformed inputs (invalid JSON, missing files, corrupt data)
- RAII pattern verified robust - all resources properly cleaned up on error

**Files Changed:**
- include/random.h: New file (Xoshiro256** and SplitMix64)
- benchmarks/configs/io_scaling_*.json: New files (scaling test configs)
- system.h, system.cc: Updated to use Xoshiro256**
- spin_model.h, spin_model.cc: Updated method signatures
- atom.h, atom.cc: Updated method signatures
- test_integration.cc: Updated tests to use new RNG

#### **2026-02-23: SoA Refactor & Template Dispatch (v2.4.0)**

**SoA (Structure of Arrays) Refactor:**
- Added flat arrays for spins: `spin_x_`, `spin_y_`, `spin_z_`
- Added flat neighbor structure: `NeighborInteraction {id, J}` with offset indexing
- Implemented Jagged Array with AoS for neighbor data (single cache line per neighbor)
- Added sync methods: `syncAtomsToSoA()` and `syncSoAToAtoms()`
- Register rollback: Old spin stored in local variables instead of memory

**Template Dispatch:**
- Created `include/spin_model_tags.h` with `HeisenbergTag` and `IsingTag`
- Added model detection: System detects Ising vs Heisenberg at initialization
- Implemented specialized Monte Carlo steps:
  - `monteCarloStepImpl<HeisenbergTag>`: Full 3D dot product
  - `monteCarloStepImpl<IsingTag>`: Z-component only, simple spin flip
- Compiler generates optimized code for each model type

**Performance Results:**
| Benchmark | Before | After | Improvement |
|-----------|--------|-------|-------------|
| Ising 400 atoms | 21.2s | 8.3s | 60% faster |
| Ising 100 atoms | 4.8s | 1.6s | 67% faster |
| Heisenberg 100 atoms | 4.7s | 4.8s | Same |

**Files Changed:**
- include/lattice.h: Added SoA arrays, NeighborInteraction struct
- src/lattice.cc: Populate flat neighbor arrays
- include/system.h: Added isIsing_ member, template declarations
- src/system.cc: Template implementations for dispatch
- include/spin_model_tags.h: New file

#### **2026-02-23: Sanitizer Integration & Documentation Fixes (v2.3.1)**

**Critical Fixes:**
1. **Exchange Convention Documentation**: Removed erroneous "inverted" references from all documentation (AGENTS.md, README.md, FORMATS.md, benchmarks/). The code uses standard physics convention: E = -J Σ S_i · S_j
2. **ASan/UBSan Integration**: Added CMake options for AddressSanitizer and UndefinedBehaviorSanitizer with HDF5 suppressions
3. **Sanitizer Testing**: Verified VEGAS passes ASan/UBSan on multiple benchmarks (B1, B3, B6) with no memory leaks

**Infrastructure Added:**
- `sanitizer_suppressions.txt`: Suppressions for known HDF5 false positives
- CMake options: `-DENABLE_ASAN=ON -DENABLE_UBSAN=ON`

**Testing Results (ASan enabled):**
- B1 (1D Ising): PASS, no leaks detected
- B3 (Ferromagnet): PASS, no leaks detected
- B6 (Heisenberg): PASS, no leaks detected
- vegas_simple_tests: PASS
- vegas_integration_tests: PASS

**Files Changed:**
- CMakeLists.txt: Added sanitizer options
- sanitizer_suppressions.txt: New file
- AGENTS.md, README.md, FORMATS.md: Documentation fixes
- generate_lattices.py: Fixed J sign conventions
- test_system/ferromagnetic_chain.txt: Fixed J value

#### **2026-02-20: Code Quality Overhaul (v2.1.0)**

**Critical Issues Fixed:**
1. **Dangling Pointer**: Replaced raw `Atom*` pointers with `Index` indices in neighbor storage (`atom.h`, `atom.cc`, `lattice.cc`)
2. **HDF5 Handle Leaks**: Fixed resource cleanup on error paths in Reporter constructor
3. **Undefined Methods**: Implemented `getAnisotropyUnit()`, `getTypeAnisotropy()`, `getKan()` in Atom class

**God Functions Refactored:**
1. **CREATE_SYSTEM()**: Now uses `ConfigParser` and `SystemBuilder` (~60% code reduction)
2. **Reporter()**: Extracted into 6 helper methods (~76% reduction in constructor)
3. **main()**: Deduplicated override warning code into `warn_about_overrides()` helper

**New Infrastructure:**
- `ConfigParser`: Clean JSON configuration parsing with validation
- `SystemBuilder`: Builder pattern for fluent System construction
- `SimulationConfig`: Structured configuration data struct
- Exception hierarchy: `VEGASException`, `FileIOException`, `InvalidInputException`, `ConfigurationException`, `NumericConversionException`, `SimulationException`, `HDF5Exception`

**Python Improvements:**
- Context managers for HDF5 file handling (eliminates resource leaks)
- Type hints added to all analyzers
- Safe division with zero temperature check
- Improved error handling throughout

**Files Changed:** 17 modified, 6 new files (+1,645 lines, -973 lines)

### **Files Changed:**
- CMakeLists.txt: Added sanitizer options
- sanitizer_suppressions.txt: New file
- AGENTS.md, README.md, FORMATS.md: Documentation fixes
- generate_lattices.py: Fixed J sign conventions
- test_system/ferromagnetic_chain.txt: Fixed J value

### **Known Issues:**

#### **Build/Configuration:**
1. **Conan Not Required**: Project uses system packages by default

#### **Code Quality (To Do):**
1. **`Atom` Class**: Still large (~320 lines), could benefit from further refactoring
2. **RNG Integration**: Xoshiro256** implemented but not yet integrated into codebase

## Build Instructions

### **Using System Packages (Recommended):**
```bash
# Install dependencies
sudo apt-get install libjsoncpp-dev libhdf5-dev libhdf5-cpp-103-1t64

# Build
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4

# Run tests
./vegas_simple_tests
```

### **Using Conan (Alternative):**
```bash
# Install Conan
pip install conan

# Build with Conan
conan profile detect --force
conan install . --output-folder=build --build=missing
cmake -S . -B build -G "Unix Makefiles" \
    -DCMAKE_TOOLCHAIN_FILE=build/conan_toolchain.cmake \
    -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

## Testing

### **Available Tests:**
1. **Simple Tests**: `./vegas_simple_tests` - Basic functionality tests (6 tests)
2. **Integration Tests**: `./vegas_integration_tests` - SpinModel factory and Atom integration (4 tests)
3. **ConfigParser Tests**: `./vegas_config_parser_tests` - Configuration parsing (8 tests)
4. **SystemBuilder Tests**: `./vegas_system_builder_tests` - Builder pattern tests (8 tests)
5. **System Tests**: `./vegas_system_tests` - System class tests (8 tests)
6. **Reporter Tests**: `./vegas_reporter_tests` - HDF5 output tests (6 tests)
7. **Validation Test**: `./validate.sh` - Complete end-to-end pipeline validation
8. **Google Tests**: If GTest found, `./vegas_tests` (not built by default)

### **Test Coverage Areas:**
- Floating point comparisons (`fp_equal`, `fp_not_equal`)
- Constants and type definitions
- Array operations
- Safe string-to-double conversion
- SpinModel factory and concrete implementations
- Atom-SpinModel integration
- HDF5 output structure validation
- Temperature validation (MIN_TEMPERATURE enforcement)
- Thermalization phases (thermalization vs measurement)
- Sigma adaptation timing

## Code Structure

### **Key Files:**
- `include/params.h` - Type definitions, constants, utility functions
- `include/system.h` - Main simulation class
- `include/atom.h` - Atom class (needs refactoring)
- `include/lattice.h` - Lattice class
- `include/reporter.h` - HDF5 output class
- `include/starter.h` - JSON configuration parser
- `include/spin_model.h` - Spin model abstraction

### **Recent Changes (v2.3.0):**
1. **`params.h`**: Added `MIN_TEMPERATURE`, `DEFAULT_THERMALIZATION_FRACTION` constants
2. **`system.h`**: Added `thermalizationSteps_`, `measurementSteps_`, `thermalizationFraction_`, `adaptSigma()`, `resetSigma()`
3. **`system.cc`**: Separated thermalization/measurement phases, sigma only adapts during thermalization
4. **`spin_model.cc`**: Fixed Cone/HN detailed balance (Rodrigues rotation), Heisenberg uses sigma, added rotation helpers
5. **`config_parser.cc`**: Added temperature validation, rejects T ≤ MIN_TEMPERATURE
6. **`vegas-analyzer-heisenberg.py`**: Fixed susceptibility (χ = ⟨M²⟩-⟨M⟩²/kT) and Cv formulas
7. **`test_config_parser.cc`**: Added T=0, T=-1.0 validation tests
8. **New test configs**: `zero_temp_config.json`, `negative_temp_config.json`

## Development Guidelines

### **Coding Standards:**
1. **Use Constants**: Never use magic numbers, use constants from `params.h`
2. **Floating Point**: Always use `fp_equal()`/`fp_not_equal()` for comparisons
3. **Error Handling**: Use `safe_stod()` for string-to-double conversion
4. **Input Validation**: Validate all inputs, check file paths for security
5. **Memory Safety**: Prefer RAII, avoid raw pointers where possible

### **Security Considerations:**
1. **File Paths**: Check for path traversal (`..`) in `CHECKFILE()`
2. **Input Validation**: Validate JSON inputs, numeric ranges
3. **Error Messages**: Don't expose system details in error messages

### **Performance:**
1. **Move Semantics**: Implement for large objects
2. **Avoid Copies**: Pass by const reference when appropriate
3. **Vectorization**: Use `std::valarray` operations where possible

### **Scientific Correctness:**
1. **Thermalization**: First 20% of MCS discarded by default (`DEFAULT_THERMALIZATION_FRACTION`)
2. **Sigma Adaptation**: Only during thermalization, fixed during measurement
3. **Temperature Validation**: Must be > `MIN_TEMPERATURE` (1e-10)
4. **Detailed Balance**: All spin proposals must be symmetric around current spin
5. **Thermodynamic Formulas**: Susceptibility and Cv must include kb factor

## Scientific Validation

### **Physical Correctness Guarantees**

The following scientific correctness properties are verified:

#### **Thermalization**
- First 20% of MCS discarded by default (`DEFAULT_THERMALIZATION_FRACTION = 0.2`)
- Thermalization implemented in C++ code, not post-processing
- HDF5 output contains only measurement phase data

#### **Detailed Balance**
- Metropolis acceptance: p = min(1, exp(-ΔE/kT)) correctly implemented
- All spin proposals are symmetric around current spin
- Cone/HN models use Rodrigues rotation to center proposals

#### **Sigma Adaptation**
- Only adapts during thermalization phase
- Fixed during measurement (preserves detailed balance)
- Reset to MAX_SIGMA at each temperature/field point

#### **Thermodynamic Formulas**
- Susceptibility: χ = (⟨M²⟩ - ⟨M⟩²) / (kT) with kb factor
- Specific Heat: Cv = σ²E / (kT²) with kb factor
- Temperature validation: T > MIN_TEMPERATURE (1e-10)

### **Testing**
- Temperature validation: Tests for T=0 (rejected), T=-1.0 (rejected), T=0.001 (accepted)
- HDF5 dimensions: Verified output is measurement_steps, not mcs
- All spin models: Verified norm preservation and proposal symmetry

## Next Development Phases

### **Phase 1: Critical Fixes (High Priority)**
1. **Fix HDF5 Resource Leaks**: Implement RAII wrappers, ensure all HDF5 handles are properly closed (Completed)
2. **Fix Duplicate EXIT() Definitions**: Resolve ODR violation, centralize error handling (Completed)
3. **Update Dockerfile**: Support both system packages and Conan, add Python dependencies (Completed)
4. **Improve Python Analyzers**: Add type hints, better error handling, dependency management

### **Phase 2: Code Quality & Testing (Medium Priority)**
1. **Replace EXIT() with Exception Hierarchy**: Create proper exception classes, convert all EXIT() calls (Completed)
2. **Expand Test Coverage**: Unit tests added for Atom, Lattice, System, Reporter, ConfigParser, and SystemBuilder (40 tests total)
3. **Refactor Atom Class**: Further extract spin logic to SpinModel, improve encapsulation
4. **Add Integration Tests**: Test full pipeline with various configurations (Completed)

### ✅ **Phase 2.5: CLI Enhancement (Completed)**
1. **Argument Parser**: Implemented proper command-line argument parsing with subcommands using cxxopts
2. **Enhanced Help System**: Comprehensive help with examples and configuration options
3. **Verbosity Control**: Added --verbose, --quiet, and --no-color options
4. **Configuration Overrides**: Fully implemented --mcs, --seed, --kb, --output, --sample, --initialstate
5. **Validate Subcommand**: Configuration validation with verbose output using ConfigParser
6. **Analyze Subcommand**: Python analyzer integration with venv detection
7. **Progress Display**: Improved progress indicators with optional detailed output (future)
8. **Version Information**: Added --version flag and build information

### **Phase 3: Features & Performance (Low Priority)**
1. **Parallelization**: Add OpenMP/MPI support for temperature/field point parallelism
2. **Performance Optimization**: Profile and optimize hot paths, improve cache locality
3. **Checkpoint/Restart**: Add simulation checkpointing for long runs
4. **GUI/Web Interface**: Create user-friendly interface for configuration and visualization

## Agent Workflow

### **When Starting Work:**
1. Read this AGENTS.md file
2. Check `git status` for uncommitted changes
3. Verify build works: `cd build && make clean && make`
4. Run tests: `./vegas_simple_tests` and `./vegas_integration_tests`
5. Run validation: `./validate.sh` (optional but recommended)

### **Before Making Changes:**
1. Understand the architecture (see README.md)
2. Check for existing constants in `params.h`
3. Review similar code patterns in the codebase
4. Consider backward compatibility

### **After Making Changes:**
1. Test compilation: `cd build && make`
2. Run basic tests: `./vegas_simple_tests`
3. Run integration tests: `./vegas_integration_tests`
4. Test main executable: `./vegas --help`
5. Consider adding/updating tests
6. Run validation script if changes affect pipeline: `./validate.sh`

### **Commit Guidelines:**
1. Use descriptive commit messages
2. Reference issue numbers if applicable
3. Keep changes focused and atomic
4. Update documentation if needed

## Troubleshooting

### **Common Issues:**

#### **Build Failures:**
```
H5Include.h not found
```
**Solution**: Use `hdf5.h` instead, ensure HDF5 development packages installed

```
json/json.h not found
```
**Solution**: Install `libjsoncpp-dev` package

#### **Linker Errors:**
```
undefined reference to H5*
```
**Solution**: Ensure linking with both HDF5 C and C++ libraries

#### **Runtime Errors:**
```
Invalid number format
```
**Solution**: Input validation failed, check data files

### **Debugging:**
1. **Compile with Debug**: `cmake .. -DCMAKE_BUILD_TYPE=Debug`
2. **Use GDB**: `gdb ./vegas`
3. **Check Logs**: Error messages go to stderr with color coding
4. **Test Components**: Use simple tests to isolate issues

## Dependencies

### **Required:**
- C++17 compiler
- CMake 3.20+
- HDF5 1.10+ (C and C++ libraries)
- jsoncpp 1.9.5+

### **Optional:**
- Google Test (for advanced testing)
- Python 3.8+ with h5py, matplotlib (for analyzers)
- Conan 2.0+ (alternative dependency management)

### **System Packages (Ubuntu/Debian):**
```bash
sudo apt-get install \
    build-essential \
    cmake \
    libjsoncpp-dev \
    libhdf5-dev \
    libhdf5-cpp-103-1t64
```

## Contact & Resources

- **Repository**: https://github.com/jdalzatec/vegas
- **Documentation**: README.md, AGENTS.md, FORMATS.md, WORKFLOW.md
- **Examples**: `example_config.json`, `test_system/` directory
- **Tests**: `tests/` directory, `validate.sh` validation script

## Version History

### **2026-02-23**: I/O Scaling & RNG Infrastructure (v2.3.2)

**I/O Scaling Profiling:**
- Profiled HDF5 parallel write scaling across 1, 4, 8, 16 processes
- Results: 1 process (0.213s), 4 processes (0.223s, ~95%), 8 processes (0.271s, ~79%), 16 processes (0.469s, ~45%)
- Conclusion: Process-level parallelism valid; 4-8 processes optimal

**RNG Infrastructure:**
- Created `include/random.h` with SplitMix64 and Xoshiro256StarStar implementations
- SplitMix64 for high-quality seeding from std::random_device
- Xoshiro256** for fast, high-quality random numbers in simulation hot path
- Not yet integrated (requires updates to atom.h, atom.cc, spin_model.h, spin_model.cc, system.h, system.cc)

**Chaos Testing:**
- Tested error paths with malformed inputs (invalid JSON, missing files, corrupt data)
- RAII pattern verified robust - all resources properly cleaned up on error

**Files Changed:**
- include/random.h: New file
- benchmarks/configs/io_scaling_*.json: New files (scaling test configs)

### **2026-02-20**: Physical Validation Suite (v2.4.0)

**Critical Bugs Fixed:**
1. **HDF5 Double-Close Bug**: Reporter class lacked move semantics - temporary objects destroyed HDF5 handles prematurely, causing zero data output. Fixed with proper move constructor/assignment operator.
2. **Lattice Bond Storage Bug**: Bonds stored unidirectionally (i→j only), breaking detailed balance for Metropolis dynamics. Fixed by storing bonds bidirectionally (i→j AND j→i).

**Physical Validation Suite:**
- 6 comprehensive benchmarks comparing against exact solutions
- B1: 1D Ising energy vs exact (1.7% relative error)
- B2: 2D Ising Tc vs Onsager (deviation: 0.138 < 0.15)
- B3: Ferromagnet ground state (exact match)
- B4: Ergodicity test (both branches converge identically)
- B5: Sigma freeze verification (code inspection confirmed)
- B6: Heisenberg isotropy (all moments within tolerance)

**New Infrastructure:**
- `benchmarks/` directory with configs, lattices, and results
- `generate_lattices.py` - Lattice file generator with bidirectional bonds
- `run_all_benchmarks.py` - Benchmark runner with automatic analysis
- `benchmarks/VALIDATION_REPORT.md` - Comprehensive validation results

**VEGAS Convention Documented:**
- J > 0 = ferromagnetic (standard physics convention)
- J < 0 = antiferromagnetic

**Files Changed:** 31 files (+7,519 lines)

### **2026-02-20**: Scientific Correctness Overhaul (v2.3.0)

**Critical Physics Fixes (P0):**
1. **Temperature Validation**: Added `MIN_TEMPERATURE` constant (1e-10), ConfigParser rejects T ≤ MIN_TEMPERATURE with clear error message
2. **Thermalization Implementation**: First 20% of MCS (`DEFAULT_THERMALIZATION_FRACTION`) properly discarded - measurements only stored from equilibrium phase
3. **Sigma Adaptation Timing**: Sigma only adapts during thermalization phase, remains FIXED during measurement (preserves detailed balance)
4. **Cone/HN Detailed Balance**: Fixed violation - proposals now centered on current spin direction using Rodrigues rotation formula
5. **Heisenberg Sigma Usage**: Model now uses sigma parameter for Gaussian perturbation magnitude (previously ignored)

**Analyzer Fixes (P1):**
1. **Susceptibility Formula**: Corrected from std(|M|)²/T to χ = (⟨M²⟩ - ⟨M⟩²) / (kT) with proper kb factor
2. **Specific Heat Formula**: Corrected to Cv = σ²E / (kT²) with proper kb factor

**New Constants (`params.h`):**
- `MIN_TEMPERATURE = 1e-10` - Minimum allowed temperature
- `DEFAULT_THERMALIZATION_FRACTION = 0.2` - 20% thermalization by default

**New System Class Members:**
- `thermalizationSteps_` - Steps for equilibration
- `measurementSteps_` - Steps for data collection (mcs - thermalization)
- `thermalizationFraction_` - Configurable fraction
- `adaptSigma()` - Adaptation method (only called during thermalization)
- `resetSigma()` - Reset to MAX_SIGMA at each T/H point

**Breaking Changes:**
- ⚠️ **HDF5 output dimension change**: Arrays now store `[n_temps, measurement_steps]` instead of `[n_temps, mcs]`
- Existing analyzers that assume full mcs length will need adjustment
- Thermalization cutoff is now in C++ code, not post-processing

**Files Changed:**
- `include/params.h`, `include/system.h`
- `src/system.cc`, `src/config_parser.cc`, `src/spin_model.cc`
- `analyzers/vegas-analyzer-heisenberg.py`
- `tests/test_config_parser.cc`, `test_system/zero_temp_config.json`, `test_system/negative_temp_config.json`

### **2026-02-20**: CLI Enhancements (v2.2.0)
- Implemented `vegas validate` subcommand using ConfigParser with verbose output
- Implemented `vegas analyze` subcommand with Python analyzer integration and venv detection
- Implemented configuration overrides (--mcs, --seed, --kb, --output, --sample, --initialstate)
- Added `closed_` flag to Reporter class to prevent HDF5 double-close
- Added SimulationConfig::applyOverrides() method for clean override handling
- Updated report_overrides() to inform users when overrides are applied

### **2026-02-20**: Code Quality Overhaul (v2.1.0)
- Fixed dangling pointer in Atom neighbor storage (indices instead of pointers)
- Fixed HDF5 handle leaks on error paths in Reporter
- Implemented undefined Atom methods (getAnisotropyUnit, getTypeAnisotropy, getKan)
- Created ConfigParser and SystemBuilder infrastructure
- Refactored CREATE_SYSTEM() and Reporter() God functions
- Added comprehensive exception hierarchy
- Added context managers and type hints to Python analyzers
- Updated version to 2.1.0

### **2026-02-17**: SpinModel Abstraction and Validation Pipeline
- Implemented SpinModel base class with factory pattern
- Added comprehensive integration tests
- Created end-to-end validation system with test examples
- Set up Python environment with uv and requirements.txt
- Added workflow documentation and automated validation script
- Fixed move semantics issues in Atom and Lattice classes
- Fixed HDF5 resource leaks in `reporter.cc` (implemented RAII)
- Fixed duplicate `EXIT()` definitions (ODR violation) with `error.h`
- Implemented exception hierarchy and replaced all EXIT() calls with specific exceptions
- Updated Dockerfile to support system packages and Python dependencies
- Enhanced command-line interface with subcommands, verbosity control, and configuration override options

### **2025-02-17**: Major Code Quality Improvements
- Fixed critical vulnerabilities (atof, division by zero, etc.)
- Improved build system and dependency management
- Added testing framework
- Enhanced documentation
- Implemented move semantics and better error handling

### **Previous**: Initial implementation
- Monte Carlo simulation engine
- Multiple spin models
- HDF5 data output
- Python analysis tools

---

*Last updated: February 2026*  
*Maintained by: VEGAS Development Team*