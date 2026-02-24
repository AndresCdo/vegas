# VEGAS: Monte Carlo Simulations for Magnetic Materials

VEGAS (VEctor General Atomistic Simulator) is a software package for simulation, graphics, and analysis tools for atomistic simulations of magnetic materials using Monte Carlo methods.

## Features

- **Multiple Spin Models**: Heisenberg, Ising, QuantumIsing, Adaptive, Cone, HN, and custom spin models
- **Monte Carlo Simulations**: Metropolis algorithm with adaptive step size
- **Magnetic Interactions**: Exchange, anisotropy, Zeeman, and Dzyaloshinskii-Moriya interactions
- **Data Output**: HDF5 format with compression and chunking
- **Analysis Tools**: Python scripts for data analysis and visualization
- **Validation Pipeline**: Automated end-to-end testing with sample system
- **Professional CLI**: Modern command-line interface with subcommands, verbosity control, and configuration options
- **Scientifically Correct**: Proper thermalization, detailed balance, and validated thermodynamic formulas
- **Parallelization**: Trivially parallelizable across temperature/field sweeps via independent processes
- **Sanitizer Support**: ASan/UBSan integration for memory safety verification
- **Physical Validation**: Comprehensive benchmark suite with exact solutions
- **Advanced Optimizations**: Checkerboard decomposition, SIMD‑ready block reordering, RNG pre‑generation (9‑19% speedup)
- **High Performance**: SoA data layout, template dispatch, Xoshiro256** RNG (60% Ising speedup), checkerboard decomposition with block reordering (9‑19% speedup)

## Installation

### Prerequisites

- C++17 compatible compiler (GCC 8+, Clang 10+, MSVC 2019+)
- CMake 3.20 or higher
- Python 3.8+ (for analysis tools)
- **Optional**: Conan 2.0+ (alternative dependency management)

### Build Instructions

#### **Using System Packages (Recommended)**

```bash
# Install dependencies (Ubuntu/Debian)
sudo apt-get install libjsoncpp-dev libhdf5-dev libhdf5-cpp-103-1t64

# Build
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4

# Run tests
./vegas_simple_tests
```

#### **Using Sanitizers (Debug/Development)**

```bash
# Build with AddressSanitizer and UndefinedBehaviorSanitizer
mkdir -p build_asan && cd build_asan
cmake .. -DCMAKE_BUILD_TYPE=Debug -DENABLE_ASAN=ON -DENABLE_UBSAN=ON
make -j4

# Run tests (suppressions needed for HDF5 false positives)
ASAN_OPTIONS=suppressions=../sanitizer_suppressions.txt \
LSAN_OPTIONS=suppressions=../sanitizer_suppressions.txt \
./vegas_simple_tests
```

#### **Using Conan (Alternative)**

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

# Run tests
cd build && ctest --output-on-failure
```

### Docker Build

The Dockerfile uses system packages (Ubuntu 24.04) and includes Python dependencies for the analyzers:

```bash
docker build -t vegas .
docker run -v $(pwd)/data:/data vegas
```

The container's entrypoint is set to the `vegas` executable, so you can run simulations directly:

```bash
docker run -v $(pwd)/config.json:/config.json vegas run /config.json
```

## Usage

### Command Line Interface

VEGAS provides a modern command-line interface with subcommands, verbosity control, and configuration options. The interface is built using the `cxxopts` library and supports both traditional single-argument usage and new subcommand-based usage.

#### Quick Examples

```bash
# Run a simulation
./vegas run config.json

# Analyze simulation results (future integration)
./vegas analyze output.h5

# Validate configuration file (future integration)
./vegas validate config.json

# Show version and information
./vegas info

# Backward compatibility mode (deprecated but supported)
./vegas config.json

# Show help
./vegas --help

# Show version
./vegas --version

# Run quietly (minimal output)
./vegas --quiet run config.json

# Disable colored output
./vegas --no-color run config.json
```

#### Subcommands

| Subcommand | Description | Required Argument |
|------------|-------------|-------------------|
| `run` | Run a Monte Carlo simulation | Configuration file (JSON) |
| `analyze` | Analyze simulation results (future) | HDF5 output file |
| `validate` | Validate configuration file (future) | Configuration file (JSON) |
| `info` | Show version and system information | None |

#### Global Options

These options can be used with any subcommand or with backward-compatible usage:

| Option | Short | Description |
|--------|-------|-------------|
| `--help` | `-h` | Show help message and exit |
| `--version` | `-v` | Show version information and exit |
| `--verbose` | `-V` | Verbose output (default) |
| `--quiet` | `-q` | Quiet mode (minimal output) |
| `--no-color` | | Disable colored terminal output |
| `--config FILE` | `-c` | Configuration file (backward compatibility) |

#### Configuration Overrides (Future Implementation)

The following options can be used to override JSON configuration values when using the `run` subcommand. **Note**: These are currently stubbed with warning messages; full implementation is planned for a future release.

```bash
./vegas run config.json --mcs 10000 --seed 42 --output results.h5
```

| Option | Description | JSON Equivalent |
|--------|-------------|-----------------|
| `--mcs N` | Monte Carlo steps | `"mcs": N` |
| `--seed N` | Random seed | `"seed": N` |
| `--kb X` | Boltzmann constant | `"kb": X` |
| `--output FILE` | Output file path | `"out": "FILE"` |
| `--sample FILE` | Sample file path | `"sample": "FILE"` |
| `--initialstate FILE` | Initial state file | `"initialstate": "FILE"` |

### ⚠️ Breaking Changes (v2.3.0)

**HDF5 Output Dimension Change**

Starting with v2.3.0, HDF5 arrays store only the **measurement phase** (80% of MCS by default), not the full MCS. This ensures equilibrium statistics but requires updating existing analysis scripts.

**Before (v2.2.0):**
```python
# magnetization_z shape: [n_temps, mcs]  (e.g., [3, 100])
```

**After (v2.3.0):**
```python
# magnetization_z shape: [n_temps, measurement_steps]  (e.g., [3, 80])
```

**Migration Guide:**
```python
# Old code (v2.2.0)
mcs = data.attrs['mcs']  # Was: 100
tau = mcs // 5           # Was: 20

# New code (v2.3.0)
measurement_steps = data['energy'].shape[1]  # Now: 80
tau = measurement_steps // 4  # Equivalent thermalization cutoff
```

The `mcs` attribute in HDF5 now stores `measurement_steps`. To get original MCS:
```python
measurement_steps = data.attrs['mcs']
original_mcs = int(measurement_steps / 0.8)  # Divide by thermalization fraction
```

### Input JSON Format

VEGAS uses JSON configuration files. Here's a minimal example:

```json
{
  "sample": "my_sample",
  "mcs": 10000,
  "kb": 1.0,
  "temperature": {
    "start": 0.001,
    "final": 10.0,
    "points": 50
  },
  "field": {
    "start": 0.0,
    "final": 5.0,
    "points": 25
  },
  "out": "output.h5"
}
```

### Complete JSON Schema

| Section | Key | Type | Default | Description |
|---------|-----|------|---------|-------------|
| General | `sample` | string | required | Sample name |
| | `mcs` | integer | 5000 | Monte Carlo steps |
| | `kb` | float | 1.0 | Boltzmann constant |
| | `out` | string | `sample.h5` | Output filename |
| Temperature | `temperature` | object/float | 0.0 | Temperature settings (must be > 1e-10) |
| | `temperature.start` | float | 0.001 | Start temperature |
| | `temperature.final` | float | 10.0 | Final temperature |
| | `temperature.points` | integer | 5 | Number of points |
| | `temperature.delta` | float | 0.1 | Temperature step |
| Field | `field` | object/float | 0.0 | Magnetic field settings |
| | `field.start` | float | 0.001 | Start field |
| | `field.final` | float | 10.0 | Final field |
| | `field.points` | integer | 5 | Number of points |
| | `field.delta` | float | 0.1 | Field step |
| Advanced | `seed` | integer | time(NULL) | Random seed |
| | `initialstate` | string | "" | Initial state file |
| | `anisotropy` | string/array | [] | Anisotropy files |

### Output Format

VEGAS outputs data in HDF5 format. **Note**: Only the measurement phase (80% of MCS) is stored; thermalization data is discarded.

```
/output.h5
├── attributes
│   ├── mcs: Measurement steps (NOT total MCS - thermalization discarded)
│   ├── seed: Random seed
│   └── kb: Boltzmann constant
├── datasets
│   ├── temperature: Temperature values [n_temps]
│   ├── field: Field values [n_temps]
│   ├── magnetization_x: X magnetization [n_temps, measurement_steps]
│   ├── magnetization_y: Y magnetization [n_temps, measurement_steps]
│   ├── magnetization_z: Z magnetization [n_temps, measurement_steps]
│   ├── energy: Total energy [n_temps, measurement_steps]
│   ├── positions: Atom positions [n_atoms, 3]
│   └── finalstates: Final spin states [n_temps, n_atoms, 3]
```

**Migration from v2.2.0:** See Breaking Changes section above.

### Analysis Tools

Python scripts for data analysis are located in the `analyzers/` directory. **Note**: The `vegas analyze` subcommand is planned for future integration with these Python analyzers, but currently the Python scripts must be run directly.

First install dependencies:

```bash
# Using uv (recommended)
uv sync

# Or using pip
pip install h5py matplotlib numpy click
```

Then run the analyzers:

```bash
# Basic analysis
python analyzers/vegas-analyzer-lite.py output.h5

# Heisenberg model analysis
python analyzers/vegas-analyzer-heisenberg.py output.h5

# XYZ trajectory analysis
python analyzers/vegas-analyzer-xyz.py output.h5
```

## Architecture

### Core Components

1. **Atom**: Represents a magnetic atom with position, spin, and interactions
2. **SpinModel**: Abstract interface for spin models (Heisenberg, Ising, etc.)
3. **Lattice**: Collection of atoms with neighbor relationships
4. **System**: Monte Carlo simulation engine
5. **Reporter**: Handles HDF5 data output
6. **Starter**: JSON configuration parser and initialization

### Code Structure

```
vegas/
├── include/           # Header files
│   ├── atom.h        # Atom class
│   ├── lattice.h     # Lattice class
│   ├── system.h      # System class
│   ├── params.h      # Type definitions and constants
│   ├── starter.h     # Configuration parser
│   └── spin_model.h  # Spin model abstraction
├── src/              # Source files
│   ├── main.cc       # Main entry point
│   ├── atom.cc       # Atom implementation
│   ├── lattice.cc    # Lattice implementation
│   ├── system.cc     # System implementation
│   ├── starter.cc    # Starter implementation
│   └── spin_model.cc # Spin model implementations
├── analyzers/        # Python analysis scripts
├── tests/            # Test files
└── CMakeLists.txt    # Build configuration
```

## Development

### Code Quality

- **Code Style**: Follows C++ Core Guidelines
- **Testing**: Google Test framework for unit tests
- **Static Analysis**: Clang-Tidy and cppcheck
- **Memory Safety**: RAII pattern, smart pointers where appropriate
- **Development Guidelines**: See [AGENTS.md](AGENTS.md) for detailed project state and development instructions

### Building and Testing

```bash
# Configure with tests enabled
cmake -S . -B build -DBUILD_TESTING=ON

# Build everything
cmake --build build -j

# Run tests
cd build && ctest --output-on-failure

# Run specific tests
./build/vegas_simple_tests
./build/vegas_integration_tests
```

### Validation

To verify the complete simulation pipeline works correctly, run the validation script:

```bash
./validate.sh
```

This will:
1. Run a simulation with the test system
2. Verify HDF5 output structure
3. Run the Python analyzer
4. Check all output files

See [WORKFLOW.md](WORKFLOW.md) for detailed workflow documentation.

### Adding New Features

1. **New Spin Models**: Implement the `SpinModel` interface
2. **New Interactions**: Add energy calculation methods to `Atom`
3. **New Output Formats**: Extend the `Reporter` class
4. **Performance Optimizations**: Use parallel algorithms where possible

## Physics

### Hamiltonian

The system Hamiltonian includes:

1. **Exchange Interaction**: \( H_{ex} = -\sum_{i<j} J_{ij} \mathbf{S}_i \cdot \mathbf{S}_j \)
2. **Zeeman Interaction**: \( H_z = -\sum_i \mathbf{H} \cdot \mathbf{S}_i \)
3. **Anisotropy**: \( H_{an} = -\sum_i K (\mathbf{n} \cdot \mathbf{S}_i)^2 \)
4. **Dzyaloshinskii-Moriya**: \( H_{dm} = -\sum_{i<j} \mathbf{D}_{ij} \cdot (\mathbf{S}_i \times \mathbf{S}_j) \)

### Monte Carlo Algorithm

The simulation follows a two-phase approach:

1. **Initialization**: Random or specified initial spin configuration

2. **Thermalization Phase** (first 20% of MCS):
   - Select random atom
   - Propose spin change (symmetric around current spin)
   - Calculate energy difference ΔE = E_new - E_old
   - Accept with probability p = min(1, exp(-ΔE/k_B T))
   - **Adapt sigma** based on acceptance rate (target ~50%)
   - Data NOT stored

3. **Measurement Phase** (remaining 80% of MCS):
   - Same Metropolis steps
   - **Sigma is FIXED** (no adaptation - preserves detailed balance)
   - Record magnetization and energy to HDF5

4. **Output**: HDF5 contains only measurement phase data

This approach ensures proper equilibrium sampling with fixed transition probabilities.

## Performance

### Optimization Features

- **Memory Layout**: Contiguous arrays for cache efficiency
- **Vectorization**: SIMD instructions for energy calculations
- **Parallelization**: Independent temperature/field points
- **I/O Optimization**: Chunked HDF5 writes with compression

### Benchmarks

Typical performance on a modern CPU (v2.5.0 with checkerboard decomposition):

| System | Model | Time | Steps/sec |
|--------|-------|------|-----------|
| 100 atoms | Ising | 1.6s | 31,250 |
| 100 atoms | Heisenberg | 4.8s | 10,400 |
| 400 atoms | Ising | 8.3s | 24,000 |
| 400 atoms | Heisenberg | ~21s | 9,500 |

**Performance Improvements (v2.5.0)**:
- Ising models: **19.4% faster** with block‑ordered SoA (cache locality)
- Heisenberg models: **9.4% faster** with block‑ordered SoA
- Graph coloring enables non‑bipartite lattices (triangular, Kagome)
- Per‑lane SIMD RNG eliminates false sharing
- Register rollback reduces memory writes

### I/O Scaling

HDF5 parallel write performance across multiple processes:

| Processes | Time (s) | Efficiency |
|-----------|----------|------------|
| 1         | 0.213    | 100%       |
| 4         | 0.223    | 95%        |
| 8         | 0.271    | 79%        |
| 16        | 0.469    | 45%        |

**Recommendation**: 4-8 processes optimal for process-level parallelism.

## Scientific Validation

### Physical Correctness

VEGAS v2.3.0 implements the following scientific correctness guarantees:

#### Thermalization
- **20% default cutoff** (`DEFAULT_THERMALIZATION_FRACTION`)
- Implemented in simulation, not post-processing
- Only equilibrium data stored

#### Detailed Balance
- Symmetric spin proposals for all models
- Cone/HN: Rodrigues rotation centers proposals on current spin
- Metropolis: p = min(1, exp(-ΔE/kT))

#### Sigma Adaptation
- **Thermalization only**: Adapts during equilibration
- **Fixed during measurement**: Preserves detailed balance
- Reset at each T/H point

#### Thermodynamic Formulas
- **Susceptibility**: χ = (⟨M²⟩ - ⟨M⟩²) / (k_B T)
- **Specific Heat**: C_V = (⟨E²⟩ - ⟨E⟩²) / (k_B T²)
- **Temperature**: Validated T > 10⁻¹⁰

### Validation Tests
```bash
# Run validation pipeline
./validate.sh

# Run unit tests
./build/vegas_simple_tests
./build/vegas_config_parser_tests  # Includes T=0, T=-1.0 validation
```

### Physical Validation Benchmarks

VEGAS v2.4.0 includes a comprehensive validation suite comparing simulations against exact solutions:

```bash
# Run all benchmarks
python run_all_benchmarks.py

# Results are stored in benchmarks/results/
# See benchmarks/VALIDATION_REPORT.md for full analysis
```

| Benchmark | Description | Result |
|-----------|-------------|--------|
| B1 | 1D Ising energy vs exact | 1.7% error (PASS) |
| B2 | 2D Ising Tc vs Onsager | deviation 0.138 (PASS) |
| B3 | Ferromagnet ground state | exact match (PASS) |
| B4 | Ergodicity test | branches converge (PASS) |
| B5 | Sigma freeze | code verified (PASS) |
| B6 | Heisenberg isotropy | moments correct (PASS) |

See [benchmarks/VALIDATION_REPORT.md](benchmarks/VALIDATION_REPORT.md) for detailed results.

## Contributing

1. Fork the repository
2. Create a feature branch
3. Write tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

### Code Review Checklist

- [ ] Code follows style guidelines
- [ ] Tests are included and pass
- [ ] Documentation is updated
- [ ] No performance regressions
- [ ] Security considerations addressed
- [ ] Validation script passes (`./validate.sh`)

## License

## Version History

| Version | Date | Description |
|---------|------|-------------|
| v2.5.0 | 2026-02-23 | **Checkerboard & Cache**: Graph coloring, block‑ordered SoA, per‑lane SIMD RNG, 9‑19% speedup via cache locality, "Gather Wall" identified |
| v2.4.0 | 2026-02-23 | SoA refactor, Template Dispatch, Xoshiro256** RNG, 60% Ising speedup |
| v2.3.2 | 2026-02-23 | I/O scaling profiling, RNG infrastructure |
| v2.3.1 | 2026-02-23 | ASan/UBSan integration, documentation fixes |
| v2.3.0 | 2026-02-20 | Physical validation suite, bug fixes |
| v2.2.0 | 2026-02-20 | CLI enhancements |
| v2.1.0 | 2026-02-20 | Code quality overhaul |

MIT License - see LICENSE file for details.

## Citation

If you use VEGAS in your research, please cite:

```bibtex
@software{vegas2026,
  title = {VEGAS: Vector General Atomistic Simulator},
  author = {VEGAS Development Team},
  year = {2026},
  url = {https://github.com/jdalzatec/vegas}
}
```

## Support

- **Issues**: GitHub Issues tracker
- **Documentation**: This README and code comments
- **Community**: GitHub Discussions

## Acknowledgments

- Developers and contributors
- Funding agencies supporting magnetic materials research
- Open source libraries: HDF5, jsoncpp, Google Test
