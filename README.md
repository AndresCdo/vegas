# VEGAS: Monte Carlo Simulations for Magnetic Materials

VEGAS (VEctor General Atomistic Simulator) is a software package for simulation, graphics, and analysis tools for atomistic simulations of magnetic materials using Monte Carlo methods.

## Features

- **Multiple Spin Models**: Heisenberg, Ising, and custom spin models
- **Monte Carlo Simulations**: Metropolis algorithm with adaptive step size
- **Magnetic Interactions**: Exchange, anisotropy, Zeeman, and Dzyaloshinskii-Moriya interactions
- **Data Output**: HDF5 format with compression and chunking
- **Analysis Tools**: Python scripts for data analysis and visualization
- **Parallelization**: Temperature and field point parallelism

## Installation

### Prerequisites

- C++17 compatible compiler (GCC 8+, Clang 10+, MSVC 2019+)
- CMake 3.20 or higher
- Conan 2.0 or higher
- Python 3.8+ (for analysis tools)

### Build Instructions

1. Install [Conan](https://conan.io/):

```bash
pip install conan
```

2. Build the project:

```bash
conan profile detect --force
conan install . --output-folder=build --build=missing
cmake -S . -B build -G "Unix Makefiles" \
    -DCMAKE_TOOLCHAIN_FILE=build/conan_toolchain.cmake \
    -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

3. Run tests (optional):

```bash
cd build
ctest --output-on-failure
```

### Docker Build

```bash
docker build -t vegas .
docker run -v $(pwd)/data:/data vegas
```

## Usage

### Basic Usage

```bash
./vegas input.json
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
| Temperature | `temperature` | object/float | 0.0 | Temperature settings |
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

VEGAS outputs data in HDF5 format with the following structure:

```
/output.h5
├── attributes
│   ├── mcs: Monte Carlo steps
│   ├── seed: Random seed
│   └── kb: Boltzmann constant
├── datasets
│   ├── temperature: Temperature values [n_temps]
│   ├── field: Field values [n_fields]
│   ├── magnetization_x: X magnetization [n_temps, n_fields, mcs]
│   ├── magnetization_y: Y magnetization [n_temps, n_fields, mcs]
│   ├── magnetization_z: Z magnetization [n_temps, n_fields, mcs]
│   ├── energy: Total energy [n_temps, n_fields, mcs]
│   └── positions: Atom positions [n_atoms, 3]
```

### Analysis Tools

Python scripts for data analysis are located in the `analyzers/` directory:

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
2. **Lattice**: Collection of atoms with neighbor relationships
3. **System**: Monte Carlo simulation engine
4. **Reporter**: Handles HDF5 data output
5. **Starter**: JSON configuration parser and initialization

### Code Structure

```
vegas/
├── include/           # Header files
│   ├── atom.h        # Atom class
│   ├── lattice.h     # Lattice class
│   ├── system.h      # System class
│   ├── params.h      # Type definitions and constants
│   └── starter.h     # Configuration parser
├── src/              # Source files
│   ├── main.cc       # Main entry point
│   ├── atom.cc       # Atom implementation
│   ├── lattice.cc    # Lattice implementation
│   ├── system.cc     # System implementation
│   └── starter.cc    # Starter implementation
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

### Building and Testing

```bash
# Configure with tests enabled
cmake -S . -B build -DBUILD_TESTING=ON

# Build everything
cmake --build build -j

# Run tests
cd build && ctest --output-on-failure

# Run specific test
./build/vegas_simple_tests
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

1. **Initialization**: Random or specified initial spin configuration
2. **Metropolis Step**:
   - Select random atom
   - Propose spin change
   - Calculate energy difference ΔE
   - Accept with probability \( p = \min(1, \exp(-ΔE/k_B T)) \)
3. **Adaptive Step Size**: Sigma adjustment based on acceptance rate
4. **Measurement**: Record magnetization and energy every N steps
5. **Thermalization**: First 20% of steps discarded for equilibration

## Performance

### Optimization Features

- **Memory Layout**: Contiguous arrays for cache efficiency
- **Vectorization**: SIMD instructions for energy calculations
- **Parallelization**: Independent temperature/field points
- **I/O Optimization**: Chunked HDF5 writes with compression

### Benchmarks

Typical performance on a modern CPU:
- 10,000 atoms: ~100,000 steps/second
- 100,000 atoms: ~10,000 steps/second
- Memory usage: ~100 MB per million atoms

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

## License

MIT License - see LICENSE file for details.

## Citation

If you use VEGAS in your research, please cite:

```bibtex
@software{vegas2025,
  title = {VEGAS: Vector General Atomistic Simulator},
  author = {VEGAS Development Team},
  year = {2025},
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
