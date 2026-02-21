# VEGAS End-to-End Workflow

This document describes the complete workflow for running simulations with VEGAS, from building the software to analyzing results.

## Prerequisites

1. **System Dependencies** (Ubuntu/Debian):
   ```bash
   sudo apt-get install build-essential cmake libjsoncpp-dev libhdf5-dev libhdf5-cpp-103-1t64
   ```

2. **Python Environment** (for analyzers):
   ```bash
   # Install uv (recommended) or use pip
   curl -LsSf https://astral.sh/uv/install.sh | sh
   uv sync  # Creates .venv with dependencies from requirements.txt
   ```

## Building VEGAS

```bash
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

## Testing the Build

```bash
# Run unit tests
./vegas_simple_tests

# Run integration tests
./vegas_integration_tests
```

## Running a Simulation

### 1. Prepare Input Files

A minimal test system is provided in `test_system/`:

- **Lattice file** (`ferromagnetic_chain.txt`): Defines atoms, interactions, and anisotropy terms
- **Initial spins** (`initial_spins.txt`): Initial spin configuration (one vector per line)
- **Configuration** (`test_config.json`): Simulation parameters (MCS, temperature range, field range)

### 2. Run Simulation

```bash
./vegas ../test_system/test_config.json
```

This produces an HDF5 output file (`test_results.h5` by default).

### 3. Verify Output

Use the provided Python script to validate HDF5 structure:

```bash
cd ../test_system
../.venv/bin/python verify_hdf5.py test_results.h5
```

Expected output shows:
- All required datasets present (temperature, field, magnetization_x/y/z, energy, positions, types)
- Correct dimensions matching configuration
- Sample values within expected ranges

## Analyzing Results

VEGAS includes Python analyzers in the `analyzers/` directory:

```bash
# Using the lightweight analyzer (produces mean values and PDF plot)
../.venv/bin/python ../analyzers/vegas-analyzer-lite.py test_results.h5

# Using the Heisenberg analyzer (more detailed analysis)
../.venv/bin/python ../analyzers/vegas-analyzer-heisenberg.py test_results.h5

# Using the XYZ analyzer (produces XYZ file for visualization)
../.venv/bin/python ../analyzers/vegas-analyzer-xyz.py test_results.h5
```

The analyzers produce:
- `.mean` files: Mean magnetization and energy values
- `.pdf` files: Plots of magnetization vs. Monte Carlo steps
- `.xyz` files: Spin configurations for visualization

## File Formats

See [FORMATS.md](FORMATS.md) for detailed specifications of:
- Lattice file format
- Initial spin files
- Anisotropy files
- JSON configuration format
- HDF5 output structure

## Troubleshooting

### Common Issues

1. **HDF5 library not found**:
   Ensure `libhdf5-dev` and `libhdf5-cpp-103-1t64` packages are installed.

2. **JSON library not found**:
   Install `libjsoncpp-dev`.

3. **Python dependencies missing**:
   Run `uv sync` to create virtual environment with all required packages.

4. **Anisotropy file format error**:
   Check that anisotropy files have either 4 (uniaxial) or 7 (cubic) numbers per line.

5. **Simulation runs but produces zero magnetization**:
   Verify initial spins are normalized and temperature range is appropriate for the system.

### Debugging Tips

- Compile with debug symbols: `cmake .. -DCMAKE_BUILD_TYPE=Debug`
- Run with `valgrind` to check for memory issues
- Check simulation logs for warnings about missing files or format errors
- Verify lattice file matches expected format (see FORMATS.md)

## Validation Test

The `test_system/` directory serves as a validation test. To run the complete validation:

```bash
# Build VEGAS
cd build && make -j$(nproc)

# Run simulation
./vegas ../test_system/test_config.json

# Verify output
cd ../test_system
../.venv/bin/python verify_hdf5.py test_results.h5

# Run analyzer (optional)
../.venv/bin/python ../analyzers/vegas-analyzer-lite.py test_results.h5
```

All steps should complete without errors. The HDF5 file should contain non-zero magnetization values indicating proper simulation.

## Physical Validation Benchmarks

For comprehensive validation against exact solutions, run the benchmark suite:

```bash
# From project root
python run_all_benchmarks.py

# This runs 6 benchmarks:
# - B1: 1D Ising energy vs exact solution
# - B2: 2D Ising critical temperature vs Onsager
# - B3: Ferromagnet ground state energy
# - B4: Ergodicity test
# - B5: Sigma freeze verification
# - B6: Heisenberg high-temperature isotropy

# Results in benchmarks/VALIDATION_REPORT.md
```

## Next Steps

- **Performance Optimization**: Profile simulation with larger systems
- **Custom Spin Models**: Implement new models in `spin_model.cc`
- **Parallel Execution**: Enable OpenMP/MPI for larger parameter sweeps
- **Visualization**: Use generated XYZ files with visualization tools (VMD, OVITO)

---

*Last updated: February 2026*