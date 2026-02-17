#!/bin/bash
# VEGAS Validation Script
# Runs a complete end-to-end validation of the simulation pipeline

set -e  # Exit on any error

echo "=== VEGAS Validation Script ==="
echo ""

# Check prerequisites
if [ ! -d "build" ]; then
    echo "ERROR: Build directory not found. Please build VEGAS first."
    echo "Run: mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make"
    exit 1
fi

if [ ! -f "build/vegas" ]; then
    echo "ERROR: vegas executable not found in build directory."
    echo "Please ensure the build completed successfully."
    exit 1
fi

if [ ! -d "test_system" ]; then
    echo "ERROR: test_system directory not found."
    echo "This directory should contain the validation test files."
    exit 1
fi

if [ ! -f ".venv/bin/python" ]; then
    echo "ERROR: Python virtual environment not found."
    echo "Please create it with: uv sync"
    exit 1
fi

# Step 1: Run simulation
echo "1. Running simulation..."
cd test_system
../build/vegas test_config.json
echo "   Simulation completed. Output: test_results.h5"

# Step 2: Verify HDF5 structure
echo ""
echo "2. Verifying HDF5 output structure..."
../.venv/bin/python verify_hdf5.py test_results.h5
if [ $? -ne 0 ]; then
    echo "   ERROR: HDF5 verification failed."
    exit 1
fi
echo "   HDF5 structure validated."

# Step 3: Run analyzer (optional)
echo ""
echo "3. Running analyzer..."
if [ -f "../analyzers/vegas-analyzer-lite.py" ]; then
    ../.venv/bin/python ../analyzers/vegas-analyzer-lite.py test_results.h5
    echo "   Analyzer completed. Output: test_results.mean and test_results.pdf"
else
    echo "   WARNING: Analyzer script not found, skipping."
fi

# Step 4: Check for expected outputs
echo ""
echo "4. Checking output files..."
required_files=("test_results.h5")
for file in "${required_files[@]}"; do
    if [ -f "$file" ]; then
        echo "   ✓ $file exists"
    else
        echo "   ✗ $file missing"
        exit 1
    fi
done

# Optional: Check file sizes
echo ""
echo "5. Output file details:"
ls -lh test_results.h5
if [ -f "test_results.mean" ]; then
    ls -lh test_results.mean
fi
if [ -f "test_results.pdf" ]; then
    ls -lh test_results.pdf
fi

echo ""
echo "=== Validation Successful ==="
echo "All validation steps completed successfully."
echo ""
echo "The VEGAS simulation pipeline is working correctly."
echo "You can now run your own simulations using the documented workflow."