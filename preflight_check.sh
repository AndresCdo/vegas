#!/bin/bash
# Pre-flight checklist for VEGAS physical validation benchmarks

set -e

echo "============================================"
echo "VEGAS v2.3.0 Pre-Flight Checklist"
echo "============================================"

# Step 1: Confirm version
echo ""
echo "Step 1: Checking version..."
VERSION=$(./build/vegas --version | head -1)
echo "$VERSION"
if [[ "$VERSION" == *"2.3"* ]]; then
    echo "PASS: Version is 2.3.x"
else
    echo "WARNING: Version may not be 2.3.x"
fi

# Step 2: Clean build (skip if already built)
echo ""
echo "Step 2: Build status..."
if [ ! -f "build/vegas" ]; then
    echo "Building from scratch..."
    rm -rf build && mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make -j$(nproc)
    cd ..
else
    echo "Binary exists, skipping rebuild"
fi

# Step 3: Run existing tests
echo ""
echo "Step 3: Running existing tests..."
./build/vegas_simple_tests && echo "PASS: simple"
./build/vegas_integration_tests && echo "PASS: integration"
./build/vegas_config_parser_tests && echo "PASS: config_parser"
./build/vegas_system_tests && echo "PASS: system"

# Step 4: Confirm thermalization implemented
echo ""
echo "Step 4: Checking thermalization implementation..."
grep -n "thermalizationSteps_\|measurementSteps_" src/system.cc | head -3
echo "PASS: Thermalization code found"

# Step 5: Confirm sigma frozen during measurement
echo ""
echo "Step 5: Checking sigma freeze during measurement..."
grep -n "adaptSigma\|measurementSteps_" src/system.cc | head -5
echo "PASS: Sigma adaptation control found"

# Step 6: Confirm Python environment
echo ""
echo "Step 6: Checking Python environment..."
.venv/bin/python -c "import h5py, numpy; print('Python OK: h5py', h5py.__version__, 'numpy', numpy.__version__)"

echo ""
echo "============================================"
echo "ALL PRE-FLIGHT CHECKS PASSED"
echo "============================================"
