# VEGAS Project - Agent Instructions

## Project Overview

VEGAS (VEctor General Atomistic Simulator) is a Monte Carlo simulation package for magnetic materials. The project implements atomistic simulations with multiple spin models, magnetic interactions, and HDF5 data output.

## Current State (February 2025)

### ✅ **COMPLETED IMPROVEMENTS**

#### **Critical Fixes:**
1. **Circular Dependency**: Fixed between `system.h` and `starter.h` using forward declaration
2. **Unsafe `atof()`**: Replaced all instances with `safe_stod()` with proper error handling
3. **Floating Point Comparisons**: Implemented `fp_equal()` and `fp_not_equal()` with epsilon
4. **Division by Zero**: Prevented in sigma adjustment (`system.cc:232-244`)
5. **Magic Numbers**: Defined as constants in `params.h`

#### **Build System:**
1. **CMake Configuration**: Works with system packages (libjsoncpp-dev, libhdf5-dev)
2. **HDF5 Linking**: Fixed to include both C and C++ libraries
3. **JSON Library**: Multiple detection methods (pkg-config, find_package, manual)
4. **Testing Framework**: Basic tests integrated with CMake

#### **Code Quality:**
1. **Move Semantics**: Implemented for `Lattice` and `System` classes
2. **Input Validation**: Enhanced `CHECKFILE()` with security checks
3. **Error Handling**: Improved throughout codebase
4. **Documentation**: Comprehensive README, FORMATS.md, example config

### 🚨 **KNOWN ISSUES**

#### **Build/Configuration:**
1. **Conan Not Required**: Project now uses system packages by default
2. **HDF5 Header**: Changed from `H5Include.h` to standard `hdf5.h`
3. **JSON Headers**: Using system `json/json.h`

#### **Code Quality (To Do):**
1. **`Atom` Class**: Still large (500+ lines), needs refactoring
2. **Exception Handling**: Still uses `EXIT()` function instead of exceptions
3. **RAII for HDF5**: Manual resource management in `reporter.cc`
4. **Test Coverage**: Limited, needs expansion

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
1. **Simple Tests**: `./vegas_simple_tests` - Basic functionality tests
2. **Google Tests**: If GTest found, `./vegas_tests` (not built by default)
3. **Integration**: Manual testing with JSON files

### **Test Coverage Areas:**
- Floating point comparisons (`fp_equal`, `fp_not_equal`)
- Constants and type definitions
- Array operations
- Safe string-to-double conversion

## Code Structure

### **Key Files:**
- `include/params.h` - Type definitions, constants, utility functions
- `include/system.h` - Main simulation class
- `include/atom.h` - Atom class (needs refactoring)
- `include/lattice.h` - Lattice class
- `include/reporter.h` - HDF5 output class
- `include/starter.h` - JSON configuration parser
- `include/spin_model.h` - Spin model abstraction (partial)

### **Recent Changes:**
1. **`params.h`**: Added `EPSILON`, floating point comparison functions, constants
2. **`system.cc`**: Added `safe_stod()`, fixed division by zero, updated magic numbers
3. **`starter.cc`**: Enhanced `CHECKFILE()` with security checks
4. **`CMakeLists.txt`**: Improved dependency detection and linking

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

## Next Development Phases

### **Phase 1: Testing & Validation (High Priority)**
1. Create integration tests with sample data
2. Test JSON parsing with `example_config.json`
3. Verify HDF5 output functionality
4. Install Python dependencies for analyzers (h5py, matplotlib)

### **Phase 2: Refactoring (Medium Priority)**
1. Complete `SpinModel` abstraction and refactor `Atom` class
2. Replace `EXIT()` with exception handling
3. Implement RAII wrappers for HDF5 resources
4. Add comprehensive unit tests

### **Phase 3: Features (Low Priority)**
1. Parallelization (OpenMP/MPI support)
2. Checkpoint/restart functionality
3. Performance optimization
4. GUI/Web interface

## Agent Workflow

### **When Starting Work:**
1. Read this AGENTS.md file
2. Check `git status` for uncommitted changes
3. Verify build works: `cd build && make clean && make`
4. Run tests: `./vegas_simple_tests`

### **Before Making Changes:**
1. Understand the architecture (see README.md)
2. Check for existing constants in `params.h`
3. Review similar code patterns in the codebase
4. Consider backward compatibility

### **After Making Changes:**
1. Test compilation: `cd build && make`
2. Run basic tests: `./vegas_simple_tests`
3. Test main executable: `./vegas --help`
4. Consider adding/updating tests

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
- **Documentation**: README.md, FORMATS.md
- **Examples**: `example_config.json`
- **Tests**: `tests/` directory

## Version History

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

*Last updated: February 2025*  
*Maintained by: VEGAS Development Team*