cmake_minimum_required(VERSION 3.20)
project(vegas LANGUAGES CXX C)

# Sanitizer options (off by default)
option(ENABLE_ASAN "Enable AddressSanitizer (memory errors)" OFF)
option(ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer (undefined behavior)" OFF)
option(ENABLE_VECTORIZATION_REPORT "Enable vectorization optimization reports" ON)
option(ENABLE_AVX2_EXPERIMENTAL "Enable experimental AVX2 intrinsics (gather-based)" OFF)

# Apply sanitizer flags if enabled
if(ENABLE_ASAN OR ENABLE_UBSAN)
    set(SANITIZER_FLAGS "-fno-omit-frame-pointer")
    
    if(ENABLE_ASAN)
        set(SANITIZER_FLAGS "${SANITIZER_FLAGS} -fsanitize=address")
        message(STATUS "AddressSanitizer enabled")
    endif()
    
    if(ENABLE_UBSAN)
        set(SANITIZER_FLAGS "${SANITIZER_FLAGS} -fsanitize=undefined")
        message(STATUS "UndefinedBehaviorSanitizer enabled")
    endif()
    
    # Set compiler and linker flags
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SANITIZER_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${SANITIZER_FLAGS}")
    
    # ASan options with suppressions
    set(ASAN_OPTIONS "suppressions=${CMAKE_SOURCE_DIR}/sanitizer_suppressions.txt")
    set(LSAN_OPTIONS "suppressions=${CMAKE_SOURCE_DIR}/sanitizer_suppressions.txt")
endif()

# Enable testing
enable_testing()

# Add cxxopts for command line parsing (header-only)
include(FetchContent)
FetchContent_Declare(
    cxxopts
    GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
    GIT_TAG        v3.3.0
)
FetchContent_MakeAvailable(cxxopts)

# Try to find packages using MODULE mode first (works with system packages)
# Then fall back to CONFIG mode (works with Conan)
find_package(HDF5 REQUIRED COMPONENTS C CXX)

# Try to find jsoncpp using pkg-config or module mode
find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    pkg_check_modules(JSONCPP QUIET jsoncpp)
endif()

if(JSONCPP_FOUND)
    # Found via pkg-config
    add_library(jsoncpp_target INTERFACE)
    target_include_directories(jsoncpp_target INTERFACE ${JSONCPP_INCLUDE_DIRS})
    target_link_libraries(jsoncpp_target INTERFACE ${JSONCPP_LIBRARIES})
    set(JSONCPP_TARGET jsoncpp_target)
else()
    # Try to find via find_package with different names
    find_package(jsoncpp CONFIG QUIET)
    find_package(JsonCpp CONFIG QUIET)
    
    if (TARGET jsoncpp::jsoncpp)
        set(JSONCPP_TARGET jsoncpp::jsoncpp)
    elseif (TARGET JsonCpp::JsonCpp)
        set(JSONCPP_TARGET JsonCpp::JsonCpp)
    else()
        # Last resort: try to find headers and library manually
        find_path(JSONCPP_INCLUDE_DIR json/json.h)
        find_library(JSONCPP_LIBRARY jsoncpp)
        
        if(JSONCPP_INCLUDE_DIR AND JSONCPP_LIBRARY)
            add_library(jsoncpp_target INTERFACE)
            target_include_directories(jsoncpp_target INTERFACE ${JSONCPP_INCLUDE_DIR})
            target_link_libraries(jsoncpp_target INTERFACE ${JSONCPP_LIBRARY})
            set(JSONCPP_TARGET jsoncpp_target)
        else()
            message(FATAL_ERROR "JsonCpp not found. Install libjsoncpp-dev package or use Conan.")
        endif()
    endif()
endif()

# AVX2 experimental intrinsics
if(ENABLE_AVX2_EXPERIMENTAL)
    add_definitions(-DVEGAS_AVX2_EXPERIMENTAL)
    message(STATUS "Experimental AVX2 intrinsics enabled (gather-based)")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -fopenmp-simd")
if(ENABLE_VECTORIZATION_REPORT)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopt-info-vec-optimized -fopt-info-vec-missed")
    message(STATUS "Vectorization optimization reports enabled")
endif()

# Add HDF5 include directories
include_directories(${HDF5_INCLUDE_DIRS})

# Main executable
add_executable(vegas ./src/main.cc)
target_sources(
    vegas PRIVATE
    ./src/atom.cc
    ./src/lattice.cc
    ./src/reporter.cc
    ./src/system.cc
    ./src/starter.cc
    ./src/spin_model.cc
    ./src/config_parser.cc
    ./src/system_builder.cc
)
target_link_libraries(vegas PRIVATE ${JSONCPP_TARGET} ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES} cxxopts)

# Test executable (if Google Test is available)
find_package(GTest QUIET)
if(GTest_FOUND)
    add_executable(vegas_tests tests/test_basic.cc)
    target_sources(
        vegas_tests PRIVATE
        ./src/atom.cc
        ./src/lattice.cc
        ./src/reporter.cc
        ./src/system.cc
        ./src/starter.cc
        ./src/spin_model.cc
        ./src/config_parser.cc
        ./src/system_builder.cc
    )
    target_link_libraries(vegas_tests PRIVATE ${JSONCPP_TARGET} ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES} GTest::gtest GTest::gtest_main)
    add_test(NAME vegas_tests COMMAND vegas_tests)
endif()

# Simple test executable (no external dependencies)
add_executable(vegas_simple_tests tests/test_simple.cc)
target_sources(
    vegas_simple_tests PRIVATE
    ./src/atom.cc
    ./src/lattice.cc
    ./src/reporter.cc
    ./src/system.cc
    ./src/starter.cc
    ./src/spin_model.cc
    ./src/config_parser.cc
    ./src/system_builder.cc
)
target_link_libraries(vegas_simple_tests PRIVATE ${JSONCPP_TARGET} ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES})
add_test(NAME vegas_simple_tests COMMAND vegas_simple_tests)

# Integration test executable
add_executable(vegas_integration_tests tests/test_integration.cc)
target_sources(
    vegas_integration_tests PRIVATE
    ./src/atom.cc
    ./src/lattice.cc
    ./src/reporter.cc
    ./src/system.cc
    ./src/starter.cc
    ./src/spin_model.cc
    ./src/config_parser.cc
    ./src/system_builder.cc
)
target_link_libraries(vegas_integration_tests PRIVATE ${JSONCPP_TARGET} ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES})
add_test(NAME vegas_integration_tests COMMAND vegas_integration_tests)

# ConfigParser test executable
add_executable(vegas_config_parser_tests tests/test_config_parser.cc)
target_sources(
    vegas_config_parser_tests PRIVATE
    ./src/config_parser.cc
    ./src/atom.cc
    ./src/lattice.cc
    ./src/reporter.cc
    ./src/system.cc
    ./src/starter.cc
    ./src/spin_model.cc
    ./src/system_builder.cc
)
target_link_libraries(vegas_config_parser_tests PRIVATE ${JSONCPP_TARGET} ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES})
add_test(NAME vegas_config_parser_tests COMMAND vegas_config_parser_tests)

# SystemBuilder test executable
add_executable(vegas_system_builder_tests tests/test_system_builder.cc)
target_sources(
    vegas_system_builder_tests PRIVATE
    ./src/config_parser.cc
    ./src/system_builder.cc
    ./src/atom.cc
    ./src/lattice.cc
    ./src/reporter.cc
    ./src/system.cc
    ./src/starter.cc
    ./src/spin_model.cc
)
target_link_libraries(vegas_system_builder_tests PRIVATE ${JSONCPP_TARGET} ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES})
add_test(NAME vegas_system_builder_tests COMMAND vegas_system_builder_tests)

# System test executable
add_executable(vegas_system_tests tests/test_system.cc)
target_sources(
    vegas_system_tests PRIVATE
    ./src/atom.cc
    ./src/lattice.cc
    ./src/reporter.cc
    ./src/system.cc
    ./src/starter.cc
    ./src/spin_model.cc
    ./src/config_parser.cc
    ./src/system_builder.cc
)
target_link_libraries(vegas_system_tests PRIVATE ${JSONCPP_TARGET} ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES})
add_test(NAME vegas_system_tests COMMAND vegas_system_tests)

# Coloring test executable
add_executable(vegas_coloring_tests tests/test_coloring.cc)
target_sources(
    vegas_coloring_tests PRIVATE
    ./src/atom.cc
    ./src/lattice.cc
    ./src/reporter.cc
    ./src/spin_model.cc
    ./src/starter.cc
    ./src/system.cc
    ./src/config_parser.cc
    ./src/system_builder.cc
)
target_link_libraries(vegas_coloring_tests PRIVATE ${JSONCPP_TARGET} ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES})
add_test(NAME vegas_coloring_tests COMMAND vegas_coloring_tests)

# Reporter test executable
add_executable(vegas_reporter_tests tests/test_reporter.cc)
target_sources(
    vegas_reporter_tests PRIVATE
    ./src/atom.cc
    ./src/lattice.cc
    ./src/reporter.cc
    ./src/spin_model.cc
    ./src/starter.cc
    ./src/system.cc
    ./src/config_parser.cc
    ./src/system_builder.cc
)
target_link_libraries(vegas_reporter_tests PRIVATE ${JSONCPP_TARGET} ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES})
add_test(NAME vegas_reporter_tests COMMAND vegas_reporter_tests)
