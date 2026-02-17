cmake_minimum_required(VERSION 3.20)
project(vegas LANGUAGES CXX C)

# Enable testing
enable_testing()

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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

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
)
target_link_libraries(vegas PRIVATE ${JSONCPP_TARGET} ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES})

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
)
target_link_libraries(vegas_simple_tests PRIVATE ${JSONCPP_TARGET} ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES})
add_test(NAME vegas_simple_tests COMMAND vegas_simple_tests)
