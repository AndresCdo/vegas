// Simple test file to verify basic functionality
// This doesn't require Google Test installation

#include "../include/params.h"
#include <iostream>
#include <cassert>

void test_floating_point_comparison() {
    std::cout << "Testing floating point comparison... ";
    assert(fp_equal(1.0, 1.0));
    assert(fp_equal(1.0, 1.0 + EPSILON/2.0));
    assert(!fp_equal(1.0, 1.0 + EPSILON*2.0));
    assert(fp_not_equal(1.0, 1.0 + EPSILON*2.0));
    std::cout << "PASS\n";
}

void test_constants() {
    std::cout << "Testing constants... ";
    assert(DEFAULT_MCS == 5000);
    assert(fp_equal(DEFAULT_TEMP_START, 0.001));
    assert(fp_equal(DEFAULT_TEMP_FINAL, 10.0));
    assert(fp_equal(DEFAULT_TEMP_DELTA, 0.1));
    assert(fp_equal(MAX_SIGMA, 60.0));
    assert(fp_equal(MIN_SIGMA, 1e-10));
    assert(fp_equal(SIGMA_ADJUSTMENT_FACTOR, 0.5));
    assert(fp_equal(MIN_REJECTION_RATE, 1e-10));
    assert(NUM_SPIN_MODELS == 5);
    assert(fp_equal(SPIN_NORM_ROUNDING_PRECISION, 10000.0));
    std::cout << "PASS\n";
}

void test_array_operations() {
    std::cout << "Testing array operations... ";
    Array a = {1.0, 2.0, 3.0};
    Array b = {4.0, 5.0, 6.0};
    
    Array sum = a + b;
    assert(fp_equal(sum[0], 5.0));
    assert(fp_equal(sum[1], 7.0));
    assert(fp_equal(sum[2], 9.0));
    
    Array product = a * b;
    assert(fp_equal(product[0], 4.0));
    assert(fp_equal(product[1], 10.0));
    assert(fp_equal(product[2], 18.0));
    
    Real dot_product = (a * b).sum();
    assert(fp_equal(dot_product, 32.0));
    std::cout << "PASS\n";
}

void test_zero_constant() {
    std::cout << "Testing ZERO constant... ";
    assert(fp_equal(ZERO[0], 0.0));
    assert(fp_equal(ZERO[1], 0.0));
    assert(fp_equal(ZERO[2], 0.0));
    std::cout << "PASS\n";
}

int main() {
    std::cout << "Running simple tests for vegas...\n";
    
    test_floating_point_comparison();
    test_constants();
    test_array_operations();
    test_zero_constant();
    
    std::cout << "\nAll tests passed!\n";
    return 0;
}