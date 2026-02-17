#include <gtest/gtest.h>
#include "../include/params.h"

// Test the floating point comparison functions
TEST(FloatingPointTest, BasicEquality) {
    EXPECT_TRUE(fp_equal(1.0, 1.0));
    EXPECT_TRUE(fp_equal(1.0, 1.0 + EPSILON/2.0));
    EXPECT_FALSE(fp_equal(1.0, 1.0 + EPSILON*2.0));
    EXPECT_TRUE(fp_not_equal(1.0, 1.0 + EPSILON*2.0));
}

// Test the safe_stod function (if we make it available in params.h)
// For now, test basic type definitions
TEST(TypeTest, BasicTypes) {
    Real test_real = 1.5;
    EXPECT_EQ(typeid(test_real), typeid(double));
    
    Array test_array = {1.0, 2.0, 3.0};
    EXPECT_EQ(test_array.size(), 3);
    EXPECT_DOUBLE_EQ(test_array[0], 1.0);
    EXPECT_DOUBLE_EQ(test_array[1], 2.0);
    EXPECT_DOUBLE_EQ(test_array[2], 3.0);
}

// Test constants are defined correctly
TEST(ConstantsTest, MagicNumbers) {
    EXPECT_EQ(DEFAULT_MCS, 5000);
    EXPECT_DOUBLE_EQ(DEFAULT_TEMP_START, 0.001);
    EXPECT_DOUBLE_EQ(DEFAULT_TEMP_FINAL, 10.0);
    EXPECT_DOUBLE_EQ(DEFAULT_TEMP_DELTA, 0.1);
    EXPECT_DOUBLE_EQ(MAX_SIGMA, 60.0);
    EXPECT_DOUBLE_EQ(MIN_SIGMA, 1e-10);
    EXPECT_DOUBLE_EQ(SIGMA_ADJUSTMENT_FACTOR, 0.5);
    EXPECT_DOUBLE_EQ(MIN_REJECTION_RATE, 1e-10);
    EXPECT_EQ(NUM_SPIN_MODELS, 5);
    EXPECT_DOUBLE_EQ(SPIN_NORM_ROUNDING_PRECISION, 10000.0);
}

// Test array operations
TEST(ArrayTest, BasicOperations) {
    Array a = {1.0, 2.0, 3.0};
    Array b = {4.0, 5.0, 6.0};
    
    Array sum = a + b;
    EXPECT_DOUBLE_EQ(sum[0], 5.0);
    EXPECT_DOUBLE_EQ(sum[1], 7.0);
    EXPECT_DOUBLE_EQ(sum[2], 9.0);
    
    Array product = a * b;
    EXPECT_DOUBLE_EQ(product[0], 4.0);
    EXPECT_DOUBLE_EQ(product[1], 10.0);
    EXPECT_DOUBLE_EQ(product[2], 18.0);
    
    Real dot_product = (a * b).sum();
    EXPECT_DOUBLE_EQ(dot_product, 32.0); // 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32
}

// Test ZERO constant
TEST(ConstantsTest, ZeroArray) {
    EXPECT_DOUBLE_EQ(ZERO[0], 0.0);
    EXPECT_DOUBLE_EQ(ZERO[1], 0.0);
    EXPECT_DOUBLE_EQ(ZERO[2], 0.0);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}