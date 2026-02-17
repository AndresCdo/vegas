#ifndef PARAMS
#define PARAMS

#include <valarray>
#include <string>
#include <cmath>

typedef double Real;
typedef std::valarray<Real> Array;
typedef unsigned int Index;
const Array ZERO = {0.0, 0.0, 0.0};
const Index AMOUNTCHUNKS = 5;
const Real EPSILON = 1e-10;

// Simulation parameters
const Index DEFAULT_MCS = 5000;
const Real DEFAULT_TEMP_START = 0.001;
const Real DEFAULT_TEMP_FINAL = 10.0;
const Real DEFAULT_TEMP_DELTA = 0.1;
const Real DEFAULT_FIELD_START = 0.001;
const Real DEFAULT_FIELD_FINAL = 10.0;
const Real DEFAULT_FIELD_DELTA = 0.1;

// Spin model parameters
const Index NUM_SPIN_MODELS = 5;
const Real SPIN_NORM_ROUNDING_PRECISION = 10000.0;

// Sigma adjustment parameters
const Real MAX_SIGMA = 60.0;
const Real MIN_SIGMA = 1e-10;
const Real SIGMA_ADJUSTMENT_FACTOR = 0.5;
const Real MIN_REJECTION_RATE = 1e-10;

template <typename T>
std::ostream & operator << (std::ostream &o, std::valarray<T> val)
{
    for (Index i = 0; i < val.size(); ++i)
    {
        o << val[i] << "\t";
    }
    return o;
}

// Floating point comparison with epsilon
inline bool fp_equal(Real a, Real b, Real epsilon = EPSILON)
{
    return std::abs(a - b) < epsilon;
}

inline bool fp_not_equal(Real a, Real b, Real epsilon = EPSILON)
{
    return std::abs(a - b) >= epsilon;
}

#endif