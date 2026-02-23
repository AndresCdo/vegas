#ifndef VEGAS_RANDOM_H
#define VEGAS_RANDOM_H

#include "params.h"

#include <cstdint>
#include <random>
#include <array>
#include <cmath>

namespace vegas {

class SplitMix64 {
public:
    using result_type = uint64_t;
    
    static constexpr uint64_t MIN = 0;
    static constexpr uint64_t MAX = UINT64_MAX;
    
    explicit SplitMix64(uint64_t seed) : state_(seed) {}
    
    uint64_t operator()() {
        uint64_t z = (state_ += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
    
    static constexpr uint64_t min() { return MIN; }
    static constexpr uint64_t max() { return MAX; }
    
private:
    uint64_t state_;
};

class Xoshiro256StarStar {
public:
    using result_type = uint64_t;
    
    static constexpr uint64_t MIN = 0;
    static constexpr uint64_t MAX = UINT64_MAX;
    
    explicit Xoshiro256StarStar(uint64_t seed) {
        SplitMix64 splitmix(seed);
        for (auto& s : state_) {
            s = splitmix();
        }
    }
    
    explicit Xoshiro256StarStar(const std::array<uint64_t, 4>& seeds) : state_(seeds) {}
    
    uint64_t operator()() {
        uint64_t* s = state_.data();
        uint64_t const t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = (s[3] << 45) | (s[3] >> 19);
        return s[1] * 5;
    }
    
    static constexpr uint64_t min() { return MIN; }
    static constexpr uint64_t max() { return MAX; }
    
    void jump() {
        uint64_t const JUMP[] = { 0x180ec6d33cfd0abaULL, 0xd5a61266f0c9392cULL,
                                   0xa958b8e01c8fd7c3ULL, 0x76a15a1eee3b49c5ULL };
        uint64_t s0 = 0, s1 = 0, s2 = 0, s3 = 0;
        for (uint64_t j : JUMP)
            for (int b = 0; b < 64; ++b) {
                if (j & (1ULL << b)) {
                    s0 ^= state_[0];
                    s1 ^= state_[1];
                    s2 ^= state_[2];
                    s3 ^= state_[3];
                }
                (*this)();
            }
        state_[0] = s0;
        state_[1] = s1;
        state_[2] = s2;
        state_[3] = s3;
    }
    
    std::array<uint64_t, 4> getState() const {
        return state_;
    }
    
    void setState(const std::array<uint64_t, 4>& newState) {
        state_ = newState;
    }
    
private:
    std::array<uint64_t, 4> state_;
};

}

#endif
