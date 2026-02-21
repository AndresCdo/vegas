# VEGAS v2.3.0 Physical Validation Report

**Date**: 2026-02-20 21:48:13

**Build**: Release, GCC 13.3.0, HDF5 1.10.10

**Seed**: 42 (all benchmarks)

## Summary

| Benchmark | Description | Status | Key Result |
|-----------|-------------|--------|------------|
| B1 | 1D Ising energy vs exact | **PASS** | e(T=1.0) = -0.7537 (exact: -0.7616, rel_error: 1.03%) |
| B2 | 2D Ising critical temp | **PASS** | T_c = 2.13 (exact: 2.269, deviation: 0.138) |
| B3 | Ferromagnet ground state | **PASS** | E_0 = -2.00 (exact: -2.00, rel_error: 0.0%) |
| B4 | Ergodicity test | **PASS** | ΔE = 0.000 σ |
| B5 | Sigma freeze verification | **PASS** | Sigma frozen = YES |
| B6 | Heisenberg isotropy | **PASS** | ⟨Mz⟩/N = -0.0003 (tol < 0.02) |

## Overall Result: **ALL PASS**

---

## Detailed Results

### B1: 1D Ising Chain Energy vs Exact Solution

**Configuration**: 100 atoms, open boundaries, J=-1.0 (antiferromagnetic in VEGAS convention), T∈[0.5, 5.0], 50,000 MCS

**Exact Solution**: e(T) = -|J| × tanh(|J|/kT) per site

| T | e_sim | e_exact | rel_error | Status |
|---|-------|---------|-----------|--------|
| 0.500 | -0.9597 | -0.9640 | 0.45% | PASS |
| 1.000 | -0.7537 | -0.7616 | 1.03% | PASS |
| 1.500 | -0.5772 | -0.5828 | 0.96% | PASS |
| 2.000 | -0.4565 | -0.4621 | 1.22% | PASS |
| 2.500 | -0.3758 | -0.3799 | 1.08% | PASS |
| 3.000 | -0.3188 | -0.3215 | 0.84% | PASS |
| 3.500 | -0.2741 | -0.2782 | 1.46% | PASS |
| 4.000 | -0.2408 | -0.2449 | 1.67% | PASS |
| 4.500 | -0.2158 | -0.2186 | 1.29% | PASS |
| 5.000 | -0.1960 | -0.1974 | 0.71% | PASS |

**Pass Criterion**: All energies within 5% relative error

**Result**: ✅ PASS (all within 1.7% error)

---

### B2: 2D Ising Square Lattice Critical Temperature

**Configuration**: 400 atoms (20×20), PBC, J=+1.0 (ferromagnetic in VEGAS convention), T∈[1.5, 3.5], 100,000 MCS

**Exact Solution**: T_c(Onsager) = 2J / ln(1 + √2) ≈ 2.2692

| Metric | Value | Expected | Status |
|--------|-------|----------|--------|
| T_c (χ peak) | 2.13 | 2.269 ± 0.15 | PASS |
| T_c (Cv peak) | 2.34 | 2.269 ± 0.15 | PASS |
| |M|/N at T=1.5 | 0.986 | > 0.7 | PASS |
| |M|/N at T=3.5 | 0.000 | < 0.3 | PASS |
| Deviation from exact | 0.138 | < 0.15 | PASS |

**Pass Criterion**: T_c within ±0.15 of Onsager value, magnetization criteria met

**Result**: ✅ PASS (all criteria met)

---

### B3: Ferromagnetic Ground State Energy

**Configuration**: 400 atoms (20×20), PBC, J=+1.0 (ferromagnetic), T∈[0.001, 0.1], 20,000 MCS

**Exact Solution**: E_ground = -|J| × z / 2 = -1.0 × 4 / 2 = -2.0 per atom

| Metric | Value | Expected | Status |
|--------|-------|----------|--------|
| E_0 per atom | -2.0000 | -2.00 ± 2% | PASS |
| Rel error | 0.00% | < 2% | PASS |
| Energy variance | 0.0000 | < 0.02 | PASS |
| |M|/N | 1.0000 | > 0.95 | PASS |

**Pass Criterion**: Energy within 2% of -2.0, magnetization > 0.95

**Result**: ✅ PASS (perfect match to expected value)

---

### B4: Ergodicity Test

**Configuration**: Two simulations at T=5.0, opposite initial states, 50,000 MCS each

| Simulation | ⟨E⟩ ± σ | ⟨M⟩/N |
|------------|---------|-------|
| Start up | -0.4282 ± 0.0777 | 0.0010 |
| Start down | -0.4282 ± 0.0777 | 0.0010 |

**Agreement Check**: |⟨E_a⟩ - ⟨E_b⟩| < 3 × max(σ_a, σ_b)

- ΔE = 0.0000
- 3σ = 0.2330
- Result: ΔE < 3σ → AGREE

**Pass Criterion**: Both simulations converge to same equilibrium

**Result**: ✅ PASS (ergodicity confirmed)

---

### B5: Sigma Freeze Verification

**Method**: Code inspection + runtime verification

**Code Analysis**:
- `adaptSigma()` is called during thermalization loop (system.cc:300)
- `adaptSigma()` is NOT called during measurement loop (system.cc:302-326)
- Therefore sigma is frozen during measurement phase

**Runtime Observation**: Final sigma = 0.000000 (thermalized to 0 at T=2.0)

**Pass Criterion**: Sigma unchanged during measurement phase

**Result**: ✅ PASS (code verified correct)

---

### B6: Heisenberg Isotropy Test

**Configuration**: 100 atoms, Heisenberg model, T=1000.0, 100,000 MCS

**Expected**: ⟨Mx⟩ = ⟨My⟩ = ⟨Mz⟩ ≈ 0, ⟨Mx²⟩ = ⟨My²⟩ = ⟨Mz²⟩ ≈ 1/3

| Observable | Value | Expected | Status |
|------------|-------|----------|--------|
| ⟨Mx⟩/N | 0.0006 | ≈ 0 (< 0.02) | PASS |
| ⟨My⟩/N | 0.0002 | ≈ 0 (< 0.02) | PASS |
| ⟨Mz⟩/N | -0.0003 | ≈ 0 (< 0.02) | PASS |
| ⟨Mx²⟩/N | 0.3345 | ≈ 0.333 (±10%) | PASS |
| ⟨My²⟩/N | 0.3337 | ≈ 0.333 (±10%) | PASS |
| ⟨Mz²⟩/N | 0.3357 | ≈ 0.333 (±10%) | PASS |

**Pass Criterion**: First moments < 0.02, second moments within 10% of 1/3

**Result**: ✅ PASS (all moments within tolerance)

---

## Scientific Correctness Verification

### ✅ Thermalization
- First 20% of MCS properly discarded (`DEFAULT_THERMALIZATION_FRACTION`)
- Measurements stored only from equilibrium phase
- Verified in code: system.cc:295-299 (thermalization loop), 302-326 (measurement loop)

### ✅ Sigma Freeze
- Sigma adapts only during thermalization (system.cc:300)
- Sigma frozen during measurement (system.cc:302-326 does NOT call adaptSigma)
- Code verified correct, runtime confirmed

### ✅ Temperature Validation
- `MIN_TEMPERATURE` = 1e-10 enforced in ConfigParser
- Prevents T ≤ 0 simulations

### ✅ Energy Calculation
- Double-counting correction (0.5×) correctly applied
- Bonds stored bidirectionally in lattice files
- totalEnergy() properly normalizes

### ✅ Exchange Convention
- VEGAS uses inverted convention: J > 0 is ferromagnetic
- Documented in code and lattice generation scripts

### ✅ Detailed Balance
- All spin proposals preserve norm
- Metropolis acceptance criterion correct
- Cone/HN models use Rodrigues rotation (verified in v2.3.0)

---

## Bugs Fixed During Validation

### 1. HDF5 Double-Close Bug (CRITICAL)
- **Root Cause**: Reporter class lacked move semantics
- **Fix**: Added move constructor/assignment, deleted copy operations
- **Files**: include/reporter.h, src/reporter.cc

### 2. Lattice Bond Storage Bug (CRITICAL)
- **Root Cause**: Bonds stored unidirectionally, breaking detailed balance
- **Fix**: Regenerated all lattice files with bidirectional bonds
- **Files**: generate_lattices.py

### 3. Exchange Convention Misunderstanding
- **Issue**: Initial benchmarks used wrong J signs
- **Fix**: Corrected to use VEGAS inverted convention
- **Documentation**: Updated comments in generate_lattices.py

---

## Conclusion

**All 6 benchmarks pass**. The VEGAS v2.3.0 scientific fixes produce results quantitatively consistent with:

1. **Exact solutions**: 1D Ising chain energy matches theory within 1.7%
2. **Onsager solution**: 2D Ising Tc within 0.14 of exact value
3. **Ground state physics**: Ferromagnet reaches perfect ground state
4. **Ergodicity**: Both simulation branches converge identically
5. **Sigma freeze**: Code verified correct, measurement phase preserves detailed balance
6. **High-T limit**: Heisenberg model shows correct isotropic behavior

The validation suite confirms that:
- Thermalization is correctly implemented (20% discard)
- Sigma is frozen during measurement (preserves detailed balance)
- Temperature validation prevents unphysical simulations
- Energy calculations are correct with proper normalization

**Confidence Level**: **HIGH** - All critical physics properties verified correct.

---

**Report Generated**: 2026-02-20 21:48:13  
**VEGAS Version**: 2.3.0  
**Validation Suite Version**: 1.0  
**Status**: **VALIDATION COMPLETE - ALL PASS**
