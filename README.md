# norm-functions

Normalization utilities for PET data analysis based on ROOT and CASToR.

Computes per-LOR normalization factors from GATE Monte Carlo simulation data and writes them in CASToR binary format ready for use in reconstruction.

## Requirements

- ROOT >= 6.26
- C++17
- CMake >= 3.15
- OpenMP (for parallelization)

## Build

```bash
mkdir build
cd build
cmake ..
make -j
```

To also build and run the regression tests:

```bash
ctest
```

## Workflow overview

```
GATE simulation
  ├── solidCyl source  →  *solidCyl*.root   (block correction + fan-sum efficiency)
  └── annular source   →  *annular*.root    (transaxial + interference correction)
          │
          ▼
  own-normFactors  -x scanner.xml  -i 'data/*.root'  -o output
          │
          ├── output_CB_df.Cdf / .Cdh          ← full normalization (use this for reconstruction)
          ├── output_effBfGAfGTrAf_df.Cdf/.Cdh ← without interference factor
          ├── output_CBsqrBfnI_df.Cdf/.Cdh     ← squared block, no interference
          ├── output_block_geom.csv            ← per-ring-pair diagnostic
          ├── output_transaxial.csv            ← per-radialID diagnostic
          ├── output_interference.csv          ← per-(radialID,trAID) diagnostic
          └── output_fan_sum_efficiency.csv    ← per-crystal efficiency diagnostic
```

**Input file routing** is automatic: any file whose path contains the string `solidCyl` is treated as a solid cylinder source; any file containing `annular` is treated as an annular source. Both types can be glob-expanded and multiple files of each type accumulate statistics before computing the final factors.

## Usage

```bash
./own-normFactors [OPTIONS]
```

### Required arguments

| Argument | Description |
|----------|-------------|
| `-x, --xml <path>` | Scanner configuration XML file (see `scanners/`) |
| `-i, --input <pattern>` | Input ROOT file path or glob pattern (quote wildcards) |
| `-o, --outputFile <name>` | Output file base name (without extension) |

### Optional arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `-d, --outputDir <path>` | `.` (current dir) | Directory for all output files |
| `-j, --threads <N>` | all cores | Number of OpenMP threads |
| `-as, --axial-sigma <σ>` | `0` (disabled) | Gaussian smoothing sigma for axial (`ringComponentVector`) |
| `-ts, --transaxial-sigma <σ>` | `0` (disabled) | Gaussian smoothing sigma for transaxial (`radialComponentVector`) |
| `-h, --help` | | Print usage and exit |

### Examples

```bash
# Minimal: single file, no smoothing
./own-normFactors -x scanners/32x16x2_4rings.xml \
                  -i data/sim_solidCyl.root \
                  -o norm_run1

# Multiple files via glob (must quote the pattern)
./own-normFactors -x scanners/32x16x2_4rings.xml \
                  -i 'data/sim_*.root' \
                  -o norm_allruns \
                  -d /results/norm/

# With recommended smoothing and 8 threads
./own-normFactors -x scanners/32x16x2_4rings.xml \
                  -i 'data/*.root' \
                  -o norm_smooth \
                  -as 0.6 -ts 1.0 \
                  -j 8

# Disable all smoothing explicitly
./own-normFactors -x scanners/CM2L_1ring.xml \
                  -i data.root \
                  -o norm_raw \
                  -as 0 -ts 0
```

### Smoothing guide

Both smoothing options are disabled by default (`sigma=0`). Enable them if artifacts are observed in reconstructed images.

**Axial smoothing** (`-as`) acts on `ringComponentVector` (block correction) *before* the geometric matrix is built, so it propagates cleanly into all subsequent components.

| `-as` value | Ripple reduction | Notes |
|-------------|-----------------|-------|
| `0` | none | Default; use when axial statistics are high |
| `0.6` | ~3× | Recommended starting point (3% → ~1% ripple) |
| `0.8` | ~2× | Moderate |
| `1.0` | ~1.5× | Mild |

**Transaxial smoothing** (`-ts`) acts on `radialComponentVector` (transaxial geometric correction) after annular processing.

| `-ts` value | Notes |
|------------|-------|
| `0` | Default; last bin (r=0 LORs) replaced with second-to-last if zero |
| `1.0` | Recommended when central LOR artifacts are visible |
| `1.5` | Stronger; use for very low-statistics annular runs |

---

## Scanner XML format

Scanner geometry is described in an XML file. See `scanners/` for working examples.

```xml
<?xml version="1.0" encoding="UTF-8"?>
<scanner name="MyScannerName">
    <geometry>
        <!-- Angular repeater: number of rsectors around the ring -->
        <nRsectorsAngPos>32</nRsectorsAngPos>
        <!-- Axial repeater: number of rsector rings along z (1 for most scanners) -->
        <nRsectorsAxial>1</nRsectorsAxial>

        <!-- Module repeater (within each rsector) -->
        <nModulesTransaxial>1</nModulesTransaxial>
        <nModulesAxial>4</nModulesAxial>

        <!-- Submodule repeater (within each module) -->
        <nSubmodulesTransaxial>1</nSubmodulesTransaxial>
        <nSubmodulesAxial>32</nSubmodulesAxial>

        <!-- Crystal repeater (within each submodule) -->
        <nCrystalsTransaxial>16</nCrystalsTransaxial>
        <nCrystalsAxial>1</nCrystalsAxial>

        <!-- DOI layers -->
        <nLayers>2</nLayers>
        <nLayersRptTransaxial>1</nLayersRptTransaxial>
        <nLayersRptAxial>1</nLayersRptAxial>

        <!-- Set true if transaxial detector order is reversed in GATE -->
        <invertDetOrder>false</invertDetOrder>

        <!-- Rsector ID ordering: 0 = transaxial-first (default), 1 = axial-first (cubic array) -->
        <rsectorIdOrder>0</rsectorIdOrder>
    </geometry>
    <physical>
        <crystalDepth>10.0</crystalDepth>   <!-- mm, used for DOI offset -->
        <axialSize>59.0</axialSize>          <!-- mm, submodule pitch along z -->
        <transAxialSize>59.0</transAxialSize><!-- mm, crystal pitch transaxially -->
        <detectorRadius>321.3</detectorRadius><!-- mm, inner face of crystal -->
    </physical>
</scanner>
```

### Provided scanner configurations

| File | Name | Rings | Crystals/ring | Layers |
|------|------|-------|---------------|--------|
| `CM2L_1ring.xml` | `CM2L_1ring_system` | 32 | 1024 | 2 |
| `16x16x2_1ring.xml` | `16x16x2_1ring_system` | 16 | 512 | 2 |
| `16x16x2_4rings.xml` | `16x16x2_4rings_system` | 64 | 512 | 2 |
| `32x16x2_4rings.xml` | `32x16x2_4rings_system` | 128 | 512 | 2 |

### Cubic array scanners

For scanners where rsectors, modules, submodules or crystals have *both* axial and transaxial counts (e.g. `nModulesTransaxial > 1`), set `rsectorIdOrder` appropriately:

- `0` — GATE repeats transaxially first: `rsectorID = rsectorAxl * nAngPos + rsectorTrs`
- `1` — GATE repeats axially first: `rsectorID = rsectorTrs * nAxial + rsectorAxl`

The code decomposes every combined ID into its axial and transaxial parts before computing ring and transaxial indices. See `NORMALISATION_NOTES.md` for the full index formulas.

---

## Output files

### Normalization matrices (CASToR binary format)

Three normalization files are produced, each a pair `.Cdh` (header) + `.Cdf` (data):

| File suffix | Formula | When to use |
|-------------|---------|-------------|
| `_CB` | `eff × B × G_ax × G_trs × I` | **Standard reconstruction** — full component-based correction |
| `_effBfGAfGTrAf` | `eff × B × G_ax × G_trs` | Reconstruction without interference correction |
| `_CBsqrBfnI` | `eff × B² × G_ax × G_trs` | Squared block factor, no interference — for specific studies |

Where:
- `eff` = intrinsic crystal efficiency (fan-sum algorithm, Pepin et al. 2011)
- `B` = block correction (per-ring intra-block coincidence counts)
- `G_ax` = geometric axial correction (per-ring-pair, from solid cylinder)
- `G_trs` = transaxial geometric correction (per-radialID, from annular source)
- `I` = interference / intra-block correction (per-(radialID, trAID))

To use in CASToR reconstruction, pass the `.Cdh` file as the normalization input:

```bash
castor-recon ... -norm output_CB_df.Cdh
```

### Diagnostic CSV files

These are small, human-readable files useful for inspecting individual normalization components.

| File | Columns | Description |
|------|---------|-------------|
| `_block_geom.csv` | `ring1, ring2, blockCorrection, geomAxCorrection` | Per-ring-pair corrections; useful to check axial uniformity |
| `_transaxial.csv` | `radialID, transaxialGeomNormFactor` | Transaxial correction vs. radial distance from scanner axis |
| `_interference.csv` | `radialID, trAID, interferenceTraFactor` | Intra-block interference correction per (radial, transaxial-in-block) bin |
| `_fan_sum_efficiency.csv` | `ringID, transaxialID, fanCount, ringAverage, efficiencyFactor` | Per-crystal fan counts and efficiency; flag crystals with low statistics |

---

## Re-running with different parameters

The computation has two stages with different costs:

1. **File processing** (slow — reads all ROOT events): builds `ringComponentVector`, `ringsComponentMatrix`, `radialComponentVector`, `blockTrAComponentMatrix`, `fanSumCounter`
2. **LOR loop** (fast, fully parallelised): iterates over all CastorID pairs and applies the components

Currently both stages run together. To iterate quickly on smoothing parameters or study different normalization variants, re-run with different `-as`/`-ts` values — the file processing dominates runtime so the same ROOT files should be kept accessible.

Smoothing is applied between the two passes of the solid cylinder processing (after block counts, before the geometric matrix), so **changing `-as` does affect both `ringComponentVector` and `ringsComponentMatrix`**. It is not a purely post-processing step.

---

## Development Log

### 2026-01-12: Geometric Mean Normalization

**Problem:** The normalization components were normalized to have mean=1 using arithmetic mean, but this caused overcorrections in the reconstructed images due to skewed distributions in the count data.

**Solution:** Implemented geometric mean normalization to reduce sensitivity to outliers and skewed distributions.

**Changes made:**

1. **Added geometric mean functions** (`src/normFunctions.cxx`):
   - `geometricMeanVector()` - computes geometric mean in log-space for vectors
   - `geometricMeanMatrix()` - computes geometric mean in log-space for matrices
   - Both functions handle non-finite and non-positive values gracefully

2. **Updated mean calculations** to use geometric mean instead of arithmetic mean for:
   - `meanRingComponentVector` (block correction)
   - `meanRingsComponentMatrix` (geometric axial correction)
   - `meanRadialComponentVector` (transaxial correction)
   - `meanBlockTrAComponentMatrix` (interference correction)

3. **Removed redundant second-pass normalizations**:
   - Removed `mean_block_global`, `mean_geomAx_global`, `mean_interf_global` computations
   - Removed the pre-computed `norm_transaxial` vector
   - Components are now used directly in the final normalization formula

4. **Simplified final normalization** to:
   ```cpp
   CBasedNF = effNormFactor * blockCorrection * geomAxCorrection * transaxialGeomNormFactor * interferenceTraFactor;
   ```

**Rationale:** Geometric mean is more appropriate for multiplicative correction factors because:
- It's less sensitive to extreme values (outliers)
- For skewed distributions, it provides a more representative central tendency
- It preserves the multiplicative relationships between factors

### 2026-01-14: 3D Fan-Sum Algorithm for Intrinsic Detector Efficiency

**Reference:** Pepin et al. 2011 IEEE NSS - "Normalization of Monte Carlo PET data using GATE" (Equation 4)

**Problem:** The previous efficiency calculation used separate per-component counts (submodule, crystal, layer), which doesn't properly capture the spatially-correlated detection efficiency variations described in the literature.

**Solution:** Implemented the 3D fan-sum algorithm as described in the paper:

```
ε_ui = (1/L × Σ_i' Σ_v Σ_j t_ui'vj) / (Σ_v Σ_j t_uivj)
```

Where:
- `u` = ring index (axial position)
- `i` = transaxial crystal position within ring
- `L` = number of transaxial crystals per ring
- `M` = number of rings
- `t_uivj` = coincidences for LOR between crystal (u,i) and crystal (v,j)

The **denominator** is the "fan" of crystal (u,i): sum of all coincidences in LORs containing that crystal.
The **numerator** is the average fan across all crystals in ring u (normalization factor).

**Implementation Details:**

1. **Ring definition:** see *2026-04-10* entry below for the full general formula
2. **Transaxial position:** Excludes layers; includes rsector_transaxial, module_transaxial, submodule_transaxial, crystal_transaxial
3. **Efficiency factor:** For a LOR between crystals (u1,i1) and (u2,i2):
   ```cpp
   effFactor = (ringAvg_u1 × ringAvg_u2) / (fanCount_u1i1 × fanCount_u2i2)
   ```

**Changes made:**

1. **New structure `FanSumCounter`** (`include/normFunctions.h`):
   - 2D array `fanCounts[ring][transaxial]` storing coincidences per crystal
   - Methods: `addCount()`, `getFanCount()`, `getRingAverageFan()`, `getEfficiencyFactor()`
   - Diagnostic methods: `getMinFanCount()`, `getMaxFanCount()`, `getGlobalAverageFan()`

2. **Modified `processFile()`** (`src/normFunctions.cxx`):
   - During **solidCyl** (cylinder) processing, accumulates fan-sum counts for both crystals in each coincidence
   - Computes transaxial ID excluding layers

3. **Modified `computeNormalizationFactors()`**:
   - Creates `FanSumCounter` with dimensions: `nRings = nRsectorsAxial * nModulesAxial * nSubmodulesAxial * nCrystalsAxial`, `nTransaxial = nRsectorsAngPos * nModulesTransaxial * nSubmodulesTransaxial * nCrystalsTransaxial`
   - Replaced old efficiency calculation with fan-sum based efficiency
   - Added diagnostic output after cylinder scan processing

4. **New CSV output:** `_fan_sum_efficiency.csv` containing per-crystal fan counts, ring averages, and efficiency factors

**Old vs New Efficiency Calculation:**

| Old (per-component) | New (fan-sum) |
|---------------------|---------------|
| `eff = (max²/count1×count2)` for submodule × crystal × layer separately | `eff = (ringAvg1 × ringAvg2) / (fan1 × fan2)` per unique crystal |
| Components counted independently | Crystals identified by (ring, transaxial) position |
| No spatial correlation | Proper fan-based normalization |

**Diagnostic Output:**
```
=== 3D Fan-Sum Counter Statistics ===
  Min fan count: <value>
  Max fan count: <value>
  Global average fan: <value>
  Poissonian error (min): <value>
```

### 2026-01-19: Gaussian Smoothing for Normalization Components

**Problem:** Two artifacts were observed in the normalization:
1. **Axial sawtooth rippling:** Period-2 oscillations (high-low-high-low) in the ring component vector, causing ~3% variations in the reconstructed image
2. **Central LOR artifacts:** Exploding normalization factors at r=0 in the transaxial (radial) component due to low statistics

**Solution:** Implemented Gaussian kernel smoothing for both axial and transaxial normalization components with configurable sigma parameters via command-line options.

**Mathematical Background:**

For a period-2 alternating pattern, a 3-point Gaussian kernel with weights `[w, 1-2w, w]` reduces amplitude by factor `(1 - 4w)`:
- `sigma=0.6` → weights `[0.17, 0.66, 0.17]` → ~3x reduction (3% → 1%)
- `sigma=1.0` → weights `[0.27, 0.46, 0.27]` → ~1.5x reduction

**Implementation Details:**

1. **New function `gaussianSmoothVector()`** (`src/normFunctions.cxx`):
   - Pre-computes normalized Gaussian kernel weights
   - Handles boundaries by only including valid indices
   - Skips non-finite and non-positive values
   - Returns smoothed vector

2. **New command-line options**:
   - `-as, --axial-sigma <value>`: Sigma for axial smoothing (default: 0)
   - `-ts, --transaxial-sigma <value>`: Sigma for transaxial smoothing (default: 0)
   - Set to 0 to disable smoothing

3. **Application points**:
   - **Axial:** Applied to `ringComponentVector` after solidCyl Pass 1, before Pass 2 (geometric matrix build)
   - **Transaxial:** Applied to `radialComponentVector` after annular processing, before geometric mean computation

**Kernel Parameters:**

| Component | Default Sigma | Kernel Half-Width | Kernel Size |
|-----------|---------------|-------------------|-------------|
| Axial | 0.6 | 1 | 3-point |
| Transaxial | 1.0 | 2 | 5-point |

**Quick Reference for Axial Smoothing:**

| Sigma | Ripple Reduction | Use Case |
|-------|------------------|----------|
| 0 | None | Disable smoothing |
| 0.45 | ~6x | Very aggressive |
| 0.6 | ~3x | Default (3% → 1%) |
| 0.8 | ~2x | Moderate (3% → 1.5%) |
| 1.0 | ~1.5x | Mild smoothing |

**Files Modified:**
- `include/normFunctions.h`: Added `gaussianSmoothVector()` declaration and updated `computeNormalizationFactors()` signature
- `src/normFunctions.cxx`: Added smoothing implementation and application
- `src/own-normFactors.cxx`: Added `-as` and `-ts` command-line options

### 2026-04-10: Cubic Array Indexing Generalisation

**Problem:** All ID computations assumed `nRsectorsAxial=1`, `nCrystalsAxial=1`, `nModulesTransaxial=1`, `nSubmodulesTransaxial=1`. In a cubic array scanner each level (rsector, module, submodule, crystal) can have both axial and transaxial components encoded in a single GATE ID, so the old simplified formulas were incorrect for those cases.

Specific bugs:
- `ringID = moduleID + nModulesAxial * submoduleID` ignored `rsectorAxlID` and the axial crystal index, and used the wrong axis ordering (submodule as outer axis rather than module)
- `transaxialID` and `ringPosID` used raw `rsectorID` instead of its transaxial component, making them wrong whenever `nRsectorsAxial > 1`
- `maxRingID` and `nFanSumRings` were missing the `nRsectorsAxial * nCrystalsAxial` factors
- The block-counting condition did not check `rsectorAxlID` equality

**Solution:** Introduced four `static inline` helper functions in `normFunctions.cxx` that centralise the decomposition logic:

```cpp
rsectorAxlComponent(rsectorID, rsectorIdOrder, nRsectorsAngPos, nRsectorsAxial)
rsectorTrsComponent(rsectorID, rsectorIdOrder, nRsectorsAngPos, nRsectorsAxial)

computeRingID(rsectorAxl, moduleID, submoduleID, crystalID,
              nModulesAxial, nSubmodulesAxial, nCrystalsAxial,
              nModulesTransaxial, nSubmodulesTransaxial, nCrystalsTransaxial)
// = rsectorAxl * (nModAxl * nSubAxl * nCryAxl)
// + (moduleID / nModTrs) * (nSubAxl * nCryAxl)
// + (submoduleID / nSubTrs) * nCryAxl
// + (crystalID / nCryTrs)

computeTransaxialID(rsectorTrs, moduleID, submoduleID, crystalID,
                    nModulesTransaxial, nSubmodulesTransaxial, nCrystalsTransaxial)
// = rsectorTrs * (nModTrs * nSubTrs * nCryTrs)
// + (moduleID % nModTrs) * (nSubTrs * nCryTrs)
// + (submoduleID % nSubTrs) * nCryTrs
// + (crystalID % nCryTrs)
```

These helpers are used uniformly across `processSolidCyl_BlockCounts`, `processSolidCyl_GeomMatrix`, `processAnnular`, and the main CastorID-pair loop.

**Also fixed:** `nRsectorsAxial` and `rsectorIdOrder` added as explicit parameters to the three `processFile` functions (previously they had no way to decompose `rsectorID`).

**Dimensions corrected:**
```cpp
maxRingID    = nRsectorsAxial * nModulesAxial * nSubmodulesAxial * nCrystalsAxial
nFanSumRings = nRsectorsAxial * nModulesAxial * nSubmodulesAxial * nCrystalsAxial
```

**Regression test:** `tests/test_indexing_32x16x2.cxx` (no ROOT dependency, run via `ctest`) verifies for the 32x16x2_4rings scanner:
- `transaxialID`, `ringPosID`, `radialID`, and `trAID` are identical between old and new formulas (the transaxial pipeline was already correct for this scanner)
- The new `ringID` is in `[0, maxRingID)` for every detector and is consistent with the `ConvertIDcylindrical` / `ReverseCastorID` axis ordering (`moduleID * nSubmodulesAxial + submoduleID`)
- `maxRingID` and `nFanSumRings` are numerically unchanged for this scanner (`nRsectorsAxial * nCrystalsAxial = 1`)

**Note on ring ordering change for existing scanners:** The old formula `moduleID + nModulesAxial*submoduleID` interleaved rings from different modules, making axial Gaussian smoothing act on physically non-adjacent rings. The new formula `moduleID * nSubmodulesAxial + submoduleID` gives the physically correct sequential z-ordering. Without smoothing, normalization factors are identical. With smoothing enabled, results will differ slightly (and are now physically correct).

**Files Modified:**
- `src/normFunctions.cxx`: Added helper functions; updated all ID computations and `processFile` function bodies
- `include/normFunctions.h`: Updated `processSolidCyl_BlockCounts`, `processSolidCyl_GeomMatrix`, `processAnnular` signatures
- `tests/test_indexing_32x16x2.cxx`: New regression test
- `CMakeLists.txt`: Added test target and `enable_testing()`
