# PET Scanner Normalisation - Development Notes

## Current Constraints

The following assumptions apply to the current implementation:

- `nCrystalsAxial = 1` (always)
- `nSubmodulesTransaxial = 1` (always)

These constraints simplify the indexing:
- Ring indexing (`ringID = moduleID + nModulesAxial * submoduleID`) is valid because `nCrystalsAxial = 1`
- Transaxial strides don't need submodule transaxial contributions since `nSubmodulesTransaxial = 1`

## Indexing Overview

### Ring Indexing (Axial Component)
```c++
ringID = moduleID + nModulesAxial * submoduleID
maxRingID = nModulesAxial * nSubmodulesAxial * nCrystalsAxial
```
Note: If `nCrystalsAxial > 1` is needed in future, the formula must be generalised (marked with TODO in code).

### Transaxial Indexing (18D Packing)
Layers are excluded from transaxial strides by design:
```c++
strideRsector = nModulesTransaxial * nSubmodulesTransaxial * nCrystalsTransaxial
totalTransaxial = nRsectorsAngPos * nModulesTransaxial * nSubmodulesTransaxial * nCrystalsTransaxial
```

### Radial ID
Computed as the folded difference of transaxial positions:
```c++
delta = abs(ringPosID1 - ringPosID2)
radialID = min(delta, totalTransaxial - delta) - 1
maxRadialID = nRsectorsAngPos * nModulesTransaxial * nSubmodulesTransaxial * nCrystalsTransaxial / 2
```

### Transaxial-in-Radial ID (trAID)
```c++
trAID = ringPosID (before adding rsector contribution)
maxTrAID = nModulesTransaxial * nSubmodulesTransaxial * nCrystalsTransaxial
```

## Normalisation Components

Four independent factors are computed:

1. **Block Correction** (`ringComponentVector`): Per-ring coincidence counts
2. **Geometric Axial Correction** (`ringsComponentMatrix`): Per-ring-pair geometric factors
3. **Transaxial Geometric Correction** (`radialComponentVector`): Per-radialID groupings
4. **Interference/Intra-block Correction** (`blockTrAComponentMatrix`): Per-(radialID, trAID) pairs

Final normalisation:
```c++
CBasedNF = effNormFactor * norm_block * norm_geomAx * transaxialGeomNormFactor * interferenceTraFactor
```

## Detector ID Swap Logic

When processing coincidences, detector pairs are ordered by CastorID:
```c++
if (castorID1 < castorID2) {
    std::swap(castorID1, castorID2);
    std::swap(layerID1, layerID2);
    std::swap(crystalID1, crystalID2);
    std::swap(submoduleID1, submoduleID2);
    std::swap(moduleID1, moduleID2);
    std::swap(rsectorID1, rsectorID2);
    std::swap(gPos1, gPos2);
}
```
All component IDs are swapped together, maintaining correspondence between CastorID and its components.

## OpenMP Parallelisation

The code uses OpenMP for parallel execution on multi-core systems.

### Build Requirements
```cmake
find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)
```

### Parallelised Sections (Phase 1 - Implemented)

1. **LUT Building** (`buildCastorIDLUT`):
   - `#pragma omp parallel for collapse(4) schedule(static)`
   - Each iteration writes to a unique index (no contention)

2. **Global Mean Reductions** (`computeNormalizationFactors`):
   - Block/geometric axial means: `collapse(2) reduction(+:sb_total,sg_total,...)`
   - Transaxial mean: `reduction(+:s_t,c_t)`
   - Interference mean: `reduction(+:si_total,ci_total)`

3. **Transaxial Normalisation Vector**:
   - `#pragma omp parallel for` (independent writes to different indices)

### Controlling Thread Count
```bash
export OMP_NUM_THREADS=4  # Use 4 threads
./own-normFactors ...
```

### Future Parallelisation Candidates

**Phase 2** (not yet implemented):
- 9-level nested loop in `computeNormalizationFactors` (main computation)
  - Requires buffered file I/O due to thread-safety
  - Expected speedup: 2-3x

**Phase 3** (not yet implemented):
- TTree processing in `processFile`
  - ROOT's TTree::GetEntry() is not thread-safe
  - Requires buffered approach: read sequentially, process in parallel
  - Expected speedup: 1.5-2x

## Future Considerations

If constraints change:
- `nCrystalsAxial > 1`: Update ring indexing formula to include crystal axial term
- `nSubmodulesTransaxial > 1`: Verify stride calculations account for submodule transaxial
- Layer inclusion in transaxial: Uncomment layer stride contribution (lines 285-288) and update maxTrAID
