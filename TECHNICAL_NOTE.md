Technical note: component-based normalisation of PET data from GATE simulations (own-normFactors)
==================================================================================================

Author: Marc Granado-Gonzalez
Status: DRAFT (skeleton), section 9 updated 2026-10-05
Code:   src/own-normFactors.cxx, src/normFunctions.cxx, src/scannerGeometry.cxx


1. Purpose
----------
- Compute a normalisation factor for every LOR of a cylindrical PET scanner, using
  GATE Monte Carlo simulations of normalisation phantoms.
- The output is a CASToR normalisation file (.Cdh/.Cdf), ready to use in reconstruction.
- TODO: motivation (why a component-based method instead of direct inversion: statistics per LOR).


2. Inputs
---------
2.1 Scanner geometry
    - CASToR cylindrical scanner file (.geom), read by ScannerGeometry.
    - Crystal centres come from the file; nothing is hardcoded (see README, "How crystal positions are derived").
    - Per layer: scanner radius, crystal depth, and the transaxial crystal count, size and gap
      (a single value applies to all layers). Everything else (axial crystal count, axial sizes
      and gaps, module/submodule/rsector counts) must be the same for all layers. See section 9.
    - GATE ID -> CASToR ID mapping (ConvertIDcylindrical / ReverseCastorID), with the options
      --invert-det-order and --rsector-id-order.
    - Transaxial FOV radius, option -f/--fov-radius (default 300 mm, must be smaller than the
      scanner radius). See 4.7.
2.2 Simulated data (GATE Coincidences tree, ROOT)
    - Solid cylinder source ("solidCyl" in the file name): block, axial geometric and crystal efficiency.
    - Annular source ("annular" in the file name): transaxial geometric and interference.
    - TODO: phantom dimensions, activity, simulated statistics, physics used.


3. Index definitions
--------------------
- ringID: axial position of a crystal (rsector axial, module, submodule, crystal axial).
- transaxialID: position of a crystal around the ring, counted in the crystals of its own layer
  (layers excluded). Used by the fan sums.
- ringPos: position of a crystal around the ring on the common transaxial grid (section 9.2).
  With the same transaxial crystal count in every layer it is equal to transaxialID.
- radialID: min(delta, L - delta) - 1, with delta = |ringPos1 - ringPos2| and L the number of
  grid points per ring. It is the transaxial distance of the LOR to the axis.
- trAID: transaxial position of crystal 1 inside its rsector, counted in the crystals of its own layer.
- TODO: figure with the indices on a sketch of the scanner.


4. Model
--------
   N(LOR) = eff x B x G_ax x G_trs x I

4.1 Crystal efficiency, eff (fan sum, Pepin et al. 2011, eq. 4)
    - Fan of crystal (u,i) = all coincidences involving it in the solidCyl data.
    - eff = (avgFan_u1 x avgFan_u2) / (fan_u1i1 x fan_u2i2).
    - Layers: with the same transaxial count in every layer, crystal i of all layers shares one fan
      (layers are merged). Otherwise each layer has its own fan-sum counter and ring average.
4.2 Block correction, B (per ring)
    - Pass 1 on the solidCyl data: count coincidences between crystals of the same block, per ring.
    - B = mean / sqrt(c_ring1 x c_ring2).
4.3 Axial geometric correction, G_ax (per ring pair)
    - Pass 2 on the solidCyl data: accumulate cos(theta) x B per (ring1, ring2).
    - G_ax = mean / M[ring1][ring2].
4.4 Transaxial geometric correction, G_trs (per radialID)
    - Annular data, each LOR weighted by eff x B x G_ax / (phantom line integral).
    - Line integral = full phantom minus inner empty cylinder; only LORs crossing the annulus are kept.
    - Bin delta = L/2 has no folding partner, so it is weighted x2.
    - G_trs = mean / R[radialID].
4.5 Interference correction, I (per radialID, trAID)
    - Same accumulation as 4.4, split by the position of the crystal in its block.
    - I = mean / (G_trs x R[radialID][trAID]).
    - Layers: one matrix shared by all layers when they have the same transaxial count, otherwise
      one matrix per layer (indexed by the layer of crystal 1). The geometric mean is taken over
      all matrices together, so differences between layers stay in I.
4.6 Means and regularisation
    - All "mean" values are geometric means (less sensitive to skewed counts than the arithmetic mean).
    - Optional Gaussian smoothing: -as on the per-ring block counts (before pass 2),
      -ts on the radial vector (after the annular data).
    - TODO: justify the default sigma values with results.
    - -ts is not recommended when the layers have different transaxial crystal counts (9.3).
4.7 FOV selection
    - A LOR is used (solidCyl pass 2, annular) and written (output loop) only when the segment
      between the two crystal centres crosses the transaxial FOV: its distance R to the axis is
      at most the FOV radius, and the closest point to the axis lies between the two crystals.
    - The second condition removes pairs on the same side of the ring (e.g. two layers of one
      block), whose line can pass close to the axis without the segment crossing the FOV.
    - The same rule is used in every pass, so the LORs written are exactly those that the
      correction factors were built from.


5. Processing pipeline
----------------------
   1. Read the .geom file and command line options.
   2. solidCyl pass 1 -> per-ring block counts + fan sums.
   3. (optional) axial smoothing.
   4. solidCyl pass 2 -> axial geometric matrix.
   5. annular -> radial vector + interference matrix.
   6. (optional) transaxial smoothing.
   7. Loop over all CASToR ID pairs (OpenMP) -> write the normalisation files and CSV diagnostics.
- TODO: running times and memory for each provided scanner.


6. Outputs
----------
- _CB             : eff x B x G_ax x G_trs x I   (standard)
- _effBfGAfGTrAf  : eff x B x G_ax x G_trs       (no interference)
- _CBsqrBfnI      : eff x B^2 x G_ax x G_trs     (squared block, no interference)
- CSV diagnostics: block_geom, transaxial, interference, fan_sum_efficiency.
  interference and fan_sum_efficiency get a leading "layer" column only when the layers have
  different transaxial crystal counts; otherwise the format is unchanged.


7. Validation
-------------
- TODO: uniform cylinder reconstructed with and without normalisation (uniformity, axial profile).
- TODO: comparison between the three outputs.
- TODO: effect of statistics and of smoothing (axial ripple, central LORs).
- Regression tests: tests/test_scanner_geometry.cxx (crystal positions, per-layer counts, common grid).
- Per-layer transaxial counts (2026-10-05), on synthetic Coincidences trees (no real GATE data yet):
    * equal counts (2 layers x 8 crystals, 16 rsectors): all .Cdf/.Cdh/CSV outputs byte-identical
      to the previous version, for the default options, --invert-det-order and --rsector-id-order 1;
    * 8 and 2 crystals: run completes, the 640 crystals give 640 distinct castorIDs and every
      castorID converts back to the same crystal (also checked on 32x16x2_4rings, 131072 crystals).
    * TODO: validate the mixed-count case on GATE data.
- FOV selection (2026-10-05), same synthetic data: the LORs written were compared with an
  independent enumeration of all crystal pairs whose segment crosses R = 300 mm: no missing and
  no extra LOR (384256 LORs with equal counts, 152832 with 8 and 2 crystals). On the LORs common
  with the previous version all factors are identical, and the CSV diagnostics are unchanged.


8. Known limitations
--------------------
- The GATE -> geometry mapping is exact only when there is 1 module and 1 submodule transaxially,
  1 crystal axially and 1 rsector axially (a warning is printed otherwise).
- solidCyl pass 1 has no FOV selection (9.4, open).
- Mixed transaxial counts give a radial vector that alternates bin to bin (section 9.3); it is
  correct as is but must not be smoothed with -ts.
- Both processing stages run together: changing a smoothing option means re-reading all the ROOT files.
- TODO: add others.


9. Layers with different transaxial crystal counts
--------------------------------------------------
9.1 Before (up to commit 437d672)
    - Every layer had to have the same number of elements; only "scanner radius" and
      "crystals size depth" could differ between layers (the parser rejected anything else).
    - All indices used one scalar nCrystalsTransaxial. Layers were excluded from the transaxial
      indices, so crystal i of layer 0 and crystal i of layer 1 were the same transaxial position:
      they shared the same fan, the same trAID column and the same radialID arithmetic.
    - The castorID layer offset was nCrystalPerLayer[layer-1] x layer, correct only for equal layers.

9.2 Now
    - "number of crystals transaxial", "crystals size trans" and "crystal gap transaxial" are read
      per layer. Crystal positions use the count and pitch of their own layer.
    - Common transaxial grid (ScannerGeometry::gridPoint). With L the least common multiple of the
      per-layer counts and r = L / n[l]:
        * all r odd: L points per submodule, crystal j of layer l at point j*r + (r-1)/2;
        * otherwise: 2L points per submodule, crystal j of layer l at point (2j+1)*r.
      In both cases every crystal centre is exactly a grid point (no rounding). Example: 10 and 2
      crystals -> 10 points, the large crystals at points 2 and 7; 10 and 5 -> 20 points.
      With equal counts the grid point is the crystal index (step 1, offset 0).
    - radialID uses ringPos on this grid, so crystals of different layers can be paired.
      trAID and transaxialID stay in the crystals of each layer, with one interference matrix and
      one fan-sum counter per layer (4.1, 4.5).
    - castorID: all crystals of a layer before the next one; the layer offset is the sum of the
      crystal counts of the previous layers, and each layer uses its own transaxial count.
    - Same code path for all scanners: with equal counts, the grid is the identity and the
      per-layer containers hold a single shared entry, so the results are unchanged (section 7).
    - Width mismatch: the grid assumes the crystals of every layer span the same transaxial width
      (count x pitch). A difference d moves the outer crystals by up to d/2, zero at the block
      centre (a stretch, not a constant offset, so it cannot be corrected by a constant shift).
      ScannerGeometry::gridMismatch gives this bound in grid spacings and a warning is printed
      above 0.25. Example: 59.0 mm against 58.8 mm with 3.6875 mm crystals -> about 0.03 spacing.
      Equal indices mean equal transaxial position in the flat block, not equal angle; the
      different layer radii were already ignored the same way before this change.

9.3 Alternating radial vector (mixed counts)
    - When the grid is doubled (some r even), the crystals of the finest layer only sit on odd
      grid points and those of the coarser layers on even points (8 and 2 crystals: points
      1, 3, ..., 15 and points 4, 12). The distance between two crystals of the same parity is
      even, between different parities odd. So odd and even radialID bins are filled by
      different populations: same-layer pairs on one parity, cross-layer pairs on the other.
    - The bins are not empty, but the two populations differ in statistics and efficiency, so
      the radial vector alternates bin to bin. On synthetic 8/2 data the mean G_trs factor was
      0.94 on even bins and 1.65 on odd bins (the size of the effect on synthetic data has no
      physical meaning; the alternation itself does).
    - Without the doubled grid (all r odd, e.g. 10 and 2) there is no strict parity split, but the
      mix of layer pairs still changes from bin to bin with the period of the coarse step.
    - Decision: keep the radial vector as it is. Neighbouring bins differ for a geometric reason
      (they are LORs of different layer pairs), so each bin is a correct factor for the LORs that
      use it; only the ordering of the bins makes the vector look irregular. Separate radial
      vectors per layer pair were considered and not needed.
    - Consequence: the -ts Gaussian smoothing must not be used for these scanners, since it
      averages neighbouring bins, i.e. it mixes the layer pairs. A warning is printed when -ts > 0
      and the layers have different transaxial counts.

9.4 LOR selection of the output loop (rsector distance and the 13-crystal cut) -- replaced
    - Before 2026-10-05 the output loop skipped pairs whose rsectors were less than 4 apart, and at exactly 4 apart it
      skips a pair when |submoduleID1 - submoduleID2| < 13, or else when |crystalID1 - crystalID2| < 13.
    - For the shipped 32-rsector scanners (radius 323.8 mm, 32 x 1.84375 mm crystals) the intent
      appears to be the |R| < 300 mm FOV of the annular stage:
        * rsectors 3 apart: R = 301.6 - 318.2 mm, all outside -> skipping them is right;
        * rsectors 5 or more apart: R <= 299.0 mm, all inside -> keeping them is right;
        * rsectors 4 apart: R = 288.2 - 310.1 mm, the FOV edge (R = 300 mm needs an angular
          separation of 44.2 deg, and 4 rsectors are 45 deg apart).
    - 13 crystals = 23.97 mm = 4.23 deg seen from the axis (layer 0), 41% of the 10.4 deg block.
    - The crystal test does not match the FOV: R depends on the signed offset between the two
      crystals (which one is closer to the other rsector), but the cut uses |crystalID1 - crystalID2|.
      For layer 0 at 4 rsectors apart: 553 of the 1024 crystal pairs are inside R < 300 mm,
      the cut keeps 380, of which 190 are outside, and drops 363 that are inside.
    - The submodule test: in all shipped scanners the submodules are axial, so it compares axial
      positions, which do not change R. (nSubmodulesTransaxial != 0 is always true.)
    - With different transaxial counts per layer, |crystalID1 - crystalID2| compared indices of
      different crystal sizes.
    - The rsector distances were also hardcoded for 32 rsectors: on a 16-rsector test scanner
      (4 rsectors = 90 deg) the old rule dropped 154880 of the 384256 LORs inside the FOV.
    - Now: replaced by the geometric FOV selection of 4.7 with a user-given FOV radius. It works
      for any crystal size, layer radius and number of rsectors, and the sign question disappears
      because R is computed from the actual crystal positions. This changes the set of LORs
      written for every scanner (boundary pairs added or removed). The solidCyl pass 2 and the
      annular pass used |R| <= 300 mm before; they now also drop pairs on the same side of the
      ring (4.7). Those pairs are mostly coincidences between layers of one block (or of two
      neighbouring blocks): crystals less than ~26 mm apart on a line pointing near the axis.
      With |R| only they were accepted (two layers of one crystal column give R = 0), and the
      annular line integral, which is computed on the infinite line and not on the segment,
      gave them a non-zero weight.
    - Effect, synthetic data with 25% of the coincidences inside one rsector (equal counts,
      16 rsectors): before, those pairs filled radialID 0, which no written LOR uses, and shifted
      the geometric mean so that every G_trs was ~1.7% low; the same-ring G_ax entries were ~4%
      low (0.956-0.964 -> 0.989-0.998 after), since same-ring pairs have cos(theta) = 1. On the
      LORs written by both versions the final factor changed by +1.8% on average (-0.8% to +5.7%).
      The size scales with the fraction of such coincidences in the data; the direction does not.
    - Open: solidCyl pass 1 (block counts and fan sums) still uses all coincidences, without FOV
      selection. Whether same-rsector coincidences bias the crystal efficiencies depends on their
      origin (inter-crystal scatter with crystal-level readout varies by crystal; randoms mostly
      cancel in the fan-sum ratio). TODO: measure their fraction per crystal and layer pair on
      the GATE data, then decide whether pass 1 uses the FOV selection too.


10. References
-------------
- Pepin et al., "Normalization of Monte Carlo PET data using GATE", IEEE NSS/MIC 2011.
- CASToR documentation (scanner file format, normalisation file format).
- GATE documentation (cylindricalPET system, Coincidences tree).
