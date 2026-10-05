#ifndef SCANNER_GEOMETRY_H
#define SCANNER_GEOMETRY_H

// Scanner geometry read from a CASToR cylindrical-PET scanner file (.geom).
//
// This module has no ROOT dependency. Field names are the CASToR ones
// ("scanner radius", "number of rsectors", "crystals size trans", ...).
//
// Conventions (same as CASToR):
//   * All lengths in mm, angles in degrees.
//   * Gaps are edge to edge (free space between two neighbouring elements).
//   * Layer-dependent fields are comma separated lists, one value per layer.
//
// In this tool a "layer" is a depth segmentation of the detector.
// "scanner radius" and "crystals size depth" must list one value per layer.
// The transaxial crystal count, size and gap may differ between layers (a single
// value is applied to all layers). All other counts, sizes and gaps are identical
// for all layers.

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

struct Position3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct ScannerGeometry {
    std::string name;

    // --- Counts (identical for all layers) ---------------------------------
    uint32_t nLayers = 0;
    uint32_t nRsectorsAngPos = 0;      // "number of rsectors"
    uint32_t nRsectorsAxial = 1;       // "number of rsectors axial"
    uint32_t nModulesTransaxial = 1;   // "number of modules transaxial"
    uint32_t nModulesAxial = 1;        // "number of modules axial"
    uint32_t nSubmodulesTransaxial = 1;// "number of submodules transaxial"
    uint32_t nSubmodulesAxial = 1;     // "number of submodules axial"

    uint32_t nCrystalsAxial = 0;       // "number of crystals axial"

    // --- Sizes and gaps (mm) -----------------------------------------------
    double crystalSizeAxial = 0.0;     // "crystals size axial"
    double crystalGapAxial = 0.0;      // "crystal gap axial"
    double submoduleGapTrans = 0.0;    // "submodule gap transaxial"
    double submoduleGapAxial = 0.0;    // "submodule gap axial"
    double moduleGapTrans = 0.0;       // "module gap transaxial"
    double moduleGapAxial = 0.0;       // "module gap axial"
    double rsectorGapAxial = 0.0;      // "rsector gap axial"

    // --- Layer-dependent values (one entry per layer) ----------------------
    std::vector<uint32_t> nCrystalsTransaxial; // "number of crystals transaxial"
    std::vector<double> crystalSizeTrans;      // "crystals size trans"
    std::vector<double> crystalGapTrans;       // "crystal gap transaxial"
    std::vector<double> layerRadius;   // "scanner radius": isocentre -> front face of the layer
    std::vector<double> layerDepth;    // "crystals size depth"

    // --- Rsector placement (degrees) ---------------------------------------
    double rsectorFirstAngleDeg = 0.0; // "rsectors first angle"
    double rsectorAngularSpanDeg = 360.0; // "rsectors angular span"

    // --- Derived values (filled by finalize()) -----------------------------
    double crystalPitchAxial = 0.0;    // crystal size + crystal gap
    double submodulePitchAxial = 0.0;  // axial distance between submodule centres
    double modulePitchAxial = 0.0;     // axial distance between module centres
    double moduleSizeAxial = 0.0;      // axial extent of one module
    double rsectorStepRad = 0.0;       // angle between two consecutive rsectors
    std::vector<double> crystalPitchTrans; // per layer: crystal size + crystal gap
    std::vector<double> moduleSizeTrans;   // per layer: transaxial extent of the crystals of one rsector
    std::vector<double> layerCentreRadius; // layerRadius + layerDepth / 2
    std::vector<double> cosRsector;
    std::vector<double> sinRsector;

    // --- Common transaxial grid (filled by finalize()) ---------------------
    // Transaxial crystal centres of all layers expressed as integer positions
    // on one evenly spaced grid, so that crystals of layers with different
    // transaxial counts can be compared (radial distance of a LOR).
    // A submodule holds gridPointsPerSubmodule grid points; the centre of
    // crystal j of layer l is grid point j * gridStep[l] + gridOffset[l].
    // With L the least common multiple of the counts:
    //   * all L / n[l] odd: L points, crystal centres fall on point centres
    //   * otherwise:        2L points, so that every crystal centre is a point
    // When all layers have the same count the grid point is the crystal index
    // itself (step 1, offset 0), so nothing changes for those scanners.
    bool uniformCrystalsTransaxial = true;
    uint32_t gridPointsPerSubmodule = 0;
    std::vector<uint32_t> gridStep;
    std::vector<uint32_t> gridOffset;
    // Upper bound, in grid spacings, of the distance between the grid point of
    // a crystal and its actual position. It is non zero when the crystals of
    // the layers do not span the same transaxial width (0 for uniform counts).
    double gridMismatch = 0.0;

    // Non fatal remarks collected while parsing / finalizing.
    std::vector<std::string> warnings;

    // Checks the values and computes the derived ones.
    // Returns false and fills `error` if the geometry is not usable.
    bool finalize(std::string &error);

    // Fills the common transaxial grid (called by finalize()).
    bool buildTransaxialGrid(std::string &error);

    // Largest transaxial crystal count over all layers.
    uint32_t maxNCrystalsTransaxial() const;

    // Number of crystals in one layer of the whole scanner.
    uint32_t nCrystalsInLayer(uint32_t layerID) const;

    // Grid point of the centre of transaxial crystal crystalTrs
    // (0 .. nCrystalsTransaxial[layerID]-1) of layer layerID.
    uint32_t gridPoint(uint32_t layerID, uint32_t crystalTrs) const {
        return crystalTrs * gridStep[layerID] + gridOffset[layerID];
    }

    // Centre of a crystal in the scanner frame (z = 0 at the axial centre of the scanner).
    //
    // The IDs are the ones stored in the GATE Coincidences tree. They are mapped as
    //   crystalID   -> transaxial crystal index  (within one rsector)
    //   submoduleID -> axial submodule index     (within one module)
    //   moduleID    -> axial module index
    //   rsectorID   -> angular rsector index
    // which is only meaningful when nModulesTransaxial = nSubmodulesTransaxial =
    // nCrystalsAxial = nRsectorsAxial = 1 (a warning is recorded otherwise).
    Position3 crystalPosition(int layerID, int crystalID, int submoduleID,
                              int moduleID, int rsectorID) const;

    void print(std::ostream &os) const;
};

// Reads a CASToR cylindrical PET scanner file.
//   * "label: value" lines, text after '#' is a comment, unknown labels are ignored
//   * label matching is case insensitive
// Returns false and fills `error` on failure. On success `geom` is finalized.
bool ParseCastorGeom(const std::string &path, ScannerGeometry &geom, std::string &error);

// Same, from a string (used by the tests).
bool ParseCastorGeomText(const std::string &text, const std::string &defaultName,
                         ScannerGeometry &geom, std::string &error);

#endif // SCANNER_GEOMETRY_H
