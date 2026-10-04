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
// In this tool a "layer" is a depth segmentation of the detector: every layer
// has the same number of elements, only "scanner radius" (front face of the
// layer) and "crystals size depth" may differ from layer to layer.

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
    uint32_t nCrystalsTransaxial = 0;  // "number of crystals transaxial"
    uint32_t nCrystalsAxial = 0;       // "number of crystals axial"

    // --- Sizes and gaps (mm) -----------------------------------------------
    double crystalSizeTrans = 0.0;     // "crystals size trans"
    double crystalSizeAxial = 0.0;     // "crystals size axial"
    double crystalGapTrans = 0.0;      // "crystal gap transaxial"
    double crystalGapAxial = 0.0;      // "crystal gap axial"
    double submoduleGapTrans = 0.0;    // "submodule gap transaxial"
    double submoduleGapAxial = 0.0;    // "submodule gap axial"
    double moduleGapTrans = 0.0;       // "module gap transaxial"
    double moduleGapAxial = 0.0;       // "module gap axial"
    double rsectorGapAxial = 0.0;      // "rsector gap axial"

    // --- Layer-dependent values (one entry per layer) ----------------------
    std::vector<double> layerRadius;   // "scanner radius": isocentre -> front face of the layer
    std::vector<double> layerDepth;    // "crystals size depth"

    // --- Rsector placement (degrees) ---------------------------------------
    double rsectorFirstAngleDeg = 0.0; // "rsectors first angle"
    double rsectorAngularSpanDeg = 360.0; // "rsectors angular span"

    // --- Derived values (filled by finalize()) -----------------------------
    double crystalPitchTrans = 0.0;    // crystal size + crystal gap
    double crystalPitchAxial = 0.0;
    double submodulePitchAxial = 0.0;  // axial distance between submodule centres
    double modulePitchAxial = 0.0;     // axial distance between module centres
    double moduleSizeTrans = 0.0;      // transaxial extent of the crystals of one rsector
    double moduleSizeAxial = 0.0;      // axial extent of one module
    double rsectorStepRad = 0.0;       // angle between two consecutive rsectors
    std::vector<double> layerCentreRadius; // layerRadius + layerDepth / 2
    std::vector<double> cosRsector;
    std::vector<double> sinRsector;

    // Non fatal remarks collected while parsing / finalizing.
    std::vector<std::string> warnings;

    // Checks the values and computes the derived ones.
    // Returns false and fills `error` if the geometry is not usable.
    bool finalize(std::string &error);

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
