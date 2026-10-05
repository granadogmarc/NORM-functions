// Tests for the CASToR .geom parser and the crystal positions (no ROOT needed).

#include "scannerGeometry.h"

#include <cmath>
#include <cstdio>
#include <string>

#ifndef SCANNER_DIR
#error "SCANNER_DIR must point to the scanners directory"
#endif

static int failures = 0;

#define CHECK(cond)                                                              \
    do {                                                                         \
        if (!(cond)) {                                                           \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);          \
            ++failures;                                                          \
        }                                                                        \
    } while (0)

#define CHECK_NEAR(a, b, tol)                                                    \
    do {                                                                         \
        const double va = (a), vb = (b);                                         \
        if (!(std::fabs(va - vb) <= (tol))) {                                    \
            std::printf("FAIL %s:%d: %s = %.9g, expected %.9g\n", __FILE__,      \
                        __LINE__, #a, va, vb);                                   \
            ++failures;                                                          \
        }                                                                        \
    } while (0)

static const char *kBase =
    "modality: PET\n"
    "scanner name: test\n"
    "number of layers: 2\n"
    "scanner radius: 321.3, 326.3\n"
    "number of rsectors: 32, 32\n"
    "number of modules axial: 4, 4\n"
    "number of submodules axial: 32, 32\n"
    "number of crystals transaxial: 32, 32\n"
    "number of crystals axial: 1, 1\n"
    "crystals size trans: 1.84375, 1.84375\n"
    "crystals size axial: 1.84375, 1.84375\n"
    "crystals size depth: 5, 5\n"
    "module gap axial: 4, 4\n";

static bool parse(const std::string &text, ScannerGeometry &g, std::string &err) {
    return ParseCastorGeomText(text, "default", g, err);
}

static void testShippedScanners() {
    const char *files[] = {"16x16x2_1ring", "16x16x2_4rings", "32x16x2_4rings",
                           "32x32x2_4rings", "CM2L_1ring"};
    for (const char *f : files) {
        ScannerGeometry g;
        std::string err;
        const std::string path = std::string(SCANNER_DIR) + "/" + f + ".geom";
        const bool ok = ParseCastorGeom(path, g, err);
        if (!ok) std::printf("  %s: %s\n", path.c_str(), err.c_str());
        CHECK(ok);
        if (!ok) continue;
        // Virtual segmentation of a 59 x 59 x 10 mm module
        CHECK_NEAR(g.moduleSizeTrans[0], 59.0, 1e-9);
        CHECK_NEAR(g.moduleSizeTrans[1], 59.0, 1e-9);
        CHECK_NEAR(g.moduleSizeAxial, 59.0, 1e-9);
        CHECK_NEAR(g.modulePitchAxial, 63.0, 1e-9);
        CHECK_NEAR(g.layerCentreRadius[0], 323.8, 1e-9);
        CHECK_NEAR(g.layerCentreRadius[1], 328.8, 1e-9);
        CHECK_NEAR(g.rsectorStepRad, 2.0 * M_PI / 32.0, 1e-12);
        CHECK(g.warnings.empty());
        // Uniform counts: the common grid is the crystal index itself
        CHECK(g.uniformCrystalsTransaxial);
        CHECK(g.gridPointsPerSubmodule == g.nCrystalsTransaxial[0]);
        for (uint32_t l = 0; l < g.nLayers; ++l)
            for (uint32_t j = 0; j < g.nCrystalsTransaxial[l]; ++j) CHECK(g.gridPoint(l, j) == j);
        CHECK_NEAR(g.gridMismatch, 0.0, 0.0);
    }
}

static void testPositions32x32x2() {
    ScannerGeometry g;
    std::string err;
    CHECK(ParseCastorGeom(std::string(SCANNER_DIR) + "/32x32x2_4rings.geom", g, err));

    // First crystal: layer 0, crystal 0, submodule 0, module 0, rsector 0
    Position3 p = g.crystalPosition(0, 0, 0, 0, 0);
    CHECK_NEAR(p.x, 323.8, 1e-9);
    CHECK_NEAR(p.y, -(59.0 / 2.0 - 1.84375 / 2.0), 1e-9);   // -28.578125
    CHECK_NEAR(p.z, -(59.0 / 2.0 - 1.84375 / 2.0) - 1.5 * 63.0, 1e-9); // -123.078125

    // Second layer is 5 mm further out
    CHECK_NEAR(g.crystalPosition(1, 0, 0, 0, 0).x, 328.8, 1e-9);

    // The whole stack is centred on z = 0 and on the rsector axis
    const Position3 last = g.crystalPosition(0, 31, 31, 3, 0);
    CHECK_NEAR(last.z, -p.z, 1e-9);
    CHECK_NEAR(last.y, -p.y, 1e-9);

    // Pitches: crystals 1.84375 mm, modules 63 mm
    CHECK_NEAR(g.crystalPosition(0, 1, 0, 0, 0).y - p.y, 1.84375, 1e-9);
    CHECK_NEAR(g.crystalPosition(0, 0, 1, 0, 0).z - p.z, 1.84375, 1e-9);
    CHECK_NEAR(g.crystalPosition(0, 0, 0, 1, 0).z - p.z, 63.0, 1e-9);

    // Rsector 8 of 32 is rotated by 90 degrees: (x, y) -> (-y, x)
    const Position3 r = g.crystalPosition(0, 0, 0, 0, 8);
    CHECK_NEAR(r.x, -p.y, 1e-9);
    CHECK_NEAR(r.y, p.x, 1e-9);
    CHECK_NEAR(r.z, p.z, 1e-12);
}

static void testGapsAreEdgeToEdge() {
    ScannerGeometry g;
    std::string err;
    const std::string text = std::string(kBase) +
        "crystal gap transaxial: 0.2, 0.2\n"
        "submodule gap axial: 0.5, 0.5\n";
    CHECK(parse(text, g, err));
    CHECK_NEAR(g.crystalPitchTrans[0], 1.84375 + 0.2, 1e-12);
    CHECK_NEAR(g.submodulePitchAxial, 1.84375 + 0.5, 1e-12);
    // Module = 32 submodules + 31 gaps; pitch = module + module gap
    CHECK_NEAR(g.moduleSizeAxial, 32 * 1.84375 + 31 * 0.5, 1e-9);
    CHECK_NEAR(g.modulePitchAxial, 32 * 1.84375 + 31 * 0.5 + 4.0, 1e-9);
    CHECK_NEAR(g.moduleSizeTrans[1], 32 * 1.84375 + 31 * 0.2, 1e-9);
}

// Replaces one "label: ..." line of kBase
static std::string withField(const std::string &label, const std::string &value) {
    std::string text = kBase;
    const size_t b = text.find(label + ":");
    const size_t e = text.find('\n', b);
    text.replace(b, e - b, label + ": " + value);
    return text;
}

static void testPerLayerTransaxialCounts() {
    ScannerGeometry g;
    std::string err;

    // 10 and 2 crystals over the same 20 mm: L = 10, ratios 1 and 5 are odd,
    // so the grid is the 10-crystal layer and the large crystals sit on points 2 and 7
    CHECK(parse(withField("number of crystals transaxial", "10, 2") +
                "crystals size trans: 2, 10\n", g, err));
    CHECK(!g.uniformCrystalsTransaxial);
    CHECK(g.maxNCrystalsTransaxial() == 10);
    CHECK(g.gridPointsPerSubmodule == 10);
    CHECK(g.gridPoint(0, 3) == 3);
    CHECK(g.gridPoint(1, 0) == 2);
    CHECK(g.gridPoint(1, 1) == 7);
    CHECK_NEAR(g.gridMismatch, 0.0, 1e-12);
    CHECK(g.warnings.empty());
    CHECK(g.nCrystalsInLayer(0) == 32u * 4 * 32 * 10);
    CHECK(g.nCrystalsInLayer(1) == 32u * 4 * 32 * 2);
    // Per-layer positions: outer crystals of both layers are centred the same
    CHECK_NEAR(g.crystalPosition(0, 0, 0, 0, 0).y, -9.0, 1e-9);
    CHECK_NEAR(g.crystalPosition(1, 0, 0, 0, 0).y, -5.0, 1e-9);
    CHECK_NEAR(g.crystalPosition(1, 1, 0, 0, 0).y, 5.0, 1e-9);

    // 10 and 5: ratio 2 is even, so the grid doubles (20 points) and every
    // crystal centre is a point: the 2 mm crystal j at 2j+1, the 4 mm crystal j at 4j+2
    CHECK(parse(withField("number of crystals transaxial", "10, 5") +
                "crystals size trans: 2, 4\n", g, err));
    CHECK(g.gridPointsPerSubmodule == 20);
    CHECK(g.gridPoint(0, 0) == 1);
    CHECK(g.gridPoint(0, 9) == 19);
    CHECK(g.gridPoint(1, 0) == 2);
    CHECK(g.gridPoint(1, 4) == 18);
    // Each grid point is half a small crystal (1 mm), and the 4 mm crystal 0
    // lies between small crystals 0 and 1 (points 1 and 3)
    CHECK(g.gridPoint(1, 0) - g.gridPoint(0, 0) == g.gridPoint(0, 1) - g.gridPoint(1, 0));

    // Different widths: 10 x 2 mm = 20 mm against 2 x 9.9 mm = 19.8 mm.
    // Edge crystals are off by 0.1 mm, i.e. 0.05 of the 2 mm grid spacing (mean width / 10)
    CHECK(parse(withField("number of crystals transaxial", "10, 2") +
                "crystals size trans: 2, 9.9\n", g, err));
    CHECK_NEAR(g.gridMismatch, 0.1 / (19.9 / 10.0), 1e-9);
    CHECK(g.warnings.empty());

    // Large width difference is reported
    CHECK(parse(withField("number of crystals transaxial", "10, 2") +
                "crystals size trans: 2, 8\n", g, err));
    CHECK(g.gridMismatch > 0.25);
    CHECK(!g.warnings.empty());

    // A single value applies to all layers
    CHECK(parse(withField("number of crystals transaxial", "16"), g, err));
    CHECK(g.nCrystalsTransaxial.size() == 2 && g.nCrystalsTransaxial[1] == 16);

    // Per-layer gap
    CHECK(parse(std::string(kBase) + "crystal gap transaxial: 0, 0.2\n", g, err));
    CHECK_NEAR(g.crystalPitchTrans[0], 1.84375, 1e-12);
    CHECK_NEAR(g.crystalPitchTrans[1], 1.84375 + 0.2, 1e-12);

    // The axial crystal count must still be the same for all layers
    CHECK(!parse(withField("number of crystals axial", "1, 2"), g, err));
    CHECK(err.find("differs between layers") != std::string::npos);

    // Wrong number of values
    CHECK(!parse(withField("number of crystals transaxial", "10, 2, 2"), g, err));
}

static void testParsingDetails() {
    ScannerGeometry g;
    std::string err;

    // Comments, extra CASToR fields, upper case labels and a single value for a
    // layer-dependent field (applied to all layers)
    const std::string text =
        "# comment line\n"
        "MODALITY: pet   # trailing comment\n"
        "Scanner Name: my scanner\n"
        "number of elements: 123456\n"
        "voxels number transaxial: 100\n"
        "number of layers: 2\n"
        "scanner radius: 321.3, 326.3\n"
        "number of rsectors: 16\n"
        "number of crystals transaxial: 8\n"
        "number of crystals axial: 1\n"
        "crystals size trans: 2\n"
        "crystals size axial: 2\n"
        "crystals size depth: 5, 5\n";
    CHECK(parse(text, g, err));
    CHECK(g.name == "my scanner");
    CHECK(g.nRsectorsAngPos == 16);
    CHECK(g.nModulesAxial == 1);              // default
    CHECK_NEAR(g.moduleGapAxial, 0.0, 0.0);   // default
    CHECK_NEAR(g.rsectorAngularSpanDeg, 360.0, 0.0);

    // Default name comes from the caller when no name is given
    std::string noName = text;
    noName.replace(noName.find("Scanner Name: my scanner\n"), 25, "");
    CHECK(parse(noName, g, err));
    CHECK(g.name == "default");
}

static void testErrors() {
    ScannerGeometry g;
    std::string err;

    // CT file
    CHECK(!parse("modality: CT\nscanner name: x\n", g, err));
    CHECK(err.find("PET") != std::string::npos);

    // Missing required field
    std::string missing = kBase;
    missing.replace(missing.find("number of rsectors: 32, 32\n"), 27, "");
    CHECK(!parse(missing, g, err));
    CHECK(err.find("number of rsectors") != std::string::npos);

    // scanner radius must list one value per layer
    std::string oneRadius = kBase;
    oneRadius.replace(oneRadius.find("scanner radius: 321.3, 326.3"), 28, "scanner radius: 321.3");
    CHECK(!parse(oneRadius, g, err));
    CHECK(err.find("scanner radius") != std::string::npos);

    // A count that differs between layers is not supported
    std::string differ = kBase;
    differ.replace(differ.find("number of rsectors: 32, 32"), 26, "number of rsectors: 32, 16");
    CHECK(!parse(differ, g, err));
    CHECK(err.find("differs between layers") != std::string::npos);

    // Non numeric value
    std::string bad = kBase;
    bad.replace(bad.find("crystals size depth: 5, 5"), 25, "crystals size depth: five");
    CHECK(!parse(bad, g, err));

    // Negative gap
    CHECK(!parse(std::string(kBase) + "crystal gap axial: -1, -1\n", g, err));

    // Unreadable file
    CHECK(!ParseCastorGeom("/nonexistent/file.geom", g, err));
}

static void testWarnings() {
    ScannerGeometry g;
    std::string err;
    CHECK(parse(std::string(kBase) + "number of modules transaxial: 2, 2\n", g, err));
    CHECK(!g.warnings.empty());
}

int main() {
    testShippedScanners();
    testPositions32x32x2();
    testGapsAreEdgeToEdge();
    testPerLayerTransaxialCounts();
    testParsingDetails();
    testErrors();
    testWarnings();
    if (failures == 0) std::printf("All scanner geometry tests passed\n");
    return failures == 0 ? 0 : 1;
}
