#include "scannerGeometry.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <map>
#include <numeric>
#include <ostream>
#include <sstream>

namespace {

const double kPi = 3.14159265358979323846;

std::string trim(const std::string &s) {
    size_t b = 0, e = s.size();
    while (b < e && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
    return s.substr(b, e - b);
}

std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

bool parseDouble(const std::string &text, double &value) {
    const std::string s = trim(text);
    if (s.empty()) return false;
    char *end = nullptr;
    value = std::strtod(s.c_str(), &end);
    return end != s.c_str() && *end == '\0';
}

// "a, b, c" -> {a, b, c}
bool parseList(const std::string &text, std::vector<double> &out) {
    out.clear();
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, ',')) {
        double v = 0.0;
        if (!parseDouble(item, v)) return false;
        out.push_back(v);
    }
    return !out.empty();
}

typedef std::map<std::string, std::string> Fields;

// Collects "label: value" pairs. Text after '#' is a comment.
Fields collectFields(const std::string &text) {
    Fields fields;
    std::istringstream in(text);
    std::string line;
    while (std::getline(in, line)) {
        const size_t hash = line.find('#');
        if (hash != std::string::npos) line.erase(hash);
        const size_t colon = line.find(':');
        if (colon == std::string::npos) continue;
        const std::string label = toLower(trim(line.substr(0, colon)));
        const std::string value = trim(line.substr(colon + 1));
        if (!label.empty()) fields[label] = value;
    }
    return fields;
}

// Reads the values of a (possibly layer-dependent) field into one entry per layer.
//   perLayerRequired: the field must list exactly nLayers values
//   otherwise a single value is applied to all layers
bool readPerLayer(const Fields &fields, const std::string &label, uint32_t nLayers,
                  bool required, bool perLayerRequired, double defaultValue,
                  std::vector<double> &out, std::string &error) {
    const Fields::const_iterator it = fields.find(label);
    if (it == fields.end() || it->second.empty()) {
        if (required) {
            error = "missing required field '" + label + "'";
            return false;
        }
        out.assign(nLayers, defaultValue);
        return true;
    }
    std::vector<double> values;
    if (!parseList(it->second, values)) {
        error = "cannot read the numeric value(s) of '" + label + "': " + it->second;
        return false;
    }
    if (values.size() == nLayers) {
        out = values;
    } else if (values.size() == 1 && !perLayerRequired) {
        out.assign(nLayers, values[0]);
    } else {
        error = "field '" + label + "' has " + std::to_string(values.size()) +
                " value(s) but the scanner has " + std::to_string(nLayers) + " layer(s)";
        return false;
    }
    return true;
}

// This tool needs the same value for every layer for these fields.
bool sameForAllLayers(const std::vector<double> &values, const std::string &label,
                      double &out, std::string &error) {
    out = values.front();
    for (size_t i = 1; i < values.size(); ++i) {
        if (std::fabs(values[i] - out) > 1e-9 * std::max(1.0, std::fabs(out))) {
            error = "field '" + label + "' differs between layers; this tool requires the same "
                    "value for all layers";
            return false;
        }
    }
    return true;
}

bool toCount(double value, const std::string &label, uint32_t &out, std::string &error) {
    if (value < 1.0 || std::fabs(value - std::round(value)) > 1e-9) {
        error = "field '" + label + "' must be a positive integer";
        return false;
    }
    out = static_cast<uint32_t>(std::lround(value));
    return true;
}

} // namespace

bool ScannerGeometry::finalize(std::string &error) {
    warnings.clear();

    if (nLayers == 0 || layerRadius.size() != nLayers || layerDepth.size() != nLayers) {
        error = "inconsistent number of layers";
        return false;
    }
    if (nCrystalsTransaxial.size() != nLayers || crystalSizeTrans.size() != nLayers ||
        crystalGapTrans.size() != nLayers) {
        error = "inconsistent number of layers";
        return false;
    }
    if (nRsectorsAngPos == 0 || nCrystalsAxial == 0 ||
        nModulesTransaxial == 0 || nModulesAxial == 0 ||
        nSubmodulesTransaxial == 0 || nSubmodulesAxial == 0 || nRsectorsAxial == 0) {
        error = "all element counts must be at least 1";
        return false;
    }
    if (crystalSizeAxial <= 0.0) {
        error = "'crystals size axial' must be positive";
        return false;
    }
    for (uint32_t l = 0; l < nLayers; ++l) {
        if (nCrystalsTransaxial[l] == 0) {
            error = "all element counts must be at least 1";
            return false;
        }
        if (crystalSizeTrans[l] <= 0.0) {
            error = "'crystals size trans' must be positive";
            return false;
        }
        if (crystalGapTrans[l] < 0) {
            error = "gaps must not be negative";
            return false;
        }
        if (layerRadius[l] <= 0.0 || layerDepth[l] <= 0.0) {
            error = "'scanner radius' and 'crystals size depth' must be positive";
            return false;
        }
    }
    if (crystalGapAxial < 0 || submoduleGapTrans < 0 ||
        submoduleGapAxial < 0 || moduleGapTrans < 0 || moduleGapAxial < 0 || rsectorGapAxial < 0) {
        error = "gaps must not be negative";
        return false;
    }
    if (rsectorAngularSpanDeg <= 0.0 || rsectorAngularSpanDeg > 360.0 + 1e-9) {
        error = "'rsectors angular span' must be in ]0, 360]";
        return false;
    }

    // Edge to edge gaps: an element is as large as its children plus the gaps between them.
    crystalPitchAxial = crystalSizeAxial + crystalGapAxial;

    const double submoduleSizeAxial = nCrystalsAxial * crystalPitchAxial - crystalGapAxial;
    submodulePitchAxial = submoduleSizeAxial + submoduleGapAxial;

    moduleSizeAxial = nSubmodulesAxial * submodulePitchAxial - submoduleGapAxial;
    modulePitchAxial = moduleSizeAxial + moduleGapAxial;

    crystalPitchTrans.resize(nLayers);
    moduleSizeTrans.resize(nLayers);
    for (uint32_t l = 0; l < nLayers; ++l) {
        crystalPitchTrans[l] = crystalSizeTrans[l] + crystalGapTrans[l];
        moduleSizeTrans[l] = nCrystalsTransaxial[l] * crystalPitchTrans[l] - crystalGapTrans[l];
    }

    if (!buildTransaxialGrid(error)) return false;

    rsectorStepRad = rsectorAngularSpanDeg * kPi / 180.0 / static_cast<double>(nRsectorsAngPos);

    layerCentreRadius.resize(nLayers);
    for (uint32_t l = 0; l < nLayers; ++l) layerCentreRadius[l] = layerRadius[l] + 0.5 * layerDepth[l];

    cosRsector.resize(nRsectorsAngPos);
    sinRsector.resize(nRsectorsAngPos);
    for (uint32_t r = 0; r < nRsectorsAngPos; ++r) {
        const double a = rsectorFirstAngleDeg * kPi / 180.0 + r * rsectorStepRad;
        cosRsector[r] = std::cos(a);
        sinRsector[r] = std::sin(a);
    }

    if (nModulesTransaxial > 1 || nSubmodulesTransaxial > 1 || nCrystalsAxial > 1 || nRsectorsAxial > 1) {
        warnings.push_back(
            "crystal positions only use crystalID (transaxial), submoduleID/moduleID (axial) and "
            "rsectorID (angular); the transaxial module/submodule counts, the axial crystal count "
            "and the axial rsector count are not used for positions (all must be 1 for exact positions)");
    }
    if (gridMismatch > 0.25) {
        std::ostringstream w;
        w << "the crystals of the layers do not span the same transaxial width: crystal "
             "positions on the common transaxial grid are off by up to "
          << gridMismatch << " grid spacing(s)";
        warnings.push_back(w.str());
    }
    if (rsectorAngularSpanDeg < 360.0 - 1e-9) {
        warnings.push_back("rsectors angular span < 360: the angle between rsectors is assumed to be "
                           "span / number of rsectors");
    }
    return true;
}

bool ScannerGeometry::buildTransaxialGrid(std::string &error) {
    uint64_t lcm = 1;
    for (uint32_t l = 0; l < nLayers; ++l) {
        lcm = std::lcm(lcm, static_cast<uint64_t>(nCrystalsTransaxial[l]));
        if (lcm > (1u << 20)) {
            error = "the transaxial crystal counts of the layers have no usable common grid "
                    "(least common multiple too large)";
            return false;
        }
    }
    const uint32_t L = static_cast<uint32_t>(lcm);

    uniformCrystalsTransaxial = true;
    bool allRatiosOdd = true;
    for (uint32_t l = 0; l < nLayers; ++l) {
        if (nCrystalsTransaxial[l] != nCrystalsTransaxial[0]) uniformCrystalsTransaxial = false;
        if ((L / nCrystalsTransaxial[l]) % 2 == 0) allRatiosOdd = false;
    }

    gridStep.resize(nLayers);
    gridOffset.resize(nLayers);
    gridPointsPerSubmodule = allRatiosOdd ? L : 2 * L;
    for (uint32_t l = 0; l < nLayers; ++l) {
        const uint32_t r = L / nCrystalsTransaxial[l];
        if (allRatiosOdd) {
            // centre of crystal j = centre of point j*r + (r-1)/2
            gridStep[l] = r;
            gridOffset[l] = (r - 1) / 2;
        } else {
            // centre of crystal j = point (2j+1)*r, points on the half pitch of the finest layer
            gridStep[l] = 2 * r;
            gridOffset[l] = r;
        }
    }

    // The grid assumes the crystals of every layer share the same transaxial
    // width (count * pitch). A difference d between two layers moves the outer
    // crystals by up to d/2 with respect to each other.
    gridMismatch = 0.0;
    if (!uniformCrystalsTransaxial) {
        double minSpan = nCrystalsTransaxial[0] * crystalPitchTrans[0];
        double maxSpan = minSpan, sumSpan = 0.0;
        for (uint32_t l = 0; l < nLayers; ++l) {
            const double span = nCrystalsTransaxial[l] * crystalPitchTrans[l];
            minSpan = std::min(minSpan, span);
            maxSpan = std::max(maxSpan, span);
            sumSpan += span;
        }
        const double spacing = sumSpan / nLayers / gridPointsPerSubmodule;
        gridMismatch = 0.5 * (maxSpan - minSpan) / spacing;
    }
    return true;
}

uint32_t ScannerGeometry::maxNCrystalsTransaxial() const {
    return nCrystalsTransaxial.empty()
               ? 0
               : *std::max_element(nCrystalsTransaxial.begin(), nCrystalsTransaxial.end());
}

uint32_t ScannerGeometry::nCrystalsInLayer(uint32_t layerID) const {
    return nRsectorsAngPos * nRsectorsAxial * nModulesTransaxial * nModulesAxial *
           nSubmodulesTransaxial * nSubmodulesAxial * nCrystalsTransaxial[layerID] * nCrystalsAxial;
}

Position3 ScannerGeometry::crystalPosition(int layerID, int crystalID, int submoduleID,
                                           int moduleID, int rsectorID) const {
    // Position in the frame of rsector 0: x radial (layer centre), y transaxial, z axial.
    const double xl = layerCentreRadius[layerID];
    const double yl = (crystalID - 0.5 * (static_cast<double>(nCrystalsTransaxial[layerID]) - 1.0)) *
                      crystalPitchTrans[layerID];
    const double zl = (moduleID - 0.5 * (static_cast<double>(nModulesAxial) - 1.0)) * modulePitchAxial +
                      (submoduleID - 0.5 * (static_cast<double>(nSubmodulesAxial) - 1.0)) * submodulePitchAxial;

    double c, s;
    if (rsectorID >= 0 && static_cast<size_t>(rsectorID) < cosRsector.size()) {
        c = cosRsector[rsectorID];
        s = sinRsector[rsectorID];
    } else {
        const double a = rsectorFirstAngleDeg * kPi / 180.0 + rsectorID * rsectorStepRad;
        c = std::cos(a);
        s = std::sin(a);
    }

    Position3 p;
    p.x = xl * c - yl * s;
    p.y = xl * s + yl * c;
    p.z = zl;
    return p;
}

void ScannerGeometry::print(std::ostream &os) const {
    os << "Scanner geometry (" << name << ")\n"
       << "  layers: " << nLayers << ", rsectors: " << nRsectorsAngPos
       << ", modules (trans x axial): " << nModulesTransaxial << " x " << nModulesAxial
       << ", submodules: " << nSubmodulesTransaxial << " x " << nSubmodulesAxial
       << ", crystals axial: " << nCrystalsAxial << "\n";
    os << "  crystal size axial: " << crystalSizeAxial
       << " mm, gaps crystal/submodule/module axial: " << crystalGapAxial << "/"
       << submoduleGapAxial << "/" << moduleGapAxial << " mm\n";
    for (uint32_t l = 0; l < nLayers; ++l) {
        os << "  layer " << l << ": " << nCrystalsTransaxial[l] << " crystals transaxial of "
           << crystalSizeTrans[l] << " mm (gap " << crystalGapTrans[l] << " mm, rsector block "
           << moduleSizeTrans[l] << " mm), scanner radius " << layerRadius[l] << " mm, depth "
           << layerDepth[l] << " mm (centre at " << layerCentreRadius[l] << " mm)\n";
    }
    if (!uniformCrystalsTransaxial) {
        os << "  common transaxial grid: " << gridPointsPerSubmodule << " points per submodule\n";
    }
    os << "  module axial size " << moduleSizeAxial << " mm, module pitch " << modulePitchAxial << " mm\n";
    os << "  rsector step: " << rsectorStepRad * 180.0 / kPi << " deg, first angle "
       << rsectorFirstAngleDeg << " deg\n";
    for (size_t i = 0; i < warnings.size(); ++i) os << "  WARNING: " << warnings[i] << "\n";
}

bool ParseCastorGeomText(const std::string &text, const std::string &defaultName,
                         ScannerGeometry &geom, std::string &error) {
    geom = ScannerGeometry();
    const Fields fields = collectFields(text);

    // Modality (PET only)
    const Fields::const_iterator mod = fields.find("modality");
    if (mod != fields.end() && toLower(mod->second) != "pet") {
        error = "modality is '" + mod->second + "' but this tool only supports PET scanners";
        return false;
    }

    const Fields::const_iterator nameIt = fields.find("scanner name");
    geom.name = (nameIt != fields.end() && !nameIt->second.empty()) ? nameIt->second : defaultName;

    // Number of layers
    {
        std::vector<double> v;
        const Fields::const_iterator it = fields.find("number of layers");
        if (it == fields.end()) { error = "missing required field 'number of layers'"; return false; }
        if (!parseList(it->second, v) || v.size() != 1 || !toCount(v[0], "number of layers", geom.nLayers, error)) {
            if (error.empty()) error = "cannot read 'number of layers'";
            return false;
        }
    }
    const uint32_t nL = geom.nLayers;

    // Reads a field that must be identical for all layers and converts it
    auto scalar = [&](const std::string &label, bool required, double def, double &out) -> bool {
        std::vector<double> v;
        return readPerLayer(fields, label, nL, required, false, def, v, error) &&
               sameForAllLayers(v, label, out, error);
    };
    auto count = [&](const std::string &label, bool required, uint32_t def, uint32_t &out) -> bool {
        double d = 0.0;
        if (!scalar(label, required, def, d)) return false;
        return toCount(d, label, out, error);
    };

    if (!count("number of rsectors", true, 1, geom.nRsectorsAngPos)) return false;
    if (!count("number of rsectors axial", false, 1, geom.nRsectorsAxial)) return false;
    if (!count("number of modules transaxial", false, 1, geom.nModulesTransaxial)) return false;
    if (!count("number of modules axial", false, 1, geom.nModulesAxial)) return false;
    if (!count("number of submodules transaxial", false, 1, geom.nSubmodulesTransaxial)) return false;
    if (!count("number of submodules axial", false, 1, geom.nSubmodulesAxial)) return false;
    if (!count("number of crystals axial", true, 1, geom.nCrystalsAxial)) return false;
    {
        std::vector<double> v;
        if (!readPerLayer(fields, "number of crystals transaxial", nL, true, false, 1.0, v, error)) return false;
        geom.nCrystalsTransaxial.resize(nL);
        for (uint32_t l = 0; l < nL; ++l)
            if (!toCount(v[l], "number of crystals transaxial", geom.nCrystalsTransaxial[l], error)) return false;
    }

    // CASToR treats the crystal sizes as optional, but they define the positions here.
    if (!readPerLayer(fields, "crystals size trans", nL, true, false, 0.0, geom.crystalSizeTrans, error)) return false;
    if (!scalar("crystals size axial", true, 0.0, geom.crystalSizeAxial)) return false;

    if (!readPerLayer(fields, "crystal gap transaxial", nL, false, false, 0.0, geom.crystalGapTrans, error)) return false;
    if (!scalar("crystal gap axial", false, 0.0, geom.crystalGapAxial)) return false;
    if (!scalar("submodule gap transaxial", false, 0.0, geom.submoduleGapTrans)) return false;
    if (!scalar("submodule gap axial", false, 0.0, geom.submoduleGapAxial)) return false;
    if (!scalar("module gap transaxial", false, 0.0, geom.moduleGapTrans)) return false;
    if (!scalar("module gap axial", false, 0.0, geom.moduleGapAxial)) return false;
    if (!scalar("rsector gap axial", false, 0.0, geom.rsectorGapAxial)) return false;

    // Layer-dependent values without a meaningful common value: one value per layer is required.
    if (!readPerLayer(fields, "scanner radius", nL, true, true, 0.0, geom.layerRadius, error)) return false;
    if (!readPerLayer(fields, "crystals size depth", nL, true, true, 0.0, geom.layerDepth, error)) return false;

    // Rsector placement
    if (!scalar("rsectors first angle", false, 0.0, geom.rsectorFirstAngleDeg)) return false;
    if (!scalar("rsectors angular span", false, 360.0, geom.rsectorAngularSpanDeg)) return false;

    error.clear();
    return geom.finalize(error);
}

bool ParseCastorGeom(const std::string &path, ScannerGeometry &geom, std::string &error) {
    std::ifstream in(path.c_str());
    if (!in) {
        error = "cannot open scanner file: " + path;
        return false;
    }
    std::stringstream buffer;
    buffer << in.rdbuf();

    // Default name: file name without directory and extension
    std::string base = path;
    const size_t slash = base.find_last_of("/\\");
    if (slash != std::string::npos) base = base.substr(slash + 1);
    const size_t dot = base.find_last_of('.');
    if (dot != std::string::npos && dot > 0) base = base.substr(0, dot);

    return ParseCastorGeomText(buffer.str(), base, geom, error);
}
