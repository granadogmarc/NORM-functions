#ifndef NORMFUNCTIONS_H
#define NORMFUNCTIONS_H

#include <iostream>
#include <map>
#include <array>
#include <cmath>
#include <TVector3.h>
#include <TTree.h>
#include <TArrayF.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "TROOT.h"
#include "TApplication.h"
#include "TGClient.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include <numeric>
#include <tuple>
#include <cstdint>
#include <unordered_map>
#include <utility>

typedef float Float_t;

struct Sinogram{

	float R,S,Phi,Cos_theta,x1,y1,z1,x2,y2,z2;

	Sinogram()
	: R(0),S(0),Phi(0),Cos_theta(1),x1(0),y1(0),z1(0),x2(0),y2(0),z2(0) {}

};




struct crystalKey {
    int rsectorID1,rsectorID2;
    int moduleID1,moduleID2;

    bool operator==(const crystalKey& other) const {
        return rsectorID1 == other.rsectorID1 &&
        		rsectorID2 == other.rsectorID2 &&
				moduleID1 == other.moduleID1 &&
				moduleID2 == other.moduleID2;
    }
};

struct crystalKeyHash {
    std::size_t operator()(const crystalKey& k) const {
        // Simple hash combine trick using shifts and XOR
        std::size_t h1 = std::hash<int>()(k.rsectorID1);
        std::size_t h2 = std::hash<int>()(k.rsectorID2);
        std::size_t h3 = std::hash<int>()(k.moduleID1);
        std::size_t h4 = std::hash<int>()(k.moduleID2);


        // Combine them (borrowed from boost::hash_combine)
        std::size_t seed = h1;
        seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= h4 + 0x9e3779b9 + (seed << 6) + (seed >> 2);

        return seed;
    }
};


struct SFKey {
    int layer1, layer2;
    int crystal1, crystal2;
    int submodule1, submodule2;
    int rSecDiff;

    uint64_t toIndex(uint32_t nLayers,
                     uint32_t nCrystalsAxial,
                     uint32_t nCrystalsTransaxial,
                     uint32_t nSubmodulesAxial,
                     uint32_t nSubmodulesTransaxial,
                     uint32_t nRsectorsAxial,
                     uint32_t nRsectorsAngPos) const

    {
        uint64_t maxLayer      = nLayers;
        uint64_t maxCrystal    = uint64_t(nCrystalsAxial)     * nCrystalsTransaxial;
        uint64_t maxSubmodule  = uint64_t(nSubmodulesAxial)   * nSubmodulesTransaxial;
        uint64_t maxModuleDiff = uint64_t(nRsectorsAxial)      * nRsectorsAngPos;

        uint64_t stride_layer      = maxLayer;
        uint64_t stride_layer2     = stride_layer      * maxLayer;
        uint64_t stride_crystal1   = stride_layer2;
        uint64_t stride_crystal2   = stride_crystal1   * maxCrystal;
        uint64_t stride_submodule1 = stride_crystal2   * maxCrystal;
        uint64_t stride_submodule2 = stride_submodule1 * maxSubmodule;
        uint64_t stride_rSecDiff = stride_submodule2 * maxSubmodule;

        return
            uint64_t(layer1)
          + uint64_t(layer2)      * stride_layer
          + uint64_t(crystal1)    * stride_layer2
          + uint64_t(crystal2)    * stride_crystal2
          + uint64_t(submodule1)  * stride_submodule1
          + uint64_t(submodule2)  * stride_submodule2
          + uint64_t(rSecDiff)  * stride_rSecDiff;
    }
};


struct DetectorIndices {
    uint32_t rsectorID;
    uint32_t moduleID;
    uint32_t submoduleID;
    uint32_t crystalID;
    uint32_t layerID;
};




std::vector<DetectorIndices> buildCastorIDLUT(
    uint32_t nRsectorsAngPos,
    uint32_t nRsectorsAxial,
    bool     invertDetOrder,
    int      rsectorIdOrder,
    uint32_t nModulesTransaxial,
    uint32_t nModulesAxial,
    uint32_t nSubmodulesTransaxial,
    uint32_t nSubmodulesAxial,
    uint32_t nCrystalsTransaxial,
    uint32_t nCrystalsAxial,
    uint8_t  nLayers,
    uint32_t *nCrystalPerLayer,
    uint32_t nLayersRptTransaxial,
    uint32_t nLayersRptAxial
);

struct SymmetryCounters{
    int totalCount = 0;
    int contributingCastorIDs = 0;
    double meanIntegral = 0;
    double meanR = 0;
};

//Change the map
//using crysLORMap = std::unordered_map<crystalKey, SymmetryCounters, crystalKeyHash>;
//into a matrix
//using matrixCrystalLOR = std::vector<std::vector<SymmetryCounters>>;  // <- your full LUT
//using matrixAxialGeometricFact =  std::vector<std::vector<SymmetryCounters>>;  // <- your full LUT
using matrixRingsComponent 		=	std::vector<std::vector<double>>;


//Change the map
//using SFMap = std::unordered_map<SFKey, SymmetryCounters, SFKeyHash>;
//into a matrix
//using matrixSymFactors = std::vector<SymmetryCounters>;  // <- your full LUT
//using matrixRadialComponent = std::vector<SymmetryCounters>;  // <- your full LUT
using vectorRingComponent			= std::vector<double>;
using vectorRadialComponent	= std::vector<double>;


struct LORValues {
    int numHits = 0;
    float integral = 0;
    float integral2 = 0;
    Sinogram sinogram{};
};

enum SimulType {UnknownSimul, EmissionSimul, TransmissionSimul,
                NormalizationSimul, AttenuationSimul, CTSimul,
                UniformSimul,LabelSimul};

//No need for key Type, we change it to a vector:
//using KeyType = std::pair<uint32_t, uint32_t>; // ID1, ID2

//using matrixLOR = std::vector<std::vector<LORValues>>;  // <- your full LUT


struct Phantom {
    std::array<float, 3> center;       // 3D or 4D position
    std::array<float, 3> half_axis;    // Half-lengths along axes
    float linf, lsup;                  // Truncation along the z-axis
    float theta, phi;                   // Orientation angles
    float ct, st, cp, sp;               // Precomputed cosines & sines
    float em_value, em_slope;           // Emission properties
    float em_polar, em_azimut;          // Emission gradient direction
    float mu_value, ct_value;           // Attenuation/CT coefficient
    int label_value, em_radial;         // Classification and emission order
    float zmin, zmax;                   // Spatial boundaries
    std::string name;                    // Object name
    float trunc_polar, trunc_azimut;    // Additional truncation angles
    float trunc_min, trunc_max;         // Additional truncation bounds

    // Constructor for easy initialization
    Phantom()
        : center{0.0f, 0.0f, 0.0f},
          half_axis{0.0f, 0.0f, 0.0f},
          linf(0), lsup(0), theta(0), phi(0),
          ct(0), st(0), cp(0), sp(0),
          em_value(0), em_slope(0), em_polar(0), em_azimut(0),
          mu_value(0), ct_value(0),
          label_value(0), em_radial(0),
          zmin(0), zmax(0),
          name(""),
          trunc_polar(0), trunc_azimut(0), trunc_min(0), trunc_max(0) {}
};
/*
struct RadialNormCounter {
    std::vector<Long64_t> radialCounter;
    int nBins;
    double rMax;
    double binSize;

    RadialNormCounter(int numBins, double radMax)
        : nBins(numBins),
          radialCounter(numBins, 0),
          rMax(radMax),
          binSize(radMax / numBins) // ✅ keep as double
    {
        std::cout << "We've created a new radialNormCounter with nBins = "
                  << nBins << " and binSize = " << binSize << std::endl;
    }

    void addRadialCount(double r)
    {
        if (std::abs(r) < rMax) {
            int binValue = static_cast<int>(r / binSize);
            if (binValue >= 0 && binValue < static_cast<int>(radialCounter.size())) {
                radialCounter[binValue]++;
            } else {
                std::cerr << "Warning: invalid binValue " << binValue << " for r=" << r << "\n";
            }
        }
    }

    Long64_t getRadialCount(double r) const {
        if (std::abs(r) < rMax) {
            int binValue = static_cast<int>(r / binSize);
            if (binValue >= 0 && binValue < static_cast<int>(radialCounter.size())) {
                return radialCounter[binValue];
            } else {
                std::cerr << "Warning: invalid binValue " << binValue << " for r=" << r << "\n";
                return 1;
            }
        } else {
            return 1;
        }
    }

	Long64_t getRadialCountInt(int r) const {
		    if (std::abs(r) < rMax) {
		        int binValue = r;

		        // Safety check: prevent out-of-range access
		        if (binValue >= 0 && binValue < static_cast<int>(radialCounter.size())) {
		            return radialCounter[binValue];
		        } else {
		        	std::cerr << "Warning: invalid binValue " << binValue << " for r=" << r << "\n";
		            return 1;  // Bin index is invalid, return default value
		        }
		    } else {
		        return 1;  // r is out of the allowed range
		    }
		}


    Long64_t getMaxRadialCount() const {
        return radialCounter.empty() ? 0 :
               *std::max_element(radialCounter.begin(), radialCounter.end());
    }
};

struct RadialNormCounter{

	std::vector<Long64_t> radialCounter;
	int					nBins;
	double 				rMax;
	double				binSize;
	//Constructor
	RadialNormCounter(int numBins,double radMax)
	: 	nBins(numBins),
		radialCounter(numBins,0),
	  rMax(radMax),
	  binSize(radMax/numBins){
		std::cout<<"We've created a new radialNormCounter with nBins = "<<nBins<<" and binSize = "<<binSize<<std::endl;
	}

	//Increment method
	void addRadialCount(double r)
	{

		if (abs(r)<rMax){
			int binValue = static_cast<int>(r / binSize);
			radialCounter[binValue]++;
		}
		else std::cout<<"Radius is larger than detector radius!"<<std::endl;
	}
	//Getters
	Long64_t getRadialCount(double r) const {
	    if (std::abs(r) < rMax) {
	        int binValue = static_cast<int>(r / binSize);

	        // Safety check: prevent out-of-range access
	        if (binValue >= 0 && binValue < static_cast<int>(radialCounter.size())) {
	            return radialCounter[binValue];
	        } else {
	        	std::cerr << "Warning: invalid binValue " << binValue << " for r=" << r << "\n";
	            return 1;  // Bin index is invalid, return default value
	        }
	    } else {
	        return -999;  // r is out of the allowed range
	    }
	}

	Long64_t getRadialCountInt(int r) const {
		    if (std::abs(r) < rMax) {
		        int binValue = r;

		        // Safety check: prevent out-of-range access
		        if (binValue >= 0 && binValue < static_cast<int>(radialCounter.size())) {
		            return radialCounter[binValue];
		        } else {
		        	std::cerr << "Warning: invalid binValue " << binValue << " for r=" << r << "\n";
		            return 1;  // Bin index is invalid, return default value
		        }
		    } else {
		        return -999;  // r is out of the allowed range
		    }
		}




    Long64_t getMaxRadialCount() const {
        return radialCounter.empty() ? 0 :
               *std::max_element(radialCounter.begin(), radialCounter.end());
    }
};
*/

struct DetectorCounters {
    std::vector<Long64_t> countsSubmodules;
    std::vector<Long64_t> countsCrystals;
    std::vector<Long64_t> countsLayers;

    // Constructor
    DetectorCounters(int nSubModules, int nCrystals, int nLayers)
        : countsSubmodules(nSubModules, 0),
          countsCrystals(nCrystals, 0),
          countsLayers(nLayers, 0) {}

    // Increment methods
    void addSubmoduleCount(int submoduleID) {
        if (submoduleID >= 0 && submoduleID < (int)countsSubmodules.size()) {
            countsSubmodules[submoduleID]++;
        }
    }

    void addCrystalCount(int crystalID) {
        if (crystalID >= 0 && crystalID < (int)countsCrystals.size()) {
            countsCrystals[crystalID]++;
        }
    }

    void addLayerCount(int layerID) {
        if (layerID >= 0 && layerID < (int)countsLayers.size()) {
            countsLayers[layerID]++;
        }
    }

    // Getters
    Long64_t getSubmoduleCount(int submoduleID) const {
        return countsSubmodules[submoduleID];
    }

    Long64_t getCrystalCount(int crystalID) const {
        return countsCrystals[crystalID];
    }

    Long64_t getLayerCount(int layerID) const {
        return countsLayers[layerID];
    }


    // Max getters
    Long64_t getMaxSubmoduleCount() const {
        return countsSubmodules.empty() ? 0 :
               *std::max_element(countsSubmodules.begin(), countsSubmodules.end());
    }

    Long64_t getMaxCrystalCount() const {
        return countsCrystals.empty() ? 0 :
               *std::max_element(countsCrystals.begin(), countsCrystals.end());
    }

    Long64_t getMaxLayerCount() const {
        return countsLayers.empty() ? 0 :
               *std::max_element(countsLayers.begin(), countsLayers.end());
    }

};






// forward declarations
std::pair<size_t, size_t> processFile( const std::string &filename,
                  //matrixLOR          &castorIDMatrix,
				  //matrixSymFactors 			&symmetryMatrix,
				  //matrixCrystalLOR 		&crystalLORMatrix,
				  //matrixAxialGeometricFact &axialGeomFactMatrix,
				  matrixRingsComponent		&ringsComponentMatrix,
				  double					meanRingsComponentMatrix,
				  matrixRingsComponent		&blockTrAComponentMatrix,
				  double						meanBlockTrAComponentMatrix,
				  vectorRadialComponent		&radialComponentVector,
				  double					meanRingComponentVector,
				  vectorRingComponent		&ringComponentVector,
				  DetectorCounters &detectorEfficyCounts,
				  //RadialNormCounter	&radialNormCounts,
				  //double 			&meanCrystalLORCounts,
				  Long64_t			&usedLOR,
				  Long64_t			&totalEvents,
                  uint32_t          nRsectorsAngPos,
                  uint32_t          nRsectorsAxial,
                  bool              invertDetOrder,
                  int               rsectorIdOrder,
                  uint32_t          nModulesTransaxial,
                  uint32_t          nModulesAxial,
                  uint32_t          nSubmodulesTransaxial,
                  uint32_t          nSubmodulesAxial,
                  uint32_t          nCrystalsTransaxial,
                  uint32_t          nCrystalsAxial,
                  uint8_t           nLayers,
                  uint32_t*   		nCrystalPerLayer,
                  uint32_t          nLayersRptTransaxial,
                  uint32_t          nLayersRptAxial,
                  const Phantom    &myPhantom,
                  const Phantom    &emptyPhantom,
                                    float			crystalDepth,
                                    float             detectorRadius,
                                    float             transAxialSize,
                                    float             axialSize
                );



void computeNormalizationFactors( const std::vector<std::string> &filenames, const std::string& scannerName,
		const std::string& outputMatrixFileName,
		const std::string& outputDir,
        uint32_t nRsectorsAngPos,
        uint32_t nRsectorsAxial,
        bool     invertDetOrder,
        int      rsectorIdOrder,
        uint32_t nModulesTransaxial,
        uint32_t nModulesAxial,
        uint32_t nSubmodulesTransaxial,
        uint32_t nSubmodulesAxial,
        uint32_t nCrystalsTransaxial,
        uint32_t nCrystalsAxial,
        uint8_t  nLayers,
        uint32_t* nCrystalPerLayer,
        uint32_t nLayersRptTransaxial,
        uint32_t nLayersRptAxial,
        const Phantom &myPhantom,
        const Phantom &emptyPhantom,
        float transAxialSize,
        float axialSize,
		float	crystalDepth,
        float detectorRadius,
        const std::string &outCSV
      );


struct ScannerConfig {
    uint32_t nRsectorsAngPos = 0;
    uint32_t nRsectorsAxial = 0;
    bool invertDetOrder = false;
    int rsectorIdOrder = 0;
    int lyr_rpt_flag = 1;
    uint32_t nModulesTransaxial = 0;
    uint32_t nModulesAxial = 0;
    uint32_t nSubmodulesTransaxial = 0;
    uint32_t nSubmodulesAxial = 0;
    uint32_t nCrystalsTransaxial = 0;
    uint32_t nCrystalsAxial = 0;
    uint8_t nLayers = 0;
    std::vector<uint32_t> nCrystalPerLayer;
    uint32_t nLayersRptTransaxial = 1;
    uint32_t nLayersRptAxial = 1;
};

inline ScannerConfig ParseGeomFile(const std::string& filename) {
    ScannerConfig config;
    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string label;
        std::getline(iss, label, ':');
        std::string values_str;
        std::getline(iss, values_str);

        std::istringstream vss(values_str);
        std::vector<uint32_t> values;
        uint32_t val;
        while (vss >> val) values.push_back(val);

        if (label.find("number of layers") != std::string::npos) {
            config.nLayers = static_cast<uint8_t>(values[0]);
        } else if (label.find("number of rsectors") != std::string::npos && label.find("axial") == std::string::npos) {
            config.nRsectorsAngPos = values[0];
        } else if (label.find("number of rsectors axial") != std::string::npos) {
            config.nRsectorsAxial = values[0];
        } else if (label.find("number of modules transaxial") != std::string::npos) {
            config.nModulesTransaxial = values[0];
        } else if (label.find("number of modules axial") != std::string::npos) {
            config.nModulesAxial = values[0];
        } else if (label.find("number of submodules transaxial") != std::string::npos) {
            config.nSubmodulesTransaxial = values[0];
        } else if (label.find("number of submodules axial") != std::string::npos) {
            config.nSubmodulesAxial = values[0];
        } else if (label.find("number of crystals transaxial") != std::string::npos) {
            config.nCrystalsTransaxial = values[0];
        } else if (label.find("number of crystals axial") != std::string::npos) {
            config.nCrystalsAxial = values[0];
        }
    }

    // Compute nCrystalPerLayer
    uint32_t totalCrystalsPerLayer =
        config.nModulesTransaxial * config.nModulesAxial *
        config.nSubmodulesTransaxial * config.nSubmodulesAxial *
        config.nCrystalsTransaxial * config.nCrystalsAxial *
        config.nRsectorsAngPos;

    config.nCrystalPerLayer.assign(config.nLayers, totalCrystalsPerLayer);

    return config;
}



TVector3 applyRotation(const TVector3& pos, double angle) ;

TVector3 convertToPosition(double x0, double y0, double z0,
                           double deltaPhi,
                           int layerID, int crystalID,
                           int submoduleID, int moduleID,
                           int rsectorID,
                           float crystalDepth,
                           float transAxialSize,
                           float axialSize,
                           uint8_t nLayers,
                           uint32_t nCrystalsTransaxial,
                           uint32_t nSubmodulesAxial);

uint32_t ReducedID(uint32_t nRsectorsAngPos,
                          uint32_t nRsectorsAxial,
                          int rsectorIdOrder,
                          uint32_t nSubmodulesTransaxial,
                          uint32_t nSubmodulesAxial,
                          uint32_t nCrystalsTransaxial,
                          uint32_t nCrystalsAxial,
                          uint8_t nLayers,
                          int layerID,
                          int crystalID,
                          int submoduleID,
                          int rsectorID);

uint32_t ConvertIDcylindrical(uint32_t  nRsectorsAngPos,
                              uint32_t  nRsectorsAxial,
                              bool      a_invertDetOrder,
                              int       a_rsectorIdOrder,
                              uint32_t  nModulesTransaxial,
                              uint32_t  nModulesAxial,
                              uint32_t  nSubmodulesTransaxial,
                              uint32_t  nSubmodulesAxial,
                              uint32_t  nCrystalsTransaxial,
                              uint32_t  nCrystalsAxial,
                              uint8_t   nLayers,
                              uint32_t	*nCrystalPerLayer,
							  uint32_t	nLayersRptTransaxial,
							  uint32_t	nLayersRptAxial,
                              int32_t	layerID,
                              int32_t	crystalID,
                              int32_t	submoduleID,
                              int32_t	moduleID,
                              int32_t	rsectorID);


std::map<uint32_t, std::array<Float_t,3>> computeLUT(TTree* Hits);


std::pair<float, float> ComputeZBounds(TVector3& pos1, TVector3& pos2);//changed


float ComputeDistance(TVector3& a, TVector3& b);//changed

float ComputePhantomLineIntegral (TVector3& a,TVector3& b,
                                 const Phantom& phantom,
                                 int mat_data,
                                 Double_t ToF);//changed

float emRad(float a, float b, float c, float s, int rad_type);


float LineIntegral(TVector3& a, TVector3& b, Phantom& ob, int mat_data, float ToF);//changed

Sinogram ConvertToSinogram(TVector3& gPos1,TVector3& gPos2);//Changed

std::vector<std::string> expandWildcard(const std::string &pattern) ;

void writeNormEntryNormMatrixFile(std::ofstream& out, uint32_t c1, uint32_t c2, float n);

void writeCdhHeaderNormMatrixFile(const std::string& outputMatrixFileName,const std::string& outputDir,
                    const std::string& scannerName,
                    size_t nEvents);

double meanVector(const std::vector<double>& v);
double meanMatrix(const std::vector<std::vector<double>>& M);

double scaledMeanMatrix(const std::vector<std::vector<double>>& mat,const std::vector<double>& weights);
inline std::ofstream openCSV(const std::string& filename);
inline int mod(int x, int m);
#endif 
