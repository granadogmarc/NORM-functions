#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <random>
#include <omp.h>
#include "TROOT.h"
#include "TApplication.h"
#include "TGClient.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include <map>
#include <numeric>
#include "TLegend.h"
#include "normFunctions.h"

#include <filesystem>
#include <vector>
#include <string>
#include <glob.h>

void printUsage(const char* programName) {
    std::cerr << "Usage: " << programName << " [OPTIONS]\n\n"
              << "PET normalization factors computation utility.\n\n"
              << "Required arguments:\n"
              << "  -s, --system <name>      Scanner system name. Available systems:\n"
              << "                             CM2L_1ring_system\n"
              << "                             16x16x2_1ring_system\n"
              << "                             16x16x2_4rings_system\n"
              << "                             32x16x2_4rings_system\n"
	      << "                             32x32x2_4rings_system\n"
              << "  -i, --input <pattern>    Input ROOT file path or glob pattern\n"
              << "                             (use quotes for wildcards: 'path/*.root')\n"
              << "  -o, --outputFile <name>  Output file name (without extension)\n\n"
              << "Optional arguments:\n"
              << "  -d, --outputDir <path>   Output directory (default: current directory)\n"
              << "  -j, --threads <N>        Number of OpenMP threads to use\n"
              << "                             (default: all available cores)\n"
              << "  -as, --axial-sigma <value>\n"
              << "                           Gaussian smoothing sigma for axial normalization\n"
              << "                             (default: 0 = disabled)\n"
              << "                             Recommended: 0.6 to reduce sawtooth ripple\n"
              << "  -ts, --transaxial-sigma <value>\n"
              << "                           Gaussian smoothing sigma for transaxial normalization\n"
              << "                             (default: 0 = disabled)\n"
              << "                             Recommended: 1.0 for central LOR smoothing\n"
              << "  -h, --help               Show this help message and exit\n\n"
              << "Examples:\n"
              << "  " << programName << " -s CM2L_1ring_system -i 'data/*.root' -o norm_output\n"
              << "  " << programName << " -s 16x16x2_4rings_system -i data.root -o out -j 4\n";
}


int main(int argc,char**argv) {

	std::string scannerName;
	std::string pattern;
	std::string outputMatrixFileName;
	std::string outputDir;
	int numThreads = 0;  // 0 means use default (all available)
	double axialSigma = 0.0;       // Gaussian smoothing sigma for axial normalization (0 = disabled by default)
	double transaxialSigma = 0.0;  // Gaussian smoothing sigma for transaxial normalization (0 = disabled by default)

	if (argc == 1) {
		printUsage(argv[0]);
		return 1;
	}

	for (int i = 1; i < argc; ++i) {
		std::string arg = argv[i];

		if (arg == "-h" || arg == "--help") {
			printUsage(argv[0]);
			return 0;
		}
		else if (arg == "-s" || arg == "--system") {
			if (i + 1 < argc) {
				scannerName = argv[++i];
			} else {
				std::cerr << "Error: missing argument after " << arg << "\n";
				return 1;
			}
		}
		else if (arg == "-i" || arg == "--input") {
			if (i + 1 < argc) {
				pattern = argv[++i];
			} else {
				std::cerr << "Error: missing argument after " << arg << "\n";
				return 1;
			}
		}
		else if (arg == "-d" || arg == "--outputDir") {
			if (i + 1 < argc) {
				outputDir = argv[++i];
			} else {
				std::cerr << "Error: missing argument after " << arg << "\n";
				return 1;
			}
		}
		else if (arg == "-o" || arg == "--outputFile") {
			if (i + 1 < argc) {
				outputMatrixFileName = argv[++i];
			} else {
				std::cerr << "Error: missing argument after " << arg << "\n";
				return 1;
			}
		}
		else if (arg == "-j" || arg == "--threads") {
			if (i + 1 < argc) {
				numThreads = std::atoi(argv[++i]);
				if (numThreads <= 0) {
					std::cerr << "Error: invalid thread count. Must be a positive integer.\n";
					return 1;
				}
			} else {
				std::cerr << "Error: missing argument after " << arg << "\n";
				return 1;
			}
		}
		else if (arg == "-as" || arg == "--axial-sigma") {
			if (i + 1 < argc) {
				axialSigma = std::atof(argv[++i]);
				if (axialSigma < 0) {
					std::cerr << "Error: axial sigma must be non-negative.\n";
					return 1;
				}
			} else {
				std::cerr << "Error: missing argument after " << arg << "\n";
				return 1;
			}
		}
		else if (arg == "-ts" || arg == "--transaxial-sigma") {
			if (i + 1 < argc) {
				transaxialSigma = std::atof(argv[++i]);
				if (transaxialSigma < 0) {
					std::cerr << "Error: transaxial sigma must be non-negative.\n";
					return 1;
				}
			} else {
				std::cerr << "Error: missing argument after " << arg << "\n";
				return 1;
			}
		}
		else {
			std::cerr << "Error: unknown argument '" << arg << "'\n\n";
			printUsage(argv[0]);
			return 1;
		}
	}

	// Set OpenMP thread count if specified
	if (numThreads > 0) {
		omp_set_num_threads(numThreads);
		std::cout << "OpenMP: Using " << numThreads << " thread(s)\n";
	} else {
		std::cout << "OpenMP: Using " << omp_get_max_threads() << " thread(s) (default)\n";
	}


	std::cout<<"outDir = "<<outputDir<<" fileName "<<outputMatrixFileName<<std::endl;
  	  std::cout<<"pattern is "<< pattern<<std::endl;
      std::vector<std::string> files = expandWildcard(pattern);

      if (files.empty()) {
          std::cerr << "No files matched pattern: " << pattern << "\n";
          return 1;
      }

      std::cout << "Matched filenames:\n";
      for (size_t i = 0; i < files.size(); ++i) {
          std::cout << "  [" << i << "] " << files[i] << "\n";
      }
      std::cout << "Total number of files: " << files.size() << "\n";


  // 2) Precompute how many unique castorIDs we expect:


  bool invertDetOrder = false;
  int rsectorIdOrder = 0;


  uint32_t nRsectorsAngPos, nRsectorsAxial;

  uint32_t nModulesTransaxial, nModulesAxial;
  uint32_t nSubmodulesTransaxial, nSubmodulesAxial;
  uint32_t nCrystalsTransaxial, nCrystalsAxial;
  uint8_t nLayers;
  uint32_t nLayersRptTransaxial, nLayersRptAxial;
  double axialSize = 59; // in mm
  double transAxialSize = 59; // in mm
  double depth = 10; // in mm	
  // Define system parameters based on scannerName

 if (scannerName=="CM2L_1ring_system"){

  nRsectorsAngPos = 32;
  nRsectorsAxial = 1;
  nModulesTransaxial = 1;
  nModulesAxial = 1;
  nSubmodulesTransaxial = 1;
  nSubmodulesAxial = 32;
  nCrystalsTransaxial = 32;
  nCrystalsAxial = 1;
  nLayers = 2;
  nLayersRptTransaxial = 1;
  nLayersRptAxial = 1;
 }
 else if (scannerName =="16x16x2_1ring_system"){

  nRsectorsAngPos = 32;
  nRsectorsAxial = 1;
  nModulesTransaxial = 1;
  nModulesAxial = 1;
  nSubmodulesTransaxial = 1;
  nSubmodulesAxial = 16;
  nCrystalsTransaxial = 16;
  nCrystalsAxial = 1;
  nLayers = 2;
  nLayersRptTransaxial = 1;
  nLayersRptAxial = 1;
 }

 else if (scannerName =="16x16x2_4rings_system"){

  nRsectorsAngPos = 32;
  nRsectorsAxial = 1;
  nModulesTransaxial = 1;
  nModulesAxial = 4;
  nSubmodulesTransaxial = 1;
  nSubmodulesAxial = 16;
  nCrystalsTransaxial = 16;
  nCrystalsAxial = 1;
  nLayers = 2;
  nLayersRptTransaxial = 1;
  nLayersRptAxial = 1;
 }

else if (scannerName =="32x16x2_4rings_system"){

  nRsectorsAngPos = 32;
  nRsectorsAxial = 1;
  nModulesTransaxial = 1;
  nModulesAxial = 4;
  nSubmodulesTransaxial = 1;
  nSubmodulesAxial = 32;
  nCrystalsTransaxial = 16;
  nCrystalsAxial = 1;
  nLayers = 2;
  nLayersRptTransaxial = 1;
  nLayersRptAxial = 1;

 }


else if (scannerName =="32x32x2_4rings_system"){

  nRsectorsAngPos = 32;
  nRsectorsAxial = 1;
  nModulesTransaxial = 1;
  nModulesAxial = 4;
  nSubmodulesTransaxial = 1;
  nSubmodulesAxial = 32;
  nCrystalsTransaxial = 32;
  nCrystalsAxial = 1;
  nLayers = 2;
  nLayersRptTransaxial = 1;
  nLayersRptAxial = 1;

 }



 else{ std::cerr << "Error: no system provided from the expected list\n";
 return 1;}


 float		crystalDepth = 10. ;// in mm
 float		detectorRadius = 321.3 + crystalDepth/nLayers*0.5 + 0.5; //in mm //Here I add half a millimeter to ensure there is no drama with the projection caused by two different lenths

 uint32_t nCrystalPerLayer[nLayers] = {nRsectorsAngPos * nRsectorsAxial *nModulesTransaxial * nModulesAxial *nSubmodulesTransaxial * nSubmodulesAxial
		 	 	 	 	 	 	 	 	 * nCrystalsTransaxial * nCrystalsAxial *nLayersRptTransaxial,
										nRsectorsAngPos * nRsectorsAxial *nModulesTransaxial * nModulesAxial *nSubmodulesTransaxial * nSubmodulesAxial
										* nCrystalsTransaxial * nCrystalsAxial *nLayersRptTransaxial};//All of them


  Phantom myPhantom;

  	myPhantom.center =			{0.,0.,0.};
  	myPhantom.half_axis =		{310.,310.,30.};
  	myPhantom.linf =			-30.;
  	myPhantom.lsup =			30.;
  	myPhantom.trunc_max = 		30;
  	myPhantom.trunc_min =		-30;
  	myPhantom.theta = 			0.;      // in degrees or radians depending on your convention
  	myPhantom.phi =				0.;

  	myPhantom.ct =				std::cos(myPhantom.theta);
  	myPhantom.st =				std::sin(myPhantom.theta);
  	myPhantom.cp =				std::cos(myPhantom.phi);
  	myPhantom.sp =				std::sin(myPhantom.phi);
  	myPhantom.em_value =		1.; // 15.9xRowSize x ColSize/1000
  	myPhantom.em_slope =		0.;
  	myPhantom.em_polar =		0.0;
  	myPhantom.em_azimut =		0.0;


  	myPhantom.em_radial =		0.;
  	myPhantom.zmin =			-30.;
  	myPhantom.zmax =			30.;
  	myPhantom.name =			"Phantom Cylinder";

      Phantom emptyPhantom;

  	emptyPhantom.center =		{0.,0.,0.};
  	emptyPhantom.half_axis =	{300.,300.,30.};
  	emptyPhantom.linf =			-30.;
  	emptyPhantom.lsup =			30.;
  	emptyPhantom.trunc_max = 	30;
  	emptyPhantom.trunc_min =	-30;
  	emptyPhantom.theta = 		0.;      // in degrees or radians depending on your convention
  	emptyPhantom.phi =			0.;
  	emptyPhantom.ct =			std::cos(emptyPhantom.theta);
  	emptyPhantom.st =			std::sin(emptyPhantom.theta);
  	emptyPhantom.cp =			std::cos(emptyPhantom.phi);
  	emptyPhantom.sp =			std::sin(emptyPhantom.phi);
  	emptyPhantom.em_value = 	1.; // it's deducted afterwards
  	emptyPhantom.em_slope =		0.;
  	emptyPhantom.em_polar =		0.0;
  	emptyPhantom.em_azimut =	0.0;
  	emptyPhantom.em_radial =	0.; //it's deducted afterwards
  	emptyPhantom.zmin =			-30.;
  	emptyPhantom.zmax =			30.;
  	emptyPhantom.name =			"empty Cylinder";

  std::cout<<"About to enter compute norm functions"<<std::endl;
  std::cout<<"Axial smoothing sigma: " << axialSigma << (axialSigma == 0 ? " (disabled)" : "") << std::endl;
  std::cout<<"Transaxial smoothing sigma: " << transaxialSigma << (transaxialSigma == 0 ? " (disabled)" : "") << std::endl;


  computeNormalizationFactors(files,scannerName,outputDir,outputMatrixFileName,
          nRsectorsAngPos,
          nRsectorsAxial,
          invertDetOrder,
          rsectorIdOrder,
          nModulesTransaxial,
          nModulesAxial,
          nSubmodulesTransaxial,
          nSubmodulesAxial,
          nCrystalsTransaxial,
          nCrystalsAxial,
          nLayers,
          nCrystalPerLayer,
          nLayersRptTransaxial,
          nLayersRptAxial,
          myPhantom,
          emptyPhantom,
		transAxialSize,
		axialSize,
		crystalDepth,
	    detectorRadius,
		  outputMatrixFileName+".csv",
		  axialSigma,
		  transaxialSigma);

  return 0;
}


