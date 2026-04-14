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
              << "  -x, --xml <path>         Scanner configuration XML file\n"
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
              << "  " << programName << " -x scanners/CM2L_1ring.xml -i 'data/*.root' -o norm_output\n"
              << "  " << programName << " -x scanners/16x16x2_4rings.xml -i data.root -o out -j 4\n";
}


int main(int argc,char**argv) {

	std::string xmlConfigFile;
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
		else if (arg == "-x" || arg == "--xml") {
			if (i + 1 < argc) {
				xmlConfigFile = argv[++i];
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


  // Load scanner configuration from XML file first (fail-fast on config errors)
  if (xmlConfigFile.empty()) {
      std::cerr << "Error: XML configuration file is required (-x/--xml)\n";
      printUsage(argv[0]);
      return 1;
  }

  ScannerConfig config = ParseScannerXML(xmlConfigFile);
  if (config.name.empty()) {
      std::cerr << "Error: Failed to load scanner configuration from: " << xmlConfigFile << "\n";
      return 1;
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

  // Compute effective detector radius (add half crystal depth offset)
  float effectiveDetectorRadius = config.detectorRadius + config.crystalDepth / config.nLayers * 0.5f + 0.5f;


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

  // Create nCrystalPerLayer array from config
  std::vector<uint32_t> nCrystalPerLayerVec = config.nCrystalPerLayer;

  computeNormalizationFactors(files, config.name, outputDir, outputMatrixFileName,
          config.nRsectorsAngPos,
          config.nRsectorsAxial,
          config.invertDetOrder,
          config.rsectorIdOrder,
          config.nModulesTransaxial,
          config.nModulesAxial,
          config.nSubmodulesTransaxial,
          config.nSubmodulesAxial,
          config.nCrystalsTransaxial,
          config.nCrystalsAxial,
          config.nLayers,
          nCrystalPerLayerVec.data(),
          config.nLayersRptTransaxial,
          config.nLayersRptAxial,
          myPhantom,
          emptyPhantom,
          config.transAxialSize,
          config.axialSize,
          config.crystalDepth,
          effectiveDetectorRadius,
          outputMatrixFileName+".csv",
          axialSigma,
          transaxialSigma);

  return 0;
}


