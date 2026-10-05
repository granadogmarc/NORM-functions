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
              << "  -g, --geom <path>        Scanner file in the CASToR cylindrical PET format (.geom),\n"
              << "                             see scanners/ for examples\n"

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
              << "                             (not recommended when the layers have different\n"
              << "                             transaxial crystal counts, see TECHNICAL_NOTE 9.3)\n"
              << "  -f, --fov-radius <mm>    Transaxial FOV radius: only LORs whose segment crosses\n"
              << "                             the circle of this radius are used and written\n"
              << "                             (default: 300)\n"
              << "  --invert-det-order       Transaxial detector order is reversed in GATE\n"
              << "                             (not part of the CASToR scanner file; default: off)\n"
              << "  --rsector-id-order <0|1> Rsector ID ordering: 0 = transaxial-first (default),\n"
              << "                             1 = axial-first (cubic array)\n"
              << "  -h, --help               Show this help message and exit\n\n"
              << "Examples:\n"
              << "  " << programName << " -g scanners/CM2L_1ring.geom -i 'data/*.root' -o norm_output\n"
              << "  " << programName << " -g scanners/16x16x2_4rings.geom -i data.root -o out -j 4\n";
}


int main(int argc,char**argv) {

	std::string geomFile;
	bool invertDetOrder = false;   // not in the CASToR scanner file
	int rsectorIdOrder = 0;        // not in the CASToR scanner file
	std::string pattern;
	std::string outputMatrixFileName;
	std::string outputDir;
	int numThreads = 0;  // 0 means use default (all available)
	double axialSigma = 0.0;       // Gaussian smoothing sigma for axial normalization (0 = disabled by default)
	double transaxialSigma = 0.0;  // Gaussian smoothing sigma for transaxial normalization (0 = disabled by default)
	double fovRadius = 300.0;      // transaxial FOV radius (mm)

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
		else if (arg == "-g" || arg == "--geom") {
			if (i + 1 < argc) {
				geomFile = argv[++i];
			} else {
				std::cerr << "Error: missing argument after " << arg << "\n";
				return 1;
			}
		}
		else if (arg == "--invert-det-order") {
			invertDetOrder = true;
		}
		else if (arg == "--rsector-id-order") {
			if (i + 1 < argc) {
				rsectorIdOrder = std::atoi(argv[++i]);
				if (rsectorIdOrder != 0 && rsectorIdOrder != 1) {
					std::cerr << "Error: --rsector-id-order must be 0 or 1.\n";
					return 1;
				}
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
		else if (arg == "-f" || arg == "--fov-radius") {
			if (i + 1 < argc) {
				fovRadius = std::atof(argv[++i]);
				if (!(fovRadius > 0)) {
					std::cerr << "Error: FOV radius must be positive.\n";
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


  // Load the scanner file (CASToR .geom) first (fail-fast on config errors)
  if (geomFile.empty()) {
      std::cerr << "Error: scanner file is required (-g/--geom)\n";
      printUsage(argv[0]);
      return 1;
  }

  ScannerGeometry geom;
  {
      std::string geomError;
      if (!ParseCastorGeom(geomFile, geom, geomError)) {
          std::cerr << "Error: cannot use scanner file " << geomFile << ": " << geomError << "\n";
          return 1;
      }
  }
  std::cout << "Loaded scanner file: " << geomFile << std::endl;
  geom.print(std::cout);

  const double innerRadius = *std::min_element(geom.layerRadius.begin(), geom.layerRadius.end());
  if (fovRadius >= innerRadius) {
      std::cerr << "Error: FOV radius (" << fovRadius << " mm) must be smaller than the scanner radius ("
                << innerRadius << " mm)\n";
      return 1;
  }
  std::cout << "Transaxial FOV radius: " << fovRadius << " mm" << std::endl;
  if (!geom.uniformCrystalsTransaxial && transaxialSigma > 0.0) {
      std::cout << "WARNING: transaxial smoothing (-ts) mixes neighbouring radial bins, which are filled by "
                   "different layer pairs when the layers have different transaxial crystal counts "
                   "(see TECHNICAL_NOTE 9.3)" << std::endl;
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

  // Number of crystals per layer (the transaxial crystal count can differ between layers)
  std::vector<uint32_t> nCrystalPerLayerVec(geom.nLayers);
  for (uint32_t l = 0; l < geom.nLayers; ++l) nCrystalPerLayerVec[l] = geom.nCrystalsInLayer(l);

  computeNormalizationFactors(files, geom.name, outputDir, outputMatrixFileName,
          geom.nRsectorsAngPos,
          geom.nRsectorsAxial,
          invertDetOrder,
          rsectorIdOrder,
          geom.nModulesTransaxial,
          geom.nModulesAxial,
          geom.nSubmodulesTransaxial,
          geom.nSubmodulesAxial,
          geom.nCrystalsTransaxial.data(),
          geom.nCrystalsAxial,
          static_cast<uint8_t>(geom.nLayers),
          nCrystalPerLayerVec.data(),
          1,   // nLayersRptTransaxial
          1,   // nLayersRptAxial
          myPhantom,
          emptyPhantom,
          geom,
          outputMatrixFileName+".csv",
          axialSigma,
          transaxialSigma,
          fovRadius);

  return 0;
}


