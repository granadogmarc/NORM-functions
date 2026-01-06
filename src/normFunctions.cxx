#include "normFunctions.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <iostream>
#include <glob.h>


// Helper: compute transaxial strides (keeps consistent packing between loops)
struct TransaxialStrides {
  int strideLayer;
  int strideCrystal;
  int strideSubmodule;
  int strideModule;
  int strideRsector;
};

static inline TransaxialStrides makeTransaxialStrides(uint32_t nModulesTransaxial,
                           uint32_t nSubmodulesTransaxial,
                           uint32_t nCrystalsTransaxial,
                           uint8_t  /*nLayers*/) {
  TransaxialStrides s;
  s.strideLayer = 1;
  // remove layers from transaxial strides (18D packing)
  s.strideCrystal = 1;
  s.strideSubmodule = int(nCrystalsTransaxial);
  s.strideModule = int(nCrystalsTransaxial) * int(nSubmodulesTransaxial);
  s.strideRsector = int(nModulesTransaxial) * int(nSubmodulesTransaxial) * int(nCrystalsTransaxial);
  return s;
}


std::pair<size_t, size_t>
processFile(const std::string &filename,
			matrixRingsComponent		&ringsComponentMatrix,
			double						meanRingsComponentMatrix,
			matrixRingsComponent		&blockTrAComponentMatrix,
			double						meanBlockTrAComponentMatrix,
			vectorRadialComponent		&radialComponentVector,
			double						meanRingComponentVector,
      vectorRingComponent		&ringComponentVector,
			DetectorCounters &detectorEfficyCounts,
      Long64_t &usedLOR,
      Long64_t &totalEvents,
            uint32_t nRsectorsAngPos,
      uint32_t nRsectorsAxial,
            bool invertDetOrder,
      int rsectorIdOrder,
            uint32_t nModulesTransaxial,
      uint32_t nModulesAxial,
            uint32_t nSubmodulesTransaxial,
      uint32_t nSubmodulesAxial,
            uint32_t nCrystalsTransaxial,
      uint32_t nCrystalsAxial,
            uint8_t nLayers,
      uint32_t *nCrystalPerLayer,
            uint32_t nLayersRptTransaxial,
      uint32_t nLayersRptAxial,
            const Phantom &myPhantom,
      const Phantom &emptyPhantom,
  float crystalDepth,
    float detectorRadius,
    float transAxialSize,
    float axialSize)
{
    size_t newKeys = 0;
    size_t localMax = 0;

    TFile *file = TFile::Open(filename.c_str());
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return {0, 0};
    }

    TTree *Coincidences = (TTree *)file->Get("Coincidences");
    if (!Coincidences) {
        std::cerr << "Cannot find TTree 'Coincidences' in file: " << filename << std::endl;
        return {0, 0};
    }

    int rsectorID1, moduleID1, submoduleID1, crystalID1, layerID1;
    int rsectorID2, moduleID2, submoduleID2, crystalID2, layerID2;
    Float_t gPosX1, gPosY1, gPosZ1, gPosX2, gPosY2, gPosZ2;



    Coincidences->SetBranchAddress("rsectorID1", &rsectorID1);
    Coincidences->SetBranchAddress("moduleID1", &moduleID1);
    Coincidences->SetBranchAddress("submoduleID1", &submoduleID1);
    Coincidences->SetBranchAddress("crystalID1", &crystalID1);
    Coincidences->SetBranchAddress("layerID1", &layerID1);
    Coincidences->SetBranchAddress("rsectorID2", &rsectorID2);
    Coincidences->SetBranchAddress("moduleID2", &moduleID2);
    Coincidences->SetBranchAddress("submoduleID2", &submoduleID2);
    Coincidences->SetBranchAddress("crystalID2", &crystalID2);
    Coincidences->SetBranchAddress("layerID2", &layerID2);


    Long64_t nEvents = Coincidences->GetEntries();
    totalEvents += nEvents;

    double PI = 3.14159265358979323846;

    int castPosID =0;

    bool first_entry = true;


	int maxRadialID = int(nRsectorsAngPos*nModulesTransaxial *nSubmodulesTransaxial * nCrystalsTransaxial*nLayers/2);
	int maxTrAID	= nModulesTransaxial *nSubmodulesTransaxial * nCrystalsTransaxial*nLayers*2;


    int maxID =   nRsectorsAngPos * nRsectorsAxial *
                nModulesTransaxial * nModulesAxial *
                nSubmodulesTransaxial * nSubmodulesAxial *
                nCrystalsTransaxial * nCrystalsAxial * nLayers;


    for (Long64_t i = 0; i < nEvents; ++i) {
        Coincidences->GetEntry(i);



        //if (i>pow(10,6)) break;


        TVector3 gPos1 = convertToPosition(323.8, -27.0,-27.0-(nModulesAxial-1)*63*0.5,2*PI/32,layerID1, crystalID1, submoduleID1, moduleID1, rsectorID1,
                    crystalDepth, transAxialSize, axialSize, nLayers, nCrystalsTransaxial, nSubmodulesAxial);
        TVector3 gPos2 = convertToPosition(323.8, -27.0,-27.0-(nModulesAxial-1)*63*0.5,2*PI/32.,layerID2, crystalID2, submoduleID2, moduleID2, rsectorID2,
                    crystalDepth, transAxialSize, axialSize, nLayers, nCrystalsTransaxial, nSubmodulesAxial);

        int castorID1 = ConvertIDcylindrical(nRsectorsAngPos, nRsectorsAxial, invertDetOrder, rsectorIdOrder,
                                             nModulesTransaxial, nModulesAxial, nSubmodulesTransaxial, nSubmodulesAxial,
                                             nCrystalsTransaxial, nCrystalsAxial, nLayers, nCrystalPerLayer,
                                             nLayersRptTransaxial, nLayersRptAxial,
                                             layerID1, crystalID1, submoduleID1, moduleID1, rsectorID1);

        int castorID2 = ConvertIDcylindrical(nRsectorsAngPos, nRsectorsAxial, invertDetOrder, rsectorIdOrder,
                                             nModulesTransaxial, nModulesAxial, nSubmodulesTransaxial, nSubmodulesAxial,
                                             nCrystalsTransaxial, nCrystalsAxial, nLayers, nCrystalPerLayer,
                                             nLayersRptTransaxial, nLayersRptAxial,
                                             layerID2, crystalID2, submoduleID2, moduleID2, rsectorID2);



        if (castorID1 > maxID) std::cout << "Error ID1 is too big ! It is: "<<castorID1<<" with: "<<layerID1<<", "<<crystalID1<<", "<< submoduleID1<<", "<<moduleID1<<", "<< rsectorID1<<std::endl;
        if (castorID2 > maxID) std::cout << "Error ID2 is too big ! It is: "<<castorID2<<" with: "<<layerID2<<", "<<crystalID2<<", "<< submoduleID2<<", "<<moduleID2<<", "<< rsectorID2<<std::endl;


        if (castorID1 < castorID2) {
            std::swap(castorID1, castorID2);
            std::swap(layerID1,layerID2);
            std::swap(crystalID1,crystalID2);
            std::swap(submoduleID1,submoduleID2);
            std::swap(moduleID1,moduleID2);
            std::swap(rsectorID1,rsectorID2);
            std::swap(gPos1, gPos2);
        }

        double R(0);

        Sinogram mySinogram;

        //This is mandatory correction for the detector radius to be adjusted to the real size of the max LOR length possible
        //The layerID goes from inner to outer so the bigger it is the bigger the correction will be


        mySinogram =  ConvertToSinogram(gPos1,gPos2);

    	//FILTER for the FIELD OF VIEW!!!!
    	if (mySinogram.R > 300) continue;
    	R=mySinogram.R;

    	if(std::isnan(R)) {std::cout<<"Error R is nan!"<<std::endl;
    	break;}




        if (filename.find("solidCyl")!=std::string::npos){

        	detectorEfficyCounts.addSubmoduleCount(submoduleID1);
            detectorEfficyCounts.addSubmoduleCount(submoduleID2);
            detectorEfficyCounts.addCrystalCount(crystalID1);
            detectorEfficyCounts.addCrystalCount(crystalID2);
            detectorEfficyCounts.addLayerCount(layerID1);
            detectorEfficyCounts.addLayerCount(layerID2);

            //creating the axial block component:
            if(moduleID1 == moduleID2){
            	if(submoduleID1 == submoduleID2){//TODO this should be generalised to identify the element with axial component

            		ringComponentVector[moduleID1+nModulesAxial*submoduleID1] +=1;
            	}
            }
            //Creating the geometric Axial component once we've filled the values of ring comonents!
            if (i>int(nEvents*0.25)){
            	meanRingComponentVector = meanVector(ringComponentVector);
            	if(first_entry) {
            		auto minIt = std::min_element(ringComponentVector.begin(), ringComponentVector.end());

            		std::cout<<"The number of coincidences up to this point is: "<<nEvents<<" we start computing the ringsComponentMatrix here with the following entries at each bin."<<std::endl;
            		std::cout<<"The bin with less entries has "<<*minIt<<" entries. Poissoinian error = "<<pow(*minIt,-0.5)<<std::endl;
            		std::cout<<"The mean value of the vector is: "<<meanRingComponentVector<<std::endl;
            		first_entry=false;
            	}

            	double blockCorrection = sqrt(meanRingComponentVector*meanRingComponentVector/(ringComponentVector[moduleID1+nModulesAxial*submoduleID1]*ringComponentVector[moduleID2+nModulesAxial*submoduleID2]));
            	ringsComponentMatrix[moduleID1+nModulesAxial*submoduleID1][moduleID2+nModulesAxial*submoduleID2] += mySinogram.Cos_theta*blockCorrection;

            }



        }


        else if (filename.find("annular")!=std::string::npos){



        	float lineIntegral = ComputePhantomLineIntegral(gPos1, gPos2, myPhantom, 1, 1);
        	float emptyIntegral = ComputePhantomLineIntegral(gPos1, gPos2, emptyPhantom, 1, 1);
          float finalIntegral = lineIntegral - emptyIntegral;

          if(!(lineIntegral>0)) std::cout<<"ERROR initial INTEGRAL IS NOT A VALUE above 0!"<<std::endl;
          if(!(emptyIntegral>0)) std::cout<<"ERROR empty INTEGRAL IS NOT A VALUE above 0!"<<std::endl;
          if(!(finalIntegral>0)) {
            std::cout<<"ERROR final INTEGRAL IS NOT A VALUE above 0! Skipping event to avoid division by zero."<<std::endl;
            continue;
          }



        	//Let's computeStrides:

              int ringID1 = moduleID1+nModulesAxial*submoduleID1; //TODO generalise it for all geometric cases
              int ringID2 = moduleID2+nModulesAxial*submoduleID2; //TODO generalise it for all geometric cases

                                          // Guard against zero ring counts or zero matrix entries to avoid NaNs/infs
                                          if (ringComponentVector[ringID1] <= 0.0 || ringComponentVector[ringID2] <= 0.0) {
                                            // missing ring statistics -> skip this event
                                            continue;
                                          }
                                          double denom_block = ringComponentVector[ringID1]*ringComponentVector[ringID2];
                                          double blockCorrection = std::sqrt((long double)meanRingComponentVector*(long double)meanRingComponentVector/ (long double)denom_block);

                                          double geomAxCorrection = 1.0;
                                          if (ringsComponentMatrix[ringID1][ringID2] != 0.0)
                                            geomAxCorrection = meanRingsComponentMatrix / (ringsComponentMatrix[ringID1][ringID2]);
                                          else {
                                            // missing geometric axial data -> skip
                                            continue;
                                          }

              //18D removing layers from transaxial elements (use central helper)
              TransaxialStrides ts = makeTransaxialStrides(nModulesTransaxial, nSubmodulesTransaxial, nCrystalsTransaxial, nLayers);
              int strideLayer     = ts.strideLayer;
              int strideCrystal   = ts.strideCrystal;
              int strideSubmodule = ts.strideSubmodule;
              int strideModule    = ts.strideModule;
              int strideRsector   = ts.strideRsector;

			  //The need to use linear repeaters implies that if there is an axial component there cannot be a transaxial component
			  //this results in no need for divisions over axial repetitions.

			  int ringPosID1 = 0;
			  int ringPosID2 = 0;

			  if (nModulesTransaxial > 1){
				  ringPosID1 += moduleID1*strideModule;
			  	  ringPosID2 += moduleID2*strideModule;
			  }
			  if(nSubmodulesTransaxial > 1){
			  	  ringPosID1 += submoduleID1 * strideSubmodule;
			  	  ringPosID2 += submoduleID2 * strideSubmodule;
			  }
			  if (nCrystalsTransaxial > 1){
				  ringPosID1 += crystalID1 * strideCrystal;
			  	  ringPosID2 += crystalID2 * strideCrystal;
		  	  }

              
/*
			  if (nLayers > 1){
				  ringPosID1 += layerID1 * strideLayer;
			  	  ringPosID2 += layerID2 * strideLayer;
		  }
		  */
			  //Since the radialID is computed as the difference between the two ringPosID,
			  //the ringPosID1 before adding the rsector value already defines the intra
			  //block position for each radialID:
			  int trAID = ringPosID1;

			  ringPosID1 += rsectorID1 * strideRsector;
			  ringPosID2 += rsectorID2 * strideRsector;

			  int totalTransaxial = nRsectorsAngPos*nModulesTransaxial*nSubmodulesTransaxial*nCrystalsTransaxial;//18D removing layers from total transaxial

			  // Compute delta between the two transaxial detector indices
			  int delta = abs(ringPosID1 - ringPosID2);



        // If delta is 0 or equals totalTransaxial then radialID becomes -1 -> debug and skip
        if (delta == 0 || delta == totalTransaxial) {
            int ringPosBefore1 = ringPosID1 - rsectorID1 * strideRsector;
            int ringPosBefore2 = ringPosID2 - rsectorID2 * strideRsector;
            std::cerr << "DEBUG: invalid delta encountered in processFile (delta=" << delta << ") -> radialID would be -1.\n"
                      << " indices: layer1=" << layerID1 << ", module1=" << moduleID1 << ", submodule1=" << submoduleID1 << ", crystal1=" << crystalID1 << ", rsector1=" << rsectorID1 << "\n"
                      << "          layer2=" << layerID2 << ", module2=" << moduleID2 << ", submodule2=" << submoduleID2 << ", crystal2=" << crystalID2 << ", rsector2=" << rsectorID2 << "\n"
                      << " ringPos before rsector: " << ringPosBefore1 << " vs " << ringPosBefore2 << "\n"
                      << " ringPos after rsector:  " << ringPosID1 << " vs " << ringPosID2 << "\n"
                      << " strideRsector=" << strideRsector << ", totalTransaxial=" << totalTransaxial << "\n";
            // Skip this invalid pairing; it should not contribute to transaxial bins
            continue;
        }

        // Fold the radial ID so that opposite orientations map to the same bin
        int radialID = std::min(delta, totalTransaxial - delta)-1;

			  // Safety: ensure radialID stays within bounds
			  // (expected size is totalTransaxial/2)
			  if (radialID < 0 || radialID >= totalTransaxial/2) {
			      std::cerr << "Error: radialID out of bounds: " << radialID << std::endl;

			  }



              // Protect detector efficiency denominators
              double sub1 = detectorEfficyCounts.getSubmoduleCount(submoduleID1);
              double sub2 = detectorEfficyCounts.getSubmoduleCount(submoduleID2);
              double effFactorSubmodule = 1.0;
              if (sub1 > 0.0 && sub2 > 0.0)
                effFactorSubmodule = pow(detectorEfficyCounts.getMaxSubmoduleCount(),2)/(sub1*sub2);

              double cry1 = detectorEfficyCounts.getCrystalCount(crystalID1);
              double cry2 = detectorEfficyCounts.getCrystalCount(crystalID2);
              double effFactorCrystal = 1.0;
              if (cry1 > 0.0 && cry2 > 0.0)
                effFactorCrystal = pow(detectorEfficyCounts.getMaxCrystalCount(),2)/(cry1*cry2);

              double lay1 = detectorEfficyCounts.getLayerCount(layerID1);
              double lay2 = detectorEfficyCounts.getLayerCount(layerID2);
              double effFactorLayer = 1.0;
              if (lay1 > 0.0 && lay2 > 0.0)
                effFactorLayer = pow(detectorEfficyCounts.getMaxLayerCount(),2)/(lay1*lay2);

			  double effNormFactor = effFactorSubmodule*effFactorCrystal*effFactorLayer;


              radialComponentVector[radialID] += blockCorrection*geomAxCorrection*effNormFactor/finalIntegral;

              blockTrAComponentMatrix[radialID][trAID] += blockCorrection*geomAxCorrection*effNormFactor/finalIntegral;


        }
      else{
        	std::cerr << "Error: The file doesn't contain information about the source." << "\n";
        	std::cerr << "Include 'annular' or 'solidCyl' in your filename to compute the norm factors accordingly" << "\n";
        	exit(EXIT_FAILURE);
        }

    }


    file->Close();
    return {newKeys, localMax};
}

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
) {
    // Compute max possible CastorID
    //uint64_t totalCrystals = 0;
    //for (int l = 0; l < nLayers; ++l) totalCrystals += nCrystalPerLayer[l];
    uint64_t maxCastorID = nRsectorsAngPos * nRsectorsAxial *
                           nModulesTransaxial * nModulesAxial *
                           nSubmodulesTransaxial * nSubmodulesAxial *
                           nCrystalsTransaxial * nCrystalsAxial *
                           nLayersRptTransaxial * nLayersRptAxial;

std::cout<<"Max castorID is "<<maxCastorID<<std::endl;
    std::vector<DetectorIndices> lut(maxCastorID);

    // Fill LUT
    uint32_t castorID;


    for (uint32_t layerID = 0; layerID < nLayers; ++layerID) {
      for (uint32_t crystalID = 0; crystalID < nCrystalsTransaxial*nCrystalsAxial; ++crystalID) {
        for (uint32_t submoduleID = 0; submoduleID < nSubmodulesTransaxial*nSubmodulesAxial; ++submoduleID) {
          for (uint32_t moduleID = 0; moduleID < nModulesTransaxial*nModulesAxial; ++moduleID) {
            for (uint32_t rsectorID = 0; rsectorID < nRsectorsAngPos*nRsectorsAxial; ++rsectorID) {

                castorID = ConvertIDcylindrical(
                    nRsectorsAngPos, nRsectorsAxial,
                    invertDetOrder, rsectorIdOrder,
                    nModulesTransaxial, nModulesAxial,
                    nSubmodulesTransaxial, nSubmodulesAxial,
                    nCrystalsTransaxial, nCrystalsAxial,
                    nLayers, nCrystalPerLayer,
                    nLayersRptTransaxial, nLayersRptAxial,
                    layerID, crystalID, submoduleID, moduleID, rsectorID);

                lut[castorID] ={
                		(int)rsectorID,
                		(int)moduleID,
                		(int)submoduleID,
						(int)crystalID,
						(int) layerID};
           	}
          }
       	}
   	  }
    }
    std::cout<<"LUT for castorIDs has been written with max castorID = "<<castorID<<std::endl;
    return lut;
}




void computeNormalizationFactors( const std::vector<std::string> &filenames, const std::string& scannerName,
    const std::string& outputDir,
    const std::string& outputMatrixFileName,
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
    float crystalDepth,
    float detectorRadius,
    const std::string &outCSV
      )
{

	double PI = 3.14159265358979323846;

	int maxCastorID = nRsectorsAngPos * nRsectorsAxial *
	                           nModulesTransaxial * nModulesAxial *
	                           nSubmodulesTransaxial * nSubmodulesAxial *
	                           nCrystalsTransaxial * nCrystalsAxial *
							   nLayers;

	int transaxialElements = 0;
	if (nModulesTransaxial>1) transaxialElements +=1;
	if (nSubmodulesTransaxial>1) transaxialElements +=1;
	if (nCrystalsTransaxial>1) transaxialElements +=1;
	if (nLayers>1) transaxialElements +=1;

	int maxSymID = maxCastorID*maxCastorID/(nRsectorsAngPos * nRsectorsAxial);
	int maxRingID	= nModulesAxial*nSubmodulesAxial*nCrystalsAxial;

	int maxRadialID = int(nRsectorsAngPos*nModulesTransaxial *nSubmodulesTransaxial * nCrystalsTransaxial*nLayers/2);
	int maxTrAID	= nModulesTransaxial *nSubmodulesTransaxial * nCrystalsTransaxial*nLayers;

	std::cout<<"maxRadial ID = "<<maxRadialID<<" maxTrAID = "<<maxTrAID<<std::endl;

	matrixRingsComponent		ringsComponentMatrix(maxRingID, std::vector<double>(maxRingID));
  matrixRingsComponent		blockTrAComponentMatrix(maxRadialID,std::vector<double>(maxTrAID));
  vectorRadialComponent		radialComponentVector(maxRadialID,0.0);
	vectorRingComponent			ringComponentVector(maxRingID,0.0);

	DetectorCounters detectorEfficyCounts(nSubmodulesTransaxial*nSubmodulesAxial, nCrystalsTransaxial*nCrystalsAxial,nLayers);





	size_t mapSize       = 0;    // track how many keys we’ve ever inserted

	size_t overallMaxHits = 0;   // track the high‐water mark of numHits

	Long64_t usedLOR(0), totalEvents(0);






  double meanRingComponentVector = 1;
  double meanRingsComponentMatrix = 1;
  double meanRadialComponentVector = 1;
  double meanBlockTrAComponentMatrix = 1;

  

  for(auto &fn : filenames) {


	  std::cout<<"Before opening the file usedLOR = "<<usedLOR<<" and totalEvents = "<<totalEvents<<std::endl;


	size_t newKeys, localMaxHits;
	std::tie(newKeys, localMaxHits)
  = processFile(fn,
        ringsComponentMatrix,
        meanRingsComponentMatrix,
        blockTrAComponentMatrix,
        meanBlockTrAComponentMatrix,
        radialComponentVector,
        meanRadialComponentVector,
        ringComponentVector,
        detectorEfficyCounts,
        usedLOR,
        totalEvents,
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
			crystalDepth,
            detectorRadius,
            transAxialSize,
            axialSize);

	   mapSize       += newKeys;
	   overallMaxHits = std::max(overallMaxHits, localMaxHits);

	   std::cout
	     << "Processed “" << fn << "”"
	     << "  →  total keys="  << mapSize
	     << ", maxHits="       << overallMaxHits
	     << "\n";

	   //We compute these values here so we can use them for the corrections for the annular source
	   if (fn.find("solidCyl")!=std::string::npos){
	   meanRingComponentVector = meanVector(ringComponentVector);
	   meanRingsComponentMatrix = meanMatrix(ringsComponentMatrix);
	   }

	   if (fn.find("annular")!=std::string::npos){
		   meanRadialComponentVector = meanVector(radialComponentVector);
		   meanBlockTrAComponentMatrix = meanMatrix(blockTrAComponentMatrix);//
	   }

  }




  // Commented out creation of most output Cdf files to avoid writing them.
  // Keep files that contain "CB" and "effBfGAfGTrAf" in their names.
  // Aggregated CSVs will be written after building the normalization arrays

  // Disabled: blockFact, geomAxFact, geomTrAFact, eff, intTrAf
  // std::ofstream normFile_blockFact(outputDir+outputMatrixFileName+"_blockFact_df.Cdf", std::ios::binary);
  // std::ofstream normFile_geomAxFact(outputDir+outputMatrixFileName+"_geomAxFact_df.Cdf", std::ios::binary);
  // std::ofstream normFile_geomTrAFact(outputDir+outputMatrixFileName+"_geomTrAFact_df.Cdf", std::ios::binary);
  // std::ofstream normFile_eff(outputDir+outputMatrixFileName+"_eff_df.Cdf", std::ios::binary);
  // std::ofstream normFile_intTrAf(outputDir+outputMatrixFileName+"_intTrAf_df.Cdf", std::ios::binary);

  // Keep CB outputs
  std::ofstream normFile_CB(outputDir+outputMatrixFileName+"_CB_df.Cdf", std::ios::binary);

  // Disabled: effBf and related intermediate products, except keep effBfGAfGTrAf
  // std::ofstream normFile_effBf(outputDir+outputMatrixFileName+"_effBf_df.Cdf", std::ios::binary);
  // std::ofstream normFile_effBfGAf(outputDir+outputMatrixFileName+"_effBfGAf_df.Cdf", std::ios::binary);
  std::ofstream normFile_effBfGAfGTrAf(outputDir+outputMatrixFileName+"_effBfGAfGTrAf_df.Cdf", std::ios::binary);
  std::ofstream normFile_CBsqrBfnI(outputDir+outputMatrixFileName+"_CBsqrBfnI_df.Cdf", std::ios::binary);




  std::cout<<"Here we still know what the outputMatrixfilename is "<<outputMatrixFileName<<std::endl;

  double meanCrystalLORcount = 0;
  int nCrystalLORcount = 0;





  // ------------------------------------------------------------------
  // Compute per-row/ per-radial means for component-wise normalization
  // blockCorrection: per ring1 mean over ring2
  // geomAxCorrection: per ring1 mean over ring2
  // transaxialGeomNormFactor: global mean over radialIDs
  // interferenceTraFactor: per radialID mean over trAID (after transaxial normalization)
  // ------------------------------------------------------------------

  // Compute global means for normalization (global per-component normalization)
  // 1) blockCorrection and geomAxCorrection: global means over all valid ring pairs
  long double sb_total = 0.0L; uint64_t cb_total = 0;
  long double sg_total = 0.0L; uint64_t cg_total = 0;
  for (int r1 = 0; r1 < maxRingID; ++r1) {
    double rc1 = ringComponentVector[r1];
    for (int r2 = 0; r2 < maxRingID; ++r2) {
      double rc2 = ringComponentVector[r2];
      if (rc1 > 0.0 && rc2 > 0.0) {
        long double bc = std::sqrt((long double)meanRingComponentVector*(long double)meanRingComponentVector/((long double)rc1*(long double)rc2));
        sb_total += bc; ++cb_total;
      }
      double rm = ringsComponentMatrix[r1][r2];
      if (rm != 0.0) {
        long double gc = (long double)meanRingsComponentMatrix / (long double)rm;
        sg_total += gc; ++cg_total;
      }
    }
  }
  double mean_block_global = (cb_total ? double(sb_total / cb_total) : 1.0);
  double mean_geomAx_global = (cg_total ? double(sg_total / cg_total) : 1.0);

  // 2) transaxialGeomNormFactor: compute normalized transaxial factor per radialID (global-normalized vector)
  size_t nRad = radialComponentVector.size();
  std::vector<double> norm_transaxial(nRad, 1.0);
  long double s_t = 0.0L; uint64_t c_t = 0;
  for (size_t r = 0; r < nRad; ++r) {
    double rv = radialComponentVector[r];
    if (rv != 0.0) {
      long double t = (long double)meanRadialComponentVector / (long double)rv;
      s_t += t; ++c_t;
    }
  }
  long double mean_t = (c_t ? s_t / c_t : 1.0L);
  for (size_t r = 0; r < nRad; ++r) {
    double rv = radialComponentVector[r];
    if (rv != 0.0)
      norm_transaxial[r] = double(((long double)meanRadialComponentVector / (long double)rv) / mean_t);
    else
      norm_transaxial[r] = 1.0;
  }

  // 3) interference: compute a single global mean over all (radial,trA) valid entries
  long double si_total = 0.0L; uint64_t ci_total = 0;
  for (size_t r = 0; r < nRad; ++r) {
    double ntr = norm_transaxial[r];
    for (size_t t = 0; t < blockTrAComponentMatrix[r].size(); ++t) {
      double b = blockTrAComponentMatrix[r][t];
      if (b != 0.0 && ntr != 0.0) {
        long double interf = (long double)meanBlockTrAComponentMatrix / ((long double)ntr * (long double)b);
        si_total += interf; ++ci_total;
      }
    }
  }
  double mean_interf_global = (ci_total ? double(si_total / ci_total) : 1.0);

  DetectorIndices d1,d2;
  int rsectorDiff = 0;
  double Navg = 0;
  int64_t LORinFOV=0;



  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////



  for (uint32_t layerID1 = 0; layerID1 < nLayers; ++layerID1) {
	  for (uint32_t layerID2 = 0; layerID2 < nLayers; ++layerID2) {
		  for (uint32_t crystalID1 = 0; crystalID1 < nCrystalsTransaxial*nCrystalsAxial; ++crystalID1) {
			  for (uint32_t crystalID2 = 0; crystalID2 < nCrystalsTransaxial*nCrystalsAxial; ++crystalID2){
				  for (uint32_t submoduleID1 = 0; submoduleID1 < nSubmodulesTransaxial*nSubmodulesAxial; ++submoduleID1) {
					  for (uint32_t submoduleID2 = 0; submoduleID2 < nSubmodulesTransaxial*nSubmodulesAxial; ++submoduleID2) {
						  for (uint32_t moduleID1 = 0; moduleID1 < nModulesTransaxial*nModulesAxial; ++moduleID1) {
							  for (uint32_t moduleID2 = 0; moduleID2 < nModulesTransaxial*nModulesAxial; ++moduleID2) {
								  for (uint32_t rsectorID1 = 0; rsectorID1 < nRsectorsAngPos*nRsectorsAxial; ++rsectorID1) {
									  for (uint32_t rsectorID2 = 0; rsectorID2 < nRsectorsAngPos*nRsectorsAxial; ++rsectorID2) {

                        //Applying an "R = 300 mm" filter

                        // Compute transaxial component of rsectorID (handles axial repeats)
                        int rsectorTrs1, rsectorTrs2;
                        if (rsectorIdOrder == 0) {
                          rsectorTrs1 = rsectorID1 % nRsectorsAngPos;
                          rsectorTrs2 = rsectorID2 % nRsectorsAngPos;
                        } else {
                          // axial-first ordering
                          rsectorTrs1 = rsectorID1 / nRsectorsAxial;
                          rsectorTrs2 = rsectorID2 / nRsectorsAxial;
                        }

                        // Circular difference on transaxial sectors (wrap-around aware)
                        int absDiff = std::abs(rsectorTrs1 - rsectorTrs2);
                        int circDiff = std::min(absDiff, int(nRsectorsAngPos) - absDiff);

                        if (circDiff < 4)
                          continue;
                        if (circDiff == 4) {
                          if ((nSubmodulesTransaxial != 0) && (abs(submoduleID1 - submoduleID2) < 13))
                            continue;
                          else if ((nCrystalsTransaxial != 0) && (abs(crystalID1 - crystalID2) < 13))
                            continue;
                        }

                        //Computing the castor ID's

                        int castorID1 = ConvertIDcylindrical(nRsectorsAngPos, nRsectorsAxial, invertDetOrder, rsectorIdOrder,
                          nModulesTransaxial, nModulesAxial, nSubmodulesTransaxial, nSubmodulesAxial,
                          nCrystalsTransaxial, nCrystalsAxial, nLayers, nCrystalPerLayer,
                          nLayersRptTransaxial, nLayersRptAxial,
                          layerID1, crystalID1, submoduleID1, moduleID1, rsectorID1);

                        int castorID2 = ConvertIDcylindrical(nRsectorsAngPos, nRsectorsAxial, invertDetOrder, rsectorIdOrder,
                          nModulesTransaxial, nModulesAxial, nSubmodulesTransaxial, nSubmodulesAxial,
                          nCrystalsTransaxial, nCrystalsAxial, nLayers, nCrystalPerLayer,
                          nLayersRptTransaxial, nLayersRptAxial,
                          layerID2, crystalID2, submoduleID2, moduleID2, rsectorID2);

                        // Only process unordered pairs once (castorID1 > castorID2)
                        if (castorID1 <= castorID2)
                          continue;

                        // Count this unique LOR only when we actually write it
                        // (increment moved down to where entries are written)




                                          TransaxialStrides ts = makeTransaxialStrides(nModulesTransaxial, nSubmodulesTransaxial, nCrystalsTransaxial, nLayers);
                                          int strideLayer     = ts.strideLayer;
                                          int strideCrystal   = ts.strideCrystal;
                                          int strideSubmodule = ts.strideSubmodule;
                                          int strideModule    = ts.strideModule;
                                          int strideRsector   = ts.strideRsector;

										  int ringPosID1 = 0;
										  int ringPosID2 = 0;

										  if (nModulesTransaxial > 1){
											  ringPosID1 += moduleID1*strideModule;
										  	  ringPosID2 += moduleID2*strideModule;
										  }
										  if(nSubmodulesTransaxial > 1){
										  	  ringPosID1 += submoduleID1 * strideSubmodule;
										  	  ringPosID2 += submoduleID2 * strideSubmodule;
										  }
										  if (nCrystalsTransaxial > 1){
											  ringPosID1 += crystalID1 * strideCrystal;
										  	  ringPosID2 += crystalID2 * strideCrystal;
									  	  }

										  /*
                       if (nLayers > 1){
											  ringPosID1 += layerID1 * strideLayer;
										  	  ringPosID2 += layerID2 * strideLayer;
									  }*/
										  //Since the radialID is computed as the difference between the two ringPosID,
										  //the ringPosID1 before adding the rsector value already defines the intra
										  //block position for each radialID:
										  int trAID = ringPosID1;
										  ringPosID1 += rsectorID1 * strideRsector;
										  ringPosID2 += rsectorID2 * strideRsector;

										  int totalTransaxial = nRsectorsAngPos*nModulesTransaxial*nSubmodulesTransaxial*nCrystalsTransaxial;// *nLayers; 18D removign nlayersfrom transaxial

										  // Compute delta between the two transaxial detector indices
										  int delta = abs(ringPosID1 - ringPosID2);

                      // If delta is 0 or equals totalTransaxial then radialID becomes -1 -> debug and skip
                      if (delta == 0 || delta == totalTransaxial) {
                          int ringPosBefore1 = ringPosID1 - rsectorID1 * strideRsector;
                          int ringPosBefore2 = ringPosID2 - rsectorID2 * strideRsector;
                          std::cerr << "DEBUG: invalid delta encountered (delta=" << delta << ") -> radialID would be -1.\n"
                                    << " indices: layer1=" << layerID1 << ", module1=" << moduleID1 << ", submodule1=" << submoduleID1 << ", crystal1=" << crystalID1 << ", rsector1=" << rsectorID1 << "\n"
                                    << "          layer2=" << layerID2 << ", module2=" << moduleID2 << ", submodule2=" << submoduleID2 << ", crystal2=" << crystalID2 << ", rsector2=" << rsectorID2 << "\n"
                                    << " ringPos before rsector: " << ringPosBefore1 << " vs " << ringPosBefore2 << "\n"
                                    << " ringPos after rsector:  " << ringPosID1 << " vs " << ringPosID2 << "\n"
                                    << " strideRsector=" << strideRsector << ", totalTransaxial=" << totalTransaxial << "\n";
                          // Skip this invalid pairing; it should not contribute to transaxial bins
                          continue;
                      }

                      // Fold the radial ID so that opposite orientations map to the same bin
                      int radialID = std::min(delta, totalTransaxial - delta)-1;
										  // Safety: ensure radialID stays within bounds
										  // (expected size is totalTransaxial/2)
										  if (radialID < 0 || radialID > totalTransaxial*0.5) {
										      std::cerr << "Error: radialID out of bounds: " << radialID << std::endl;
										  }
										  else if(radialID == totalTransaxial*0.5-1){
											  radialID = totalTransaxial*0.5-2; //Attempt to try to repair the wierd last bin effect.
										  }
										  if (trAID < 0 || trAID >= maxTrAID) {
										      std::cerr << "BAD trAID = " << trAID
										                << " (max=" << maxTrAID << ")\n";
										      std::abort();
										  }


										  int ringID1 = moduleID1+nModulesAxial*submoduleID1; //TODO generalise it for all geometric cases
										  int ringID2 = moduleID2+nModulesAxial*submoduleID2; //TODO generalise it for all geometric cases

										  double blockCorrection = sqrt(meanRingComponentVector*meanRingComponentVector/(ringComponentVector[ringID1]*ringComponentVector[ringID2]));
										  double geomAxCorrection = meanRingsComponentMatrix/(ringsComponentMatrix[ringID1][ringID2]);


										  double effFactorSubmodule = pow(detectorEfficyCounts.getMaxSubmoduleCount(),2)/(detectorEfficyCounts.getSubmoduleCount(submoduleID1)*detectorEfficyCounts.getSubmoduleCount(submoduleID2));
										  double effFactorCrystal = pow(detectorEfficyCounts.getMaxCrystalCount(),2)/(detectorEfficyCounts.getCrystalCount(crystalID1)*detectorEfficyCounts.getCrystalCount(crystalID2));
										  double effFactorLayer = pow(detectorEfficyCounts.getMaxLayerCount(),2)/(detectorEfficyCounts.getLayerCount(layerID1)*detectorEfficyCounts.getLayerCount(layerID2));

										  double effNormFactor = effFactorSubmodule*effFactorCrystal*effFactorLayer;



                      double transaxialGeomNormFactor = 1.0;
                      if (radialID >= 0 && size_t(radialID) < norm_transaxial.size())
                        transaxialGeomNormFactor = norm_transaxial[radialID];

                      double norm_block = (mean_block_global != 0.0) ? blockCorrection / mean_block_global : blockCorrection;
                      double norm_geomAx = (mean_geomAx_global != 0.0) ? geomAxCorrection / mean_geomAx_global : geomAxCorrection;

                      double interferenceTraFactor = 0.0;
                      if (transaxialGeomNormFactor != 0.0 && blockTrAComponentMatrix[radialID][trAID] != 0.0) {
                        interferenceTraFactor = meanBlockTrAComponentMatrix / (transaxialGeomNormFactor * blockTrAComponentMatrix[radialID][trAID]);
                        if (mean_interf_global != 0.0)
                          interferenceTraFactor /= mean_interf_global;
                      }

                      double CBasedNF = effNormFactor * norm_block * norm_geomAx * transaxialGeomNormFactor * interferenceTraFactor;

                      // Increment written-LOR counter and write normalized entries
                      LORinFOV++;
                      writeNormEntryNormMatrixFile(normFile_CBsqrBfnI, castorID2, castorID1, effNormFactor*norm_block*norm_block*norm_geomAx*transaxialGeomNormFactor); // disabled
                      writeNormEntryNormMatrixFile(normFile_effBfGAfGTrAf, castorID2, castorID1, effNormFactor*norm_block*norm_geomAx*transaxialGeomNormFactor);
                      writeNormEntryNormMatrixFile(normFile_CB, castorID2, castorID1, CBasedNF);



									  }
								  }
							  }
						  }
					  }
				  }
			  }
		  }
	  }
  }


  //std::cout<<"Average number of coincidences per LOR is "<<Navg<<" with usedLOR "<<usedLOR<<" while LOR in the FOV are "<< LORinFOV <<std::endl;
  std::cout<<"Total number of events so far is: "<<totalEvents<<std::endl;

  // Diagnostic: report radial-vector occupancy and key geometry/stride params
  {
    size_t firstNonZero = radialComponentVector.size();
    size_t lastNonZero = 0;
    size_t nonZeroCount = 0;
    for (size_t i = 0; i < radialComponentVector.size(); ++i) {
      if (radialComponentVector[i] != 0.0) {
        ++nonZeroCount;
        if (firstNonZero == radialComponentVector.size()) firstNonZero = i;
        lastNonZero = i;
      }
    }

    int64_t totalTransaxial = int64_t(nRsectorsAngPos) * nModulesTransaxial * nSubmodulesTransaxial * nCrystalsTransaxial;
    int64_t strideRsector_with_layers = int64_t(nModulesTransaxial) * nSubmodulesTransaxial * nCrystalsTransaxial * nLayers;
    int64_t strideRsector_no_layers = int64_t(nModulesTransaxial) * nSubmodulesTransaxial * nCrystalsTransaxial;

    std::cout << "Radial diagnostic: nonZeroCount=" << nonZeroCount
              << ", firstNonZero=" << (firstNonZero == radialComponentVector.size() ? -1 : (int)firstNonZero)
              << ", lastNonZero=" << (nonZeroCount ? (int)lastNonZero : -1)
              << ", maxRadialID=" << maxRadialID
              << ", totalTransaxial=" << totalTransaxial
              << ", nRsectorsAngPos=" << nRsectorsAngPos
              << ", nModulesTransaxial=" << nModulesTransaxial
              << ", nSubmodulesTransaxial=" << nSubmodulesTransaxial
              << ", nCrystalsTransaxial=" << nCrystalsTransaxial
              << ", nLayers=" << (int)nLayers
              << std::endl;

    std::cout << "StrideRsector (with layers)=" << strideRsector_with_layers
              << ", (no layers)=" << strideRsector_no_layers << std::endl;

    if (nonZeroCount) {
      std::cout << "Radial sample around first non-zero index (up to 8 entries):\n";
      size_t start = firstNonZero;
      size_t end = std::min(radialComponentVector.size(), start + 8);
      for (size_t i = start; i < end; ++i)
        std::cout << "  idx=" << i << ", val=" << radialComponentVector[i] << std::endl;
    }
  }

  // --- aggregated CSV outputs (indexed access) ---------------------------------
  // These are much smaller than per-LOR logging and can be indexed by user.
  // 1) block_geom: per ring pair -> blockCorrection, geomAxCorrection
  {
    std::ofstream csvBlockGeom(outputDir + outputMatrixFileName + "_block_geom.csv");
    csvBlockGeom << "ring1,ring2,blockCorrection,geomAxCorrection\n";
    for (size_t i = 0; i < ringsComponentMatrix.size(); ++i) {
      for (size_t j = 0; j < ringsComponentMatrix[i].size(); ++j) {
        double blockCorr = 0.0;
        if (ringComponentVector[i] != 0.0 && ringComponentVector[j] != 0.0)
          blockCorr = sqrt(meanRingComponentVector * meanRingComponentVector / (ringComponentVector[i] * ringComponentVector[j]));
        double geomAxCorr = 0.0;
        if (ringsComponentMatrix[i][j] != 0.0)
          geomAxCorr = meanRingsComponentMatrix / (ringsComponentMatrix[i][j]);
        csvBlockGeom << i << "," << j << "," << blockCorr << "," << geomAxCorr << "\n";
      }
    }
    csvBlockGeom.close();
  }

  // 2) transaxial: per radialID -> transaxialGeomNormFactor
  {
    std::ofstream csvTrans(outputDir + outputMatrixFileName + "_transaxial.csv");
    csvTrans << "radialID,transaxialGeomNormFactor\n";
    for (size_t r = 0; r < radialComponentVector.size(); ++r) {
      double tgnf = 0.0;
      if (radialComponentVector[r] != 0.0)
        tgnf = meanRadialComponentVector / radialComponentVector[r];
      csvTrans << r << "," << tgnf << "\n";
    }
    csvTrans.close();
  }

  // 3) interference: per (radialID, trAID) -> interferenceTraFactor
  {
    std::ofstream csvInterf(outputDir + outputMatrixFileName + "_interference.csv");
    csvInterf << "radialID,trAID,interferenceTraFactor\n";
    for (size_t r = 0; r < blockTrAComponentMatrix.size(); ++r) {
      for (size_t t = 0; t < blockTrAComponentMatrix[r].size(); ++t) {
        double tgnf = 0.0;
        if (radialComponentVector[r] != 0.0)
          tgnf = meanRadialComponentVector / radialComponentVector[r];
        double denom = tgnf * blockTrAComponentMatrix[r][t];
        double interf = 0.0;
        if (denom != 0.0)
          interf = meanBlockTrAComponentMatrix / denom;
        csvInterf << r << "," << t << "," << interf << "\n";
      }
    }
    csvInterf.close();
  }
  // ---------------------------------------------------------------------------


  std::cout<<"Do we know know what the outputMatrixfilename is "<<outputMatrixFileName<<std::endl;



    // (headers and binary closes handled below)

    std::cout<<"Creating header for "<<outputMatrixFileName<<"_effBfGAfGTrAf"<<std::endl;
    writeCdhHeaderNormMatrixFile(outputDir,outputMatrixFileName+"_effBfGAfGTrAf",
                        scannerName,
                        LORinFOV);

  std::cout<<"Creating header for "<<outputMatrixFileName<<"_CBsqrBfnI"<<std::endl;
   writeCdhHeaderNormMatrixFile(outputDir,outputMatrixFileName+"_CBsqrBfnI",
                         scannerName,
                         LORinFOV);

  std::cout<<"Creating header for "<<outputMatrixFileName<<"_CB"<<std::endl;
  writeCdhHeaderNormMatrixFile(outputDir,outputMatrixFileName+"_CB",
                      scannerName,
                      LORinFOV);



             normFile_CB.close();
             normFile_effBfGAfGTrAf.close();
             normFile_CBsqrBfnI.close();          
}


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
                              int32_t	rsectorID)
{
  // Castor ID definition
  uint32_t castorID = 0;
  uint8_t  layer = layerID;

    // layerID represents the actual layer level

    // add the number of crystals contained in previous layers as
    // CASToR indexes all crystals of a layer ring before the next layer


  	if (layer > 0) castorID += nCrystalPerLayer[layer-1]*layer;


    int32_t nTrsCrystalsPerSubmodule = nCrystalsTransaxial;
    int32_t nTrsCrystalsPerModule = nTrsCrystalsPerSubmodule * nSubmodulesTransaxial;
    int32_t nTrsCrystalsPerRsector = nTrsCrystalsPerModule * nModulesTransaxial;
    int32_t nCrystalsPerRing = nTrsCrystalsPerRsector * nRsectorsAngPos;


    // Rsector axial(=ring) and transaxial(=angular) ID
    // Fastest ordering orientation (axial or transaxial) depends on the repeaters
    int32_t rsectorAxlID = 0 ;
    int32_t rsectorTrsID = 0 ;

    // standard, transaxial first
    if(a_rsectorIdOrder == 0)
    {
      rsectorAxlID = rsectorID/nRsectorsAngPos ;
      rsectorTrsID = (int32_t)(rsectorID%nRsectorsAngPos) ;
    }
    else // using cubic array, axial first
    {
      rsectorAxlID = rsectorID%nRsectorsAxial ;
      rsectorTrsID = (int32_t)(rsectorID/nRsectorsAxial) ;
    }


    // Compute axial ID
    int32_t ringID = rsectorAxlID * nModulesAxial * nSubmodulesAxial * nCrystalsAxial
                   + (int32_t)(moduleID/nModulesTransaxial) * nSubmodulesAxial * nCrystalsAxial
                   + (int32_t)(submoduleID/nSubmodulesTransaxial) * nCrystalsAxial
                   + (int32_t)(crystalID/nCrystalsTransaxial);

    // Recover transaxial ID for each element
    moduleID = moduleID % nModulesTransaxial;
    submoduleID = submoduleID % nSubmodulesTransaxial;
    crystalID = crystalID % nCrystalsTransaxial;

    // Reverse transaxial ordering
    if( a_invertDetOrder )
    {
      moduleID    = nModulesTransaxial-1 - moduleID;
      submoduleID = nSubmodulesTransaxial-1 - submoduleID;
      crystalID   = nCrystalsTransaxial-1 - crystalID;
    }

    // Compute final ID
    castorID += nCrystalsPerRing * ringID
             +  nTrsCrystalsPerRsector * rsectorTrsID
             +  nTrsCrystalsPerModule * moduleID
             +  nTrsCrystalsPerSubmodule * submoduleID
             +  crystalID;




  return castorID;
}





//---------------------------------------------------------------------
// Helper Function 2: ComputeZBounds
//---------------------------------------------------------------------
/// Given two 3D positions, returns a pair (zmin, zmax)
std::pair<float, float> ComputeZBounds(TVector3& pos1,TVector3& pos2)
{
    float zmin = std::min(pos1.Z(), pos2.Z());
    float zmax = std::max(pos1.Z(), pos2.Z());
    return {zmin, zmax};
}

//---------------------------------------------------------------------
// Helper Function 3: ComputeDistance
//---------------------------------------------------------------------
/// Computes the Euclidean distance between two 3D points.
float ComputeDistance(TVector3&  a,
		TVector3&  b)
{/*
    float sumSq = 0.0f;
    for (size_t i = 0; i < 3; ++i)
        sumSq += (b[i] - a[i]) * (b[i] - a[i]);
        */
	TVector3 AB = b-a;
    return AB.Mag();
}


//---------------------------------------------------------------------
// Helper Function 4: ComputePhantomLineIntegral
//---------------------------------------------------------------------
/// Computes the line integral contribution through a phantom object,
/// given the start (a) and end (b) positions, the phantom parameters,
/// the material type, and ToF (if applicable).

float ComputePhantomLineIntegral (TVector3&  a,
								  TVector3& b,
                                 const Phantom& phantom,
                                 int mat_data,
                                 Double_t ToF)
{
    float dels = ComputeDistance(a, b);

    double c0 = 0.0, c1 = 0.0, c2 = 0.0;

    // 2D ellipse in XY-plane
    c0 = (a.X() * a.X()) + (a.Y() * a.Y()) - (phantom.half_axis[0] * phantom.half_axis[0]);
    c1 = (b.X() - a.X()) * a.X() + (b.Y() - a.Y()) * a.Y();
    c2 = (b.X() - a.X()) * (b.X() - a.X()) + (b.Y() - a.Y()) * (b.Y() - a.Y());

    double determ = c1 * c1 - c0 * c2;

    if (determ <= 0)
        return 0.0f;

    double path1 = (-c1 + std::sqrt(determ)) / c2;
    double path2 = (-c1 - std::sqrt(determ)) / c2;

    if (std::isnan(path1) || std::isnan(path2)) {
        std::cerr << "NaN detected! c0=" << c0 << ", c1=" << c1 << ", c2=" << c2 << ", determ=" << determ << "\n";
    }

    if (path1 < path2)
        path1 = path2;
    if (path2 < 0)
        path2 = 0.0;

    float value = 0.0f;
    if (mat_data == EmissionSimul)
        value = phantom.em_value;
    else if (mat_data == AttenuationSimul || mat_data == TransmissionSimul)
        value = phantom.mu_value;

    if (value == 0.0f)
        return 0.0f;

    return static_cast<float>(dels * std::fabs(path1 - path2) * value);
}

float emRad(float a, float b, float c, float s, int rad_type)
{
	int em_radial_max = 4;
  float delta,f,u,v,R;
  R = a+b*s+c*pow(s,2.0);
  delta = 4.0*a*c - pow(b,2.0);
  f = 0.0;
  if (rad_type < 0 || rad_type > em_radial_max) return f;
  else if (rad_type == 1 || rad_type == 3) {
    u = 2*c*s + b;
    if (delta > 0.0) u /= sqrt(delta);
    else if (delta < 0.0) u /= sqrt(-delta);
    f = 0;
    v = 0;
    if (rad_type == 1) {
      f = (2.0*c*s+b)/(4.0*c);
      v = delta/(8.0*c);
    } else if (rad_type == 3) {
      f = (R/(8.0*c)+3.0*delta/pow(8.0*c,2.0))*(2.0*c*s+b);
      v = 3.0*pow(delta,2.0)/(128.0*pow(c,2.0));
    }
    f = s - f*sqrt(R);
    if (c > 0.0 && delta > 0.0) f -= v/sqrt(c)*log(u+sqrt(1.0+pow(u,2.0)));
    else if (c < 0.0 && delta < 0.0) f += v/sqrt(-c)*asin(u);
    else if (c > 0.0 && delta == 0.0) f -= v/sqrt(c)*log(u);
  } else if (rad_type == 2) {
    f = (1.0 - a)*s - 0.5*b*pow(s,2.0) - 1.0/3.0*c*pow(s,3.0);
  } else if (rad_type == 4) {
    f = (1.0 - pow(a,2.0))*s - a*b*pow(s,2.0) - 1.0/3.0*(2.0*a*c+pow(b,2.0))*pow(s,3.0);
    f += -0.5*b*c*pow(s,4.0) - 0.2*pow(c,2.0)*pow(s,5.0);
  }
  return f;
}




float LineIntegral(const std::array<Float_t, 3>& a, const std::array<Float_t, 3>& b, Phantom& ob, int mat_data, float ToF)
{
	float away = 1.0e6;
	int em_radial_max = 4;
  /*
    integral of object ob from point a to point b
    path1 is bigger than path2
  */
  int debug=0;
  std::array<Float_t, 3> del,delt,ap,ar,at,br,bp,bt;
  double dels,determ,c1,c2,c0,m;
  float l2inf,l2sup,value,value_grad[3];
  float phi_trunc,theta_trunc,linf_trunc,lsup_trunc;
  double path1,path2,z1,z2;
  int i;
  float as,bs,cs;
  double x1,x2;
  float degree_to_rad = 3.141592/180;

  l2inf = ob.linf;
  l2sup = ob.lsup;
  theta_trunc = ob.trunc_polar*degree_to_rad;
  phi_trunc = ob.trunc_azimut*degree_to_rad;
  linf_trunc = ob.trunc_min;
  lsup_trunc = ob.trunc_max;
  value = 0.0;
  for (i=0;i<3;i++) value_grad[i] = 0.0;
  if (mat_data == EmissionSimul) {
    value = ob.em_value;
    if (ob.em_slope != 0.0) {
      value_grad[0] = ob.em_slope *sin(ob.em_polar)*cos(ob.em_azimut);
      value_grad[1] = ob.em_slope *sin(ob.em_polar)*sin(ob.em_azimut);
      value_grad[2] = ob.em_slope *cos(ob.em_polar);
    }
  } else if (mat_data == AttenuationSimul || mat_data == TransmissionSimul) value = ob.mu_value;
  if (value == 0.0) return 0;

  /*shift coordinates*/
  for(i=0;i<3;++i) {
     ap[i] = (a[i] - (ob.center[i]));
     bp[i] = (b[i] - (ob.center[i]));
  }

  /*rotate only if needed to speed-up*/
  if(((ob.theta) == 0) && ((ob.phi) == 0)) {
     for(i=0;i<3;++i) {
	ar[i] = ap[i];
	br[i] = bp[i];
     }
  } else {  /*rotation needed*/
     ar[1] = -ob.sp*ap[0] + ob.cp*ap[1];
     ar[0] = -ob.st*ap[2] + ob.ct*(ob.cp*ap[0] + ob.sp*ap[1]);
     ar[2] =  ob.ct*ap[2] + ob.st*(ob.cp*ap[0] + ob.sp*ap[1]);
     br[1] = -ob.sp*bp[0] + ob.cp*bp[1];
     br[0] = -ob.st*bp[2] + ob.ct*(ob.cp*bp[0] + ob.sp*bp[1]);
     br[2] =  ob.ct*bp[2] + ob.st*(ob.cp*bp[0] + ob.sp*bp[1]);
  }

  dels = 0;
  for(i=0;i<3;++i) {
     del[i] = br[i]-ar[i];
     dels += del[i]*del[i];
  }
  dels = sqrt(dels);
  /*get the coefficient of the equation for the intersection of the
   line ab with the ellipsoid*/
  c0 = -1;
  c1 = 0;
  c2 = 0;
  for(i=0;i<3;++i) {
     if (ob.half_axis[i] > 0) {
       c0 += ar[i]*ar[i]  /(ob.half_axis[i]*ob.half_axis[i]);
       c1 += ar[i]*del[i] /(ob.half_axis[i]*ob.half_axis[i]);
       c2 += del[i]*del[i]/(ob.half_axis[i]*ob.half_axis[i]);
     }
  }
  determ = c1*c1 - c0*c2;
  if(determ <= 0) return 0;
  else {
    path1 = (-c1 + sqrt(determ))/c2;
    path2 = (-c1 - sqrt(determ))/c2;
    if ((lsup_trunc < away) || (linf_trunc > -away)) {
      at[2] = (ar[0]*cos(phi_trunc) + ar[1]*sin(phi_trunc))*sin(theta_trunc) + ar[2]*cos(theta_trunc);
      bt[2] = (br[0]*cos(phi_trunc) + br[1]*sin(phi_trunc))*sin(theta_trunc) + br[2]*cos(theta_trunc);
      delt[2] = bt[2]-at[2];
      z1    = at[2] + path1*delt[2];
      z2    = at[2] + path2*delt[2];
      if ((z1 > lsup_trunc) && (z2 > lsup_trunc)) return 0;
      if ((z1 < linf_trunc) && (z2 < linf_trunc)) return 0;
      if (z1 > lsup_trunc) path1 = (lsup_trunc - at[2])/delt[2];
      if (z2 < linf_trunc) path2 = (linf_trunc - at[2])/delt[2];
      if (z1 < linf_trunc) path1 = (linf_trunc - at[2])/delt[2];
      if (z2 > lsup_trunc) path2 = (lsup_trunc - at[2])/delt[2];
    }
    if (l2sup < away || l2inf > -away) {
      z1    = ar[2] + path1*del[2];
      z2    = ar[2] + path2*del[2];
      if ((z1 > l2sup) && (z2 > l2sup)) return 0;
      if ((z1 < l2inf) && (z2 < l2inf)) return 0;
      if (z1 > l2sup) path1 = (l2sup - ar[2])/del[2];
      if (z2 < l2inf) path2 = (l2inf - ar[2])/del[2];
      if (z1 < l2inf) path1 = (l2inf - ar[2])/del[2];
      if (z2 > l2sup) path2 = (l2sup - ar[2])/del[2];
    }
    if (mat_data == EmissionSimul) {
      if (ob.em_slope != 0.0) {
	m = 0.0;
	for (i=0;i<3;i++) m += value_grad[i]*(ar[i] + 0.5*del[i]*(path2+path1));
	return dels*fabs(path1-path2)*(value+m);
      } else if (ob.em_radial) {
	as = 0.0;
	bs = 0.0;
	cs = 0.0;
	for (i=0;i<3;i++) if (ob.half_axis[i] > 0.0) {
	  as += pow(ar[i]/ob.half_axis[i],2.0);
	  bs += (2.0/dels)*(ar[i]*del[i]/pow(ob.half_axis[i],2.0));
	  cs += pow(del[i]/(dels*ob.half_axis[i]),2.0);
	}
	if (path1 > path2) return
	  value*(emRad(as,bs,cs,path1*dels,ob.em_radial)-emRad(as,bs,cs,path2*dels,ob.em_radial));
	else return
	  value*(emRad(as,bs,cs,path2*dels,ob.em_radial)-emRad(as,bs,cs,path1*dels,ob.em_radial));
      } else {
        // Uniform objects
        if (path1 < path2) {
          if (fabs(path1-path2) > 1.0e-7) {
            printf("oAnalyticalPhantom::LineIntegral> Error: path1 = %0g < path2 = %0g, difference = %0g > 1.0e-7\n",path1,path2,path1-path2);
            //exit(EXIT_FAILURE);
          } else {
            printf("oAnalyticalPhantom::LineIntegral> Warning: path1 = %0g < path2 = %0g, difference = %0g > 0.0 (but < 1.0e-7)\n",path1,path2,path1-path2);
            path1 = path2;
          }
        } else if (path2 < 0) {
          if (fabs(path2) > 1.0e-7) {
            printf("oAnalyticalPhantom::LineIntegral> Error: path2 = %0g < -1.0e-7\n",path2);
            //exit(EXIT_FAILURE);
          } else {
            printf("oAnalyticalPhantom::LineIntegral> Warning: path2 = %0g < 0 (but > -1.0e-7)\n",path2);
            path2 = 0.0;
          }
        }
        /*
        if (m_TOFInfoFlag) {
          x1 = erf((dels*(path1-0.5)-ToF)/(m_sigmaToF*sqrt(2.0)));
          x2 = erf((dels*(path2-0.5)-ToF)/(m_sigmaToF*sqrt(2.0)));
          if (x1 < x2) {
            printf(" line_integral> Error: x1 = %0g < x2 = %0g, difference = %0g\n",x1,x2,x1-x2);
            //if (fabs(x1-x2) > 1.0e-7) exit(EXIT_FAILURE);
          }
          if (debug) {
            printf("                       ToF = %0g\n",ToF+0.5*dels);
          }
          return value*0.5*fabs(x1-x2);

        }
        */
      //else {


          if (debug) {
            printf(" line_integral> DEBUG\n");
            printf("                       detector 1 = <%0g ; %0g ; %0g>\n",ar[0],ar[1],ar[2]);
            printf("                       detector 2 = <%0g ; %0g ; %0g>\n",br[0],br[1],br[2]);
            //printf("                       dels = %0g\n",dels);
            printf("                       path1 = %0g\n",path1);
            printf("                       patt2 = %0g\n",path2);
            printf("                       length = %0g [mm]\n",dels*fabs(path1-path2));
          }
          if (path1 < path2) {
            printf(" line_integral> Error: path1 = %0g < path2 = %0g, difference = %0g\n",path1,path2,path1-path2);
            //if (fabs(path1-path2) > 1.0e-7) exit(EXIT_FAILURE);
          } else if (path2 < 0) {
            printf(" line_integral> Error: path2 = %0g < 0\n",path2);
            //if (fabs(path2) > 1.0e-7) exit(EXIT_FAILURE);
          }
          return dels*fabs(path1-path2)*value;
        //}
      }
    } else {
      return dels*fabs(path1-path2)*value;
    }
  }
}

Sinogram ConvertToSinogram(TVector3& gPos1, TVector3& gPos2) {


    double dx, dy, dz, d;
    float  r, phi, tan_theta, z;
    float pi=(float)M_PI;

    dz = gPos2.Z() - gPos1.Z();
    dy = gPos2.Y() - gPos1.Y();
    dx = gPos2.X() - gPos1.X();
    d = sqrt( dx*dx + dy*dy);
    r = (float)((gPos1.Y()*gPos2.X()-gPos1.X()*gPos2.Y() )/d);
//    if (dx == 0.0) phi = 0.0; else
    phi = (float)(atan2( dx, dy));
    tan_theta = (float)(dz / d);
    z = (float)(gPos1.Z()+(gPos1.X()*dx+gPos1.Y()*dy)*dz/(d*d));
    if (phi < 0.0) {
        phi += pi;
        r *= -1.0;
        tan_theta *= -1.0;
    }
    if (phi == pi){
        phi=0.0;
        r *= -1.0;
        tan_theta *= -1.0;
    }

    double cos_theta = 1.0/sqrt(1.0+tan_theta*tan_theta);

	Sinogram sinogram;


	float S = (gPos1.Z() + gPos2.Z()) / 2.0;


    // Store results
    sinogram.R = r;
    sinogram.S = S;
    sinogram.Phi = phi;
    sinogram.Cos_theta = cos_theta;

    return sinogram;
}




std::vector<std::string> expandWildcard(const std::string &pattern) {
    glob_t glob_result;
    std::vector<std::string> files;

    glob(pattern.c_str(), GLOB_TILDE, nullptr, &glob_result);
    for (size_t i = 0; i < glob_result.gl_pathc; ++i)
        files.emplace_back(glob_result.gl_pathv[i]);
    globfree(&glob_result);

    return files;
}


void writeNormEntryNormMatrixFile(std::ofstream& out, uint32_t c1, uint32_t c2, float n)
{
	out.write(reinterpret_cast<const char*>(&n),  sizeof(float));
    out.write(reinterpret_cast<const char*>(&c1), sizeof(uint32_t));
    out.write(reinterpret_cast<const char*>(&c2), sizeof(uint32_t));

}


void writeCdhHeaderNormMatrixFile(const std::string& outputDir,const std::string& outputMatrixFileName,
                    const std::string& scannerName,
                    size_t nEvents)
{
    std::ofstream header(outputDir+outputMatrixFileName+"_df.Cdh");
    header << "Data filename: " << outputMatrixFileName+"_df.Cdf" << "\n";
    header << "Number of events: " << nEvents << "\n";
    header << "Data mode: normalization\n";
    header << "Data type: PET\n";
    header << "Scanner name: " << scannerName << "\n";
    header << "Calibration factor: 1\n";
    header << "Isotope: unknown\n";
    header << "Normalization correction flag: 1";
    //header << "!END OF INTERFILE :=\n";
    header.close();
}

TVector3 applyRotation(const TVector3& pos, double angle) {
    // Rotation around Z-axis by given angle


    double cosA = std::cos(angle);
    double sinA = std::sin(angle);

    double R = pow(pos.X()*pos.X()+pos.Y()*pos.Y(),0.5);
    double xRot = pos.X() * cosA - pos.Y() * sinA;
    double yRot =  pos.X() * sinA + pos.Y() * cosA;
    double zRot = pos.Z();

    return TVector3(xRot, yRot, zRot);
}

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
               uint32_t nSubmodulesAxial) {

  // Compute variable displacements according to TODOs:
  // - layer thickness = crystalDepth / nLayers
  // - crystal transaxial pitch = transAxialSize / nCrystalsTransaxial
  // - submodule axial pitch = axialSize / nSubmodulesAxial
  double layerThickness = (nLayers > 0) ? (crystalDepth / static_cast<double>(nLayers)) : 0.0;
  double crystalPitch = (nCrystalsTransaxial > 0) ? (transAxialSize / static_cast<double>(nCrystalsTransaxial)) : 0.0;
  double submodulePitch = (nSubmodulesAxial > 0) ? (axialSize / static_cast<double>(nSubmodulesAxial)) : 0.0;

  double x = x0 + layerID * layerThickness;
  double y = y0 + crystalID * crystalPitch;
  double z = z0 + submoduleID * submodulePitch;

  TVector3 localPos(x, y, z);

  // Apply rotation around Z (to place it in the correct angular sector)
  TVector3 rotated = applyRotation(localPos, rsectorID * deltaPhi);

  // Preserve original module Z-shift behavior (use existing constant)
  rotated.SetZ(rotated.Z() + moduleID * 63.0);

  return rotated;
}

double meanVector(const std::vector<double>& v)
{
  if (v.empty()) return 0.0;
  double sum = 0.0;
  size_t count = 0;
  for (double x : v) {
    if (std::isfinite(x)) { sum += x; ++count; }
  }
  return count ? (sum / double(count)) : 0.0;
}
double meanMatrix(const std::vector<std::vector<double>>& M)
{
  double sum = 0.0;
  size_t count = 0;
  for (const auto& row : M) {
    for (double x : row) {
      if (std::isfinite(x)) { sum += x; ++count; }
    }
  }
  return count ? sum / double(count) : 0.0;
}



double scaledMeanMatrix(
    const std::vector<std::vector<double>>& mat,
    const std::vector<double>& weights
){
    const int N = mat.size();
    if (N == 0) return 0.0;

    if (weights.size() != static_cast<size_t>(N)) {
        throw std::runtime_error("scaledMeanMatrix(): weight size mismatch");
    }

    const int M = mat[0].size();

    double sum = 0.0;

    for (int i = 0; i < N; ++i) {
        const double w = weights[i];

        if (mat[i].size() != static_cast<size_t>(M)) {
            throw std::runtime_error("scaledMeanMatrix(): jagged matrix");
        }

        for (int j = 0; j < M; ++j) {
            sum += mat[i][j] * w;
        }
    }

    return sum / (N * M);
}



inline std::ofstream openCSV(const std::string& filename) {
    std::ofstream f(filename);
    if (!f.is_open())
        throw std::runtime_error("ERROR: Cannot open CSV file: " + filename);
    return f;
}

inline int mod(int x, int m)
{
    return ( (x % m) + m ) % m;
}

