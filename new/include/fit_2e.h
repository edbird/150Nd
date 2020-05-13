#ifndef FIT_H
#define FIT_H

#include "TH1F.h"
#include "TFractionFitter.h"
#include "TArray.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TPaveStats.h"
#include "TPaveText.h"

#include "THStack.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2.h"
#include "TStyle.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TMath.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include "math.h"
#include <string>
#include "TLegend.h"
#include "TColor.h"
#include "TVector.h"
#include <vector>


// TODO: enable this at some point
//#include "InputFiles.h"


//#include "fit_2e.h"
#include "InputNumberDef.h"
#include "InputFileDef.h"
#include "InputNameDef.h"
#include "InputColorDef.h"


#define ND_ON 1
#define NEIGHBOURS_ON 1
#define INTERNALS_ON 1
#define EXTERNALS_ON 1
#define RADON_ON 1


//Int_t nIntBkg, nExtBkg, nSignals;
// I don't remember what this is for
Int_t nHistograms;

//static const long TotalTime = 167629292; // P1+P2
const double TotalTime = 167629292.; // P1+P2

//long AcceptedTime[3];
double AcceptedTime[3];

// set to 1 to enable truth file for 150Nd
#define TRUTH_ENABLE 1

// MODE 1 = 2e
// MODE 2 = 2eNg_29Sep2015
//#define MODE 2
// the path of my nd150 MC with truth information
// /unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/nd150/nd150_rot_2b2n_m4/Nd150_2eNg_output_truth.root

// enable/disable filepath changing for external disk at home
// default 0 is path at UCL server
#define MODE_AT_HOME 1
// enable/disable reading from ramdisk
#define MODE_AT_HOME_RAMDISK 0

#if MODE == 1
    //TString filePath="/unix/nemo3/users/sblot/Nd150Analysis/newAnalysis/2eNg_29Sep2015/";

    // 2e directory
    #if MODE_AT_HOME
        #if MODE_AT_HOME_RAMDISK
            TString filePath = "/mnt/ramdisknd150/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
        #else
            TString filePath = "/media/ecb/Maxtor/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
        #endif
    #else
        TString filePath = "/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
    #endif
#elif MODE == 2
    //TString filePath="/unix/nemo3/users/sblot/Nd150Analysis/newAnalysis/2e/";
    
    // 2e directory
    //TString filePath="/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
    
    // 2eNg directory
    #if MODE_AT_HOME
        #if MODE_AT_HOME_RAMDISK
            TString filePath = "/mnt/ramdisknd150/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2eNg_29Sep2015/";
        #else
            TString filePath = "/media/ecb/Maxtor/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2eNg_29Sep2015/";
        #endif
    #else
        TString filePath = "/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2eNg_29Sep2015/";
    #endif
    
    // TODO email summer and ask which she thinks should we be using
    // TODO: really don't remember if it should be 2e or 2eNg_29Sep2015
#endif

// 2019-01-16 I changed the to phase 1
// 2019-01-22: Phase can be either 0 for "Phase 1" or 1 for "Phase 2"
// 2020-03-11: Think that Phase is "1" or "2", corresponding to value of
// thePhase of 0 or 1 respectively
Int_t thePhase = 1;
TString Phase = "2";
//Phase.Form("%i", thePhase);


// MC / data info
//
// Group Name:
// Activities file - activities are pulled from this file for scaling
//
// Rn222:
// radonActivities.txt
//
// Rn220:
// radonActivities.txt
//
// Note: There are no duplicate names in the list of MC samples for Rn222 and Rn220
//
// External:
// externalActivities.txt
//
// Internal:
// internalActivities.txt
//
// Nd150:
// internalActivities2.txt
//
// Neighbours:
// internalActivities2.txt

// TODO: some of these use links from the 2eNg... directory to the 2e
// directory
// check if the 2e directory actually has any significant number of events
// passing cuts, if not then remove them and remove the links


// processed
double ExternalBkgActivity[2][nExternalBkgs];
double ExternalBkgEfficiency[2][nExternalBkgs];
double ExternalBkgPassCut[2][nExternalBkgs];
double ExternalBkgNGenMC[nExternalBkgs];

double Rn222BkgActivity[2][nRn222Bkgs];
double Rn222BkgEfficiency[2][nRn222Bkgs];
double Rn222BkgPassCut[2][nRn222Bkgs];
double Rn222BkgNGenMC[nRn222Bkgs];

double Rn220BkgActivity[2][nRn220Bkgs];
double Rn220BkgEfficiency[2][nRn220Bkgs];
double Rn220BkgPassCut[2][nRn220Bkgs];
double Rn220BkgNGenMC[nRn220Bkgs];

double InternalBkgActivity[2][nInternalBkgs];
double InternalBkgEfficiency[2][nInternalBkgs];
double InternalBkgPassCut[2][nInternalBkgs];
double InternalBkgNGenMC[nInternalBkgs];

double NeighbourActivity[2][nNeighbours];
double NeighbourEfficiency[2][nNeighbours];
double NeighbourPassCut[2][nNeighbours];
double NeighbourNGenMC[nNeighbours];

double Nd150Activity[2][nNd150Samples];
double Nd150Efficiency[2][nNd150Samples];
double Nd150PassCut[2][nNd150Samples];
double Nd150NGenMC[nNd150Samples];


// raw
double ExternalBkgActivity_rawdata[2][nExternalBkgs];
double ExternalBkgEfficiency_rawdata[2][nExternalBkgs];
double ExternalBkgPassCut_rawdata[2][nExternalBkgs];
double ExternalBkgNGenMC_rawdata[nExternalBkgs];

double Rn222BkgActivity_rawdata[2][nRn222Bkgs];
double Rn222BkgEfficiency_rawdata[2][nRn222Bkgs];
double Rn222BkgPassCut_rawdata[2][nRn222Bkgs];
double Rn222BkgNGenMC_rawdata[nRn222Bkgs];

double Rn220BkgActivity_rawdata[2][nRn220Bkgs];
double Rn220BkgEfficiency_rawdata[2][nRn220Bkgs];
double Rn220BkgPassCut_rawdata[2][nRn220Bkgs];
double Rn220BkgNGenMC_rawdata[nRn220Bkgs];

double InternalBkgActivity_rawdata[2][nInternalBkgs];
double InternalBkgEfficiency_rawdata[2][nInternalBkgs];
double InternalBkgPassCut_rawdata[2][nInternalBkgs];
double InternalBkgNGenMC_rawdata[nInternalBkgs];

double NeighbourActivity_rawdata[2][nNeighbours];
double NeighbourEfficiency_rawdata[2][nNeighbours];
double NeighbourPassCut_rawdata[2][nNeighbours];
double NeighbourNGenMC_rawdata[nNeighbours];

double Nd150Activity_rawdata[2][nNd150Samples];
double Nd150Efficiency_rawdata[2][nNd150Samples];
double Nd150PassCut_rawdata[2][nNd150Samples];
double Nd150NGenMC_rawdata[nNd150Samples];


///////////////////////////////////////////////////////////////////////////////




//its annoying to have to keep writing the histogram names. If we always work in the right order, should be OK.
static const int numHistograms = 49; //40; // 44;
TString histogramNames[numHistograms] = {
    "hRun_",                            // 00
    "hNElectrons_",                     // 01
    "hTotalE_",                         // 02
    "hEeMax_",                          // 03
    "hElectronLengthMax_",              // 04

    "hVertexZMax_",
    "hVertexSectorMax_",
    "hVertexRMax_",
    "hElectronFirstGgMax_",
    "hElectronLastGgMax_",              // 09

    "hVertexMinDistPromptGgMax_",
    "hElectronLDCorrMax_",
    "hEeMin_",
    "hElectronLengthMin_",
    "hVertexZMin_",                     // 14

    "hVertexSectorMin_",
    "hVertexRMin_",
    "hElectronFirstGgMin_",
    "hElectronLastGgMin_",
    "hVertexMinDistPromptGgMin_",       // 19

    "hElectronLDCorrMin_",
    "hNLowEGammas_",
    "hLowEGammaEnergy_",
    "hSummedLowEGammaE_",
    "hLowEMinDistPromptGg_",            // 24

    "hInternalPullee_",
    "hInternalProbee_",
    "hExternalPullee_",
    "hExternalProbee_",
    "hCosee_",                          // 29

    "hCoseeWeighted_",
    "hVertexDZ_",
    "hVertexDR_",
    "hVertexDRPhi_",
    "hNAPromptGgHits_",                 // 34

    "hNAfterCuts_",
    "hVertexZSecMax_",
    "hVertexZSecMin_",
    "hEeMaxVEeMin_",
    "hNAPromptGgHitsDist2VertexMin_",   // 39

    "hTrackSignMax_",
    "hTrackSignMin_",
    "hnGammaClusters_",
    "hnInCluster_",
    "hclusterHitEnergy_",               // 44

    "hclusterHitEnergyMin_",
    "hclusterHitEnergyMax_",
    "hnLowEnergyHits_",                  // 47
    "hclusterEnergy_"

};

// TODO: what to do about hMCElectronEnergy and hMCTrueElectronEnergy

Int_t histogramDrawFlag[numHistograms] =
{
    1,    // hRun_
    0,    // hNElectrons_
    1,    // hTotalE_
    1,    // hEeMax_
    1,    // hElectronLengthMax_

    1,    // hVertexZMax_
    1,    // hVertexSectorMax_
    1,    // hVertexRMax_
    1,    // hElectronFirstGgMax_
    1,    // hElectronLastGgMax_

    1,//1,    // hVertexMinDistPromptGgMax_
    1,    // hElectronLDCorrMax_
    1,    // hEeMin_
    1,    // hElectronLengthMin_
    1,    // hVertexZMin_

    1,    // hVertexSectorMin_
    1,    // hVertexRMin_
    1,    // hElectronFirstGgMin_
    1,    // hElectronLastGgMin_
    1,//1,    // hVertexMinDistPromptGgMin_

    1,    // hElectronLDCorrMin_
    1,    // hNLowEGammas_
    1,    // hLowEGammaEnergy_
    1,    // hSummedLowEGammaE_
    1,//1,    // hLowEMinDistPromptGg_

    1,    // hInternalPullee_
    1,    // hInternalProbee_
    1,    // hExternalPullee_
    1,    // hExternalProbee_
    0,    // hCosee_

    0,    // hCoseeWeighted_
    1,//1,    // hVertexDZ_
    1,//1,    // hVertexDR_
    1,//1,    // hVertexDRPhi_
    1,//1,    // hNAPromptGgHits_

    1,    // hNAfterCuts_
    0,    // hVertexZSecMax_
    0,    // hVertexZSecMin_
    0,    // hEeMaxVEeMin_
    1,//1    // hNAPromptGgHitsDist2VertexMin_

    0,    //hTrackSignMax_
    0,    //hTrackSignMin_
    1,      // hnGammaClusters_
    1,      // hnInCluster_
    1,      // hclusterHitEnergy_

    1,      // hclusterHitEnergyMin_
    1,      // hclusterHitEnergyMax_
    1,      // hnLowEnergyHits_
    1       // hClusterEnergy_
};

// draw raw data histograms as well?
Int_t histogramDrawFlag_rawdata = 1;


// map for data structure for raw data
std::map<TString, TH1*> histogramPointers_rawdata;
// map for data structure for processed data
std::map<TString, TH1*> histogramPointers;

#endif
