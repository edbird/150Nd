#ifndef FIT_H
#define FIT_H

#include "TH1D.h"
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
#include "TH1.h"
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

#if 0
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
            //TString filePath = "/media/ecb/Maxtor/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
            //TString filePath = "/media/ecb/backup/users/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
            // change to using md0 over network
            //TString filePath = "/mnt/md0/users/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
            TString filePath = "/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
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
            //TString filePath = "/media/ecb/Maxtor/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2eNg_29Sep2015/";
            //TString filePath = "/media/ecb/backup/users/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2eNg_29Sep2015/";
            // change to using md0 over network
            //TString filePath = "/mnt/md0/users/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2eNg_29Sep2015/";
            TString filePath = "/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2eNg_29Sep2015/";
        #endif
    #else
        TString filePath = "/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2eNg_29Sep2015/";
    #endif
    
    // TODO email summer and ask which she thinks should we be using
    // TODO: really don't remember if it should be 2e or 2eNg_29Sep2015
#endif
#endif

TString filePath = "/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
//TString filePath = "/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2eNg_29Sep2015/";

// 2020-01-16 I changed the to phase 1
// 2020-01-22: Phase can be either 0 for "Phase 1" or 1 for "Phase 2"
// 2020-03-11: Think that Phase is "1" or "2", corresponding to value of
// thePhase of 0 or 1 respectively
Int_t nocutonphase = 0;
// delete this comment
Int_t thePhase = 0;
TString Phase = "1";
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
static const int numHistograms = 61; //59; // 57; //49; //40; // 44;
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
    "hnLowEnergyHits_",
    "hclusterEnergy_",
    "hFoilSide_",                       //49

    "hTrueVertexR_",    
    "hFoilSideTrueVertexR_",
    "hTrueVertexLayer_",
    "hSingleEnergy_",
    "hLowEnergy_",                      // 54

    "hHighEnergy_",                     
    "hHighLowEnergy_",
    "hTotalClusterEnergy_",
    "hHotSpotLocation_",
    "hEnergySum_",                      //59

    "hEnergyDiff_"

};

// TODO: what to do about hMCElectronEnergy and hMCTrueElectronEnergy

Int_t histogramDrawFlag[numHistograms] =
{
    0,    // hRun_
    0,    // hNElectrons_
    1,    // hTotalE_
    1,    // hEeMax_
    0,    // hElectronLengthMax_

    0,    // hVertexZMax_
    0,    // hVertexSectorMax_
    0,    // hVertexRMax_
    0,    // hElectronFirstGgMax_
    0,    // hElectronLastGgMax_

    0,//1,    // hVertexMinDistPromptGgMax_
    0,    // hElectronLDCorrMax_
    1,    // hEeMin_
    0,    // hElectronLengthMin_
    0,    // hVertexZMin_

    0,    // hVertexSectorMin_
    0,    // hVertexRMin_
    0,    // hElectronFirstGgMin_
    0,    // hElectronLastGgMin_
    0,//1,    // hVertexMinDistPromptGgMin_

    0,    // hElectronLDCorrMin_
    0,    // hNLowEGammas_
    0,    // hLowEGammaEnergy_
    0,    // hSummedLowEGammaE_
    0,//1,    // hLowEMinDistPromptGg_

    0,    // hInternalPullee_
    0,    // hInternalProbee_
    0,    // hExternalPullee_
    0,    // hExternalProbee_
    0,    // hCosee_

    0,    // hCoseeWeighted_
    0,//1,    // hVertexDZ_
    0,//1,    // hVertexDR_
    0,//1,    // hVertexDRPhi_
    0,//1,    // hNAPromptGgHits_

    0,    // hNAfterCuts_
    0,    // hVertexZSecMax_
    0,    // hVertexZSecMin_
    0,    // hEeMaxVEeMin_
    0,//1    // hNAPromptGgHitsDist2VertexMin_

    0,    //hTrackSignMax_
    0,    //hTrackSignMin_
    0,      // hnGammaClusters_
    0,      // hnInCluster_
    0,      // hclusterHitEnergy_

    0,      // hclusterHitEnergyMin_
    0,      // hclusterHitEnergyMax_
    0,      // hnLowEnergyHits_
    0,      // hClusterEnergy_
    0,      // hFoilSide_

    0,      // hTrueVertexR_
    0,      // hFoilSideTrueVertexR_
    0,      // hTrueVertexLayer_
    1,      // hSingleEnergy_
    1,      // hLowEnergy_

    1,      // hHighEnergy_
    0,      // hHighLowEnergy_
    0,      // hTotalClusterEnergy_
    0,      // hHotSpotLocation_
    1,      // hEnergySum_

    1       // hEnergyDiff_
};

// draw raw data histograms as well?
Int_t histogramDrawFlag_rawdata = 1;
#define RAWENABLE 0

// map for data structure for raw data
//std::map<TString, TH1*> histogramPointers_rawdata;
// map for data structure for processed data
//std::map<TString, TH1*> histogramPointers;

const int N_HOTSPOT_ELIPSE = 13;
double hotspot_elipse[13][4] = {
    {5.799247,  107.017857, 0.07 / 2.0,        5.0 / 2.0},
    {5.851,     7.65,       0.029573 / 2.0,    8.81 / 2.0},
    {5.765,     88.0,       0.0191346 / 2.0,   4.78448 / 2.0},
    {5.792,     80.98,      0.0143 / 2.0,      3.5846 / 2.0},
    {5.744,     54.5,       0.026385 / 2.0,    6.92231 / 2.0},
    {5.83,      48.,        0.042556 / 2.0,    10.2898 / 2.0},
    {5.805,     37.5,       0.034152 / 2.0,    12.4745 / 2.0},
    {5.76,      15.44,      0.0111 / 2.0,      5.1315 / 2.0},
    {5.785,     1.0,        0.0212464 / 2.0,   4.24026 / 2.0},
    {5.7388,    -12.829,    0.01413 / 2.0,     6.8421 / 2.0},
    {5.745,     -25.0,      0.029942 / 2.0,    6.2256 / 2.0},
    {5.798,     -80.955,    0.030184 / 2.0,    12.471 / 2.0},
    {5.8045,    -113.702,   0.023622 / 2.0,    6.852 / 2.0}
};


#endif
