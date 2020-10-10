#include "TH1F.h"
#include "TFractionFitter.h"
#include "TArray.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TROOT.h"

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
#include <vector>

#include "../include/fit_2e.h"

//-----------------------------------------
//    Functions
//----------------------------------------
void loadFiles();

void  makeHistograms(TString thePath, TString histName);
//void  makeHistograms(TTree *tree, TString histName, Int_t isotope);
void  fitHistograms();
void drawPlots(TH1F *data,TH1F *mc, THStack *hs);
double getActErr(double Npass, double Ngen, double Ndata, double sf_err );
double getChi2(TH1F *data, TH1F* mc, Int_t &ndof);

void loadFiles()
{

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);

    //First, create a root file to hold all of the histograms
    TFile *myFile = TFile::Open("Nd150_2e_P" + Phase + ".root", "RECREATE");
    TDirectory *dir = gDirectory;
    TDirectory *dir_sh = dir->mkdir("singleHistos"); 
    dir_sh->mkdir("unscaled");
    myFile->Close();

    // Read in all files and make the cuts you want.
    //First we will read in the data
    std::cout << ">>>>> Reading in data file...";
    makeHistograms("betabeta/", DataFile);
    std::cout << "... DONE." << std::endl;

    //Now we'll read in all the backgrounds: externals, Rn222, Rn220, internals and neighbouring foils
    std::cout << std::endl;
    std::cout << ">>>>> Reading in external backgrounds..." << std::endl;
    for(int i = 0; i < nExternalBkgs; i++)
    {
        makeHistograms("externals/", ExternalBkgFiles[i]);
    }
    std::cout << "... DONE." << std::endl;

    std::cout << std::endl;
    std::cout << ">>>>> Reading in Rn220 backgrounds..." << std::endl;
    for(int i = 0; i < nRn220Bkgs; i++)
    {
        makeHistograms("externals/", Rn220BkgFiles[i]);
    }
    std::cout << "... DONE." << std::endl;

    std::cout << std::endl;
    std::cout << ">>>>> Reading in Rn222 backgrounds..." << std::endl;
    for(int i = 0; i < nRn222Bkgs; i++)
    {
        makeHistograms("externals/", Rn222BkgFiles[i]);
    }
    std::cout << "... DONE." << std::endl;

    std::cout << std::endl;
    std::cout << ">>>>> Reading in internal backgrounds..." << std::endl;
    for(int i = 0; i < nInternalBkgs; i++)
    {
        makeHistograms("internals/", InternalBkgFiles[i]);
    }
    std::cout << "... DONE." << std::endl;

    std::cout << std::endl;
    std::cout << ">>>>> Reading in Nd150 backgrounds..." << std::endl;
    for(int i = 0; i < nNd150Samples; i++)
    {
        makeHistograms("nd150/", Nd150Files[i]);
    }
    std::cout << "... DONE." << std::endl;

    std::cout << std::endl;
    std::cout << ">>>>> Reading in neighbouring foil backgrounds..." << std::endl;
    for(int i = 0; i < nNeighbours; i++)
    {
        makeHistograms("neighbours/", NeighbourFiles[i]);
    }
    std::cout << "... DONE." << std::endl;

}


void makeHistograms(TString thePath, TString histName)
{


    hRun                    = new TH1F("hRun_" + histName,
                                       "Phase " + Phase + " " + histName + " runs; run numbers",
                                       8000, 0, 10000); // 1000, 8000

    hNElectrons             = new TH1F("hNElectrons_" + histName,
                                       "Phase " + Phase + " " + histName + " N tracks; N tracks/event",
                                       20, -0.5, 19.5); // limits? all events are 2 for nd150

    hNAPromptGgHits         = new TH1F("hNAPromptGgHits_" + histName,
                                       "Phase " + Phase + " " + histName + " N unassoc. prompt gg hits; N unassoc prompt gg hits",
                                       20, -0.5, 19.5); // limits? blank histogram for nd150?

    hTotalE                 = new TH1F("hTotalE_" + histName,
                                       "Phase " + Phase + " " + histName + " total energy; #SigmaE_{e} (MeV)",
                                       50, 0.0, 4.0); // 0.0, 2.2

    hInternalPullee         = new TH1F("hInternalPullee_" + histName,
                                       "Phase " + Phase + " " + histName + " internal hypothesis; ee Pull",
                                       50, -40., 40.); // limits? all events within -10, 10 for nd150

    hInternalProbee         = new TH1F("hInternalProbee_" + histName,
                                       "Phase " + Phase + " " + histName + " internal hypothesis; ee Probability",
                                       50, 0.0, 1.); // limits?

    hExternalPullee         = new TH1F("hExternalPullee_" + histName,
                                       "Phase " + Phase + " " + histName + " external hypothesis; ee Pull",
                                       50, -40., 40.); // limits? all events within -30, 10 for nd150

    hExternalProbee         = new TH1F("hExternalProbee_" + histName,
                                       "Phase " + Phase + " " + histName + " external hypothesis; ee Probability",
                                       50, 0.0, 1.); // limits?

    hCosee                  = new TH1F("hCosee_" + histName,
                                       "Phase " + Phase + " " + histName + " angular correlation of electron tracks; cos(#theta)_{ee}",
                                       50, -1., 1.); // limits?

    hCoseeWeighted          = new TH1F("hCoseeWeighted_" + histName,
                                       "Phase " + Phase + " " + histName + "corrected angular correlation of electron tracks; cos(#theta)_{ee}",
                                       50, -1., 1.); // limits?

    hVertexDZ               = new TH1F("hVertexDZ_" + histName,
                                       "Phase " + Phase + " " + histName + "distance between recon vertices; #DeltaZ (cm)",
                                       50, -20., 20.); // limits? all events within -10, 10 for nd150

    hVertexDR               = new TH1F("hVertexDR_" + histName,
                                       "Phase " + Phase + " " + histName + "distance between recon vertices; #DeltaR (cm)",
                                       50, -1., 1.); // limits?

    hVertexDRPhi            = new TH1F("hVertexDRPhi_" + histName,
                                       "Phase " + Phase + " " + histName + "distance between recon vertices; #DeltaR#phi (cm*rad)",
                                       50,-10.,10.); // limits? all events witnin -10, 10 for nd150

    hEeMax                  = new TH1F("hEeMax_" + histName,
                                       "Phase " + Phase + " " + histName + " higher energy Ee; E_{e} (MeV)"              ,
                                       50, 0.0, 4.0); // limits ok

    hElectronLengthMax      = new TH1F("hElectronLengthMax_" + histName,
                                       "Phase " + Phase + " " + histName + " higher energy e^{-} track length;  track length (cm)",
                                       50, 0, 600); // limits?

    hVertexZMax             = new TH1F("hVertexZMax_" + histName,
                                       "Phase " + Phase + " " + histName + " higher energy vertex position;  Z (cm)"                  ,
                                       50, -120, 120); // limits?

    hVertexSectorMax        = new TH1F("hVertexSectorMax_" + histName,
                                       "Phase " + Phase + " " + histName + " higher energy vertex position;  sector"                 ,
                                       50, 5.7, 5.9); // limits?
                                       //50, 5.5, 6.1); // limits?

    hVertexRMax             = new TH1F("hVertexRMax_" + histName,
                                       "Phase " + Phase + " " + histName + " higher energy vertex position;  R (cm)"                  ,
                                       50, 154.7, 155.0); // limits? all events within 154.85, 154.9 for nd150

    hElectronFirstGgMax     = new TH1F("hElectronFirstGgMax_" + histName,
                                       "Phase " + Phase + " " + histName + " higher energy e^{-} first gg hit location;first gg layer",
                                       9, -0.5, 8.5); // all in 0 for nd150

    hElectronLastGgMax      = new TH1F("hElectronLastGgMax_" + histName,
                                       "Phase " + Phase + " " + histName + " higher energy e^{-} last gg hit location;last gg layer",
                                       9, -0.5, 8.5); // all in 7/8 for nd150

    hVertexMinDistPromptGgMax  = new TH1F("hVertexMinDistPromptGgMax_" + histName,
                                       "Phase " + Phase + " " + histName + " distance from higher energy vertex to closest, prompt unassoc. gg hit;(v_{x}, v_{y}, v_{z}) - (gg_{x},gg_{y},gg_{z}) (cm)",
                                       50, 0, 600); // blank?

    hElectronLDCorrMax      = new TH1F("hElectronLDCorrMax_" + histName,
                                       "Phase "+Phase+" "+histName+" higher energy e^{-} LD correction;LD corr", 
                                       50, 0.9,1.1); // limits? don't know what this histogram is


    hVertexZSecMax          = new TH2F("hVertexZSecMax_" + histName,
                                       "Phase "+Phase+" "+histName+" higher energy vertex location; sector; Z (cm)",
                                       100, 5.7, 5.9,
                                       //100, 5.5, 6.1,
                                       100, -120, 120); 

    hEeMin                  = new TH1F("hEeMin_" + histName,
                                       "Phase "+Phase+" "+histName+" lower energy Ee; E_{e} (MeV)",
                                       50, 0.0, 4);

    hElectronLengthMin      = new TH1F("hElectronLengthMin_" + histName,
                                       "Phase "+Phase+" "+histName+" lower energy e^{-} track length;  track length (cm)",
                                       50, 0, 600);

    hVertexZMin             = new TH1F("hVertexZMin_" + histName,
                                       "Phase "+Phase+" "+histName+" lower energy vertex position;  Z (cm)"                  ,
                                       50, -120, 120);

    hVertexSectorMin        = new TH1F("hVertexSectorMin_" + histName,
                                       "Phase "+Phase+" "+histName+" lower energy vertex position;  sector"                 ,
                                       50, 5.7, 5.9);

    hVertexRMin             = new TH1F("hVertexRMin_" + histName,
                                       "Phase "+Phase+" "+histName+" lower energy vertex position;  R (cm)"                  ,
                                       50, 154.7, 155.0);

    hElectronFirstGgMin     = new TH1F("hElectronFirstGgMin_" + histName,
                                       "Phase "+Phase+" "+histName+" lower energy e^{-} first gg hit location;first gg layer",
                                       9, -0.5, 8.5);

    hElectronLastGgMin      = new TH1F("hElectronLastGgMin_" + histName,
                                       "Phase "+Phase+" "+histName+" lower energy e^{-} last gg hit location;last gg layer",
                                       9, -0.5, 8.5);

    hVertexMinDistPromptGgMin   = new TH1F("hVertexMinDistPromptGgMin_" + histName,
                                           "Phase " + Phase + " " + histName + " distance from lower energy vertex to closest, prompt unassoc. gg hit;(v_{x},v_{y},v_{z}) - (gg_{x},gg_{y},gg_{z}) (cm)",
                                           50, 0, 600);

    hElectronLDCorrMin      = new TH1F("hElectronLDCorrMin_" + histName,
                                       "Phase "+Phase+" "+histName+" lower energy e^{-} LD correction;LD corr"            ,
                                       50, 0.9, 1.1);

    hVertexZSecMin          = new TH2F("hVertexZSecMin_" + histName,
                                       "Phase "+Phase+" "+histName+" lower energy vertex location; sector; Z (cm)"              ,
                                       150, 5.7, 5.9,
                                       100, -120, 120);


    hNLowEGammas            = new TH1F("hNLowEGammas_" + histName,
                                       "Phase "+Phase+" "+histName+" unassoc. scint hits E < 200 keV; N hits"     ,
                                       11, -0.5, 10.5);

    hLowEGammaEnergy        = new TH1F("hLowEGammaEnergy_" + histName,
                                       "Phase "+Phase+" "+histName+" unassoc. scint hits E < 200 keV; Energy (indiv. hits) (MeV)" ,
                                       50, 0, 0.2);

    hSummedLowEGammaE       = new TH1F("hSummedLowEGammaE_" + histName,
                                       "Phase "+Phase+" "+histName+" unassoc. scint hits E < 200 keV; Energy (indiv. hits) (MeV)" ,
                                       50, 0, 0.5);

    hLowEMinDistPromptGg    = new TH1F("hLowEMinDistPromptGg_" + histName,
                                       "Phase "+Phase+" "+histName+" distance from unassoc. calo hits to closest, prompt unassoc. gg hit;(v_{x},v_{y},v_{z}) - (gg_{x},gg_{y},gg_{z}) (cm)",
                                       100, 0, 600.);//all low E hits


    hEeMaxVEeMin = new TH2F("hEeMaxVEeMin_" + histName,
                                       "Phase " + Phase + " " + histName + " electron energies; Ee^{Max} (MeV); Ee^{min} (MeV)",
                                       50, 0, 4, 50, 0, 4);

    hNAPromptGgHitsDist2VertexMin = new TH1F("hNAPromptGgHitsDist2VertexMin_" + histName,
                                       "Phase " + Phase + " " + histName + " min vertex distance to unassociated prompt geiger hit; min vertex distance (cm); y label",
                                       50, 0., 100.);
                                        



    hNAfterCuts      = new TH1F("hNAfterCuts_" + histName,
                                       "Phase " + Phase + " " + histName + " Analysis Cut Flow; ; N events",
                                       50, 0, 50); //never know how many we'll make... :)

    hEffAfterCuts      = new TH1F("hEffAfterCuts_" + histName,
                                       "Phase " + Phase + " " + histName + " Analysis Cut Flow; ; efficiency",
                                       50, 0, 50); //never know how many we'll make... :)



    TFile *aFile;
    if((histName.CompareTo("nd150_rot_2b2n_m4") == 0)   ||
       (histName.CompareTo("nd150_rot_2n2b_m4") == 0))
    {
        aFile = TFile::Open(filePath+thePath+histName+"/Nd150_2eNg_output_truth.root");
    }
    else
    {
        aFile = TFile::Open(filePath+thePath+histName+"/Nd150_2eNg_output.root");
    }
    TTree *theTree = (TTree*)aFile->Get("Nd150_2eNg/Nd150_2eNg");

    // general
    int         event;
    int         run;
    int         runStatus; // TODO
    int         nElectrons; // TODO

    double      radonWeight;
    double      bi210Weight;
    int         foilSide;
    double      eventTime;
    double      trueVertexR;
    double      trueVertexZ;
    double      trueVertexSector;
    int         trueVertexLayer;
    // 2020-02-19 I am now confused. The truth info appears to be in all of my nd150_rot_? files ?
    // I'm not confused. This info is here, it was trueElectronEnergy which was missing

    double      electronEnergy[2];
    double      eTrackLength[2];
    int         electronSide[2];
    double      trackSign[2];
    //double      radiusOfCurvature[2]; // TODO - branch not in Nd150 tree
    double      electronMeasTime[2];
    double      electronDMeasTime[2];
    int         electronBlockType[2];


    double      eeInternalPull;
    double      eeInternalProb;
    double      eeExternalPull;
    double      eeExternalProb;
    double      cosee;
    double      cosee_weight;

    int         electronPMT[2];
    int         electronLDFlag[2];
    double      electronLDCorr[2];
    double      electronLDCorrErr[2];

    double      vertexZ[2];
    double      vertexSec[2];
    double      vertexR[2];

    bool        vertexInHotSpot[2];

    int         firstGgHitLayer[2];
    int         lastGgHitLayer[2];

    int         NAPromptGgHits; // TODO - max value is 20 ? EDIT: APPEARS TO BE 5
    int         NAPromptGgHitsSide[5];
    double      NAPromptGgHitsDist2Vertex[5]; // TODO check number 20 - TODO: NUMBER IS 5 NOT 20
    double      NAPromptGgHitsDist2Calo[5];
    
    
    int         nGammaClusters;
    int         nInCluster[10];                          
    double      clusterEnergy[10];
    double      clusterTimeSpan[10];


    //double clusterLength : clusterLength[nGammaClusters]/D                    *
    //double eMingInternalVertexPull : eMingInternalVertexPull[nGammaClusters]/D*
    //double eMingInternalVertexProb : eMingInternalVertexProb[nGammaClusters]/D*
    //double coseMing  : coseMing[nGammaClusters]/D                             *
    //double eMaxgInternalVertexPull : eMaxgInternalVertexPull[nGammaClusters]/D*
    //double eMaxgInternalVertexProb : eMaxgInternalVertexProb[nGammaClusters]/D*
    //double coseMaxg  : coseMaxg[nGammaClusters]/D                             *
    

    int         nTotalClusterHits; // TODO - check max value is 10 ?
    // TODO: MAX VALUE IS 20 NOT 10!
    double      clusterHitEnergy[20];
    int         clusterHitPMT[20];
    int         clusterHitLDFlag[20];
    double      clusterHitLDCorr[20];
    double      clusterHitLDCorrErr[20];
    double      clusterHitSec[20];
    double      clusterHitZ[20];

    double      trueElectronEnergy[2]; // TODO


    //bool        mcFlag;
    //radon related stuff

    //counters 
    //int         nGammas, nGammaClusters
    //nElectrons, nParticles;
    //double      electronCaloZ[2], electronCaloSector[2], electronCaloR[2];
    // re-enabled, but not in tree?
    //double      electronCaloZ[2];
    //double      electronCaloSector[2];
    //double      electronCaloR[2];

    //double      dist2VertexPromptGg0[20];
    //double      dist2VertexPromptGg1[20];

    // these are the old variable names
    //int         nLowEClusters;
    //double      lowEClusterE[10];
    //double      lowEClusterSector[10];
    //double      lowEClusterZ[10];
    //double      lowEClusterR[10];
    //double      lowEMinDistPromptGg[10];
    //int         lowEClusterPMTs[10];
    // these are the branches they were mapped to
    //theTree->SetBranchAddress("nTotalClusterHits"  , &nLowEClusters);
    //theTree->SetBranchAddress("clusterHitEnergy"   , lowEClusterE);
    //theTree->SetBranchAddress("clusterHitPMT"      , lowEClusterPMTs);
    //theTree->SetBranchAddress("clusterHitSec"      , lowEClusterSector);
    //theTree->SetBranchAddress("clusterHitZ"        , lowEClusterZ);
    // theTree->SetBranchAddress( "lowEClusterR"      , lowEClusterR );
    // theTree->SetBranchAddress( "lowEMinDistPromptGg", lowEMinDistPromptGg );


    theTree->SetBranchAddress("Event"               , &event);
    theTree->SetBranchAddress("Run"                 , &run);
    theTree->SetBranchAddress("runStatus"           , &runStatus);  
    theTree->SetBranchAddress("nElectrons"          , &nElectrons);  
    //theTree->SetBranchAddress( "nGammaClusters"   , &nGammas );  
    //  theTree->SetBranchAddress( "nParticles"     , &nParticles );  

    theTree->SetBranchAddress("radonWeight"         , &radonWeight);
    theTree->SetBranchAddress("bi210Weight"         , &bi210Weight);
    theTree->SetBranchAddress("foilSide"            , &foilSide);
    theTree->SetBranchAddress("eventTime"           , &eventTime);

    theTree->SetBranchAddress("trueVertexR"         , &trueVertexR);
    theTree->SetBranchAddress("trueVertexZ"         , &trueVertexZ);
    theTree->SetBranchAddress("trueVertexSector"    , &trueVertexSector);
    theTree->SetBranchAddress("trueVertexLayer"     , &trueVertexLayer);

    theTree->SetBranchAddress("electronEnergy"      , electronEnergy);
    theTree->SetBranchAddress("eTrackLength"        , eTrackLength);
    theTree->SetBranchAddress("electronSide"        , electronSide);
    theTree->SetBranchAddress("trackSign"           , trackSign);
    //theTree->SetBranchAddress("radiusOfCurvature" , radiusOfCurvature);
    theTree->SetBranchAddress("electronMeasTime"    , electronMeasTime);
    theTree->SetBranchAddress("electronDMeasTime"   , electronDMeasTime);
    theTree->SetBranchAddress("electronBlockType"   , electronBlockType);

    theTree->SetBranchAddress("internalPull"        , &eeInternalPull);
    theTree->SetBranchAddress("internalProb"        , &eeInternalProb);
    theTree->SetBranchAddress("externalPull"        , &eeExternalPull);
    theTree->SetBranchAddress("externalProb"        , &eeExternalProb);
    theTree->SetBranchAddress("cosee"               , &cosee);
    theTree->SetBranchAddress("cosee_weight"        , &cosee_weight);

    theTree->SetBranchAddress("electronPMT"         , electronPMT);
    theTree->SetBranchAddress("electronLDFlag"      , electronLDFlag);
    theTree->SetBranchAddress("electronLDCorr"      , electronLDCorr);
    theTree->SetBranchAddress("electronLDCorrErr"   , electronLDCorrErr);

    theTree->SetBranchAddress("vertexZ"             , vertexZ);
    theTree->SetBranchAddress("vertexSec"           , vertexSec);
    theTree->SetBranchAddress("vertexR"             , vertexR);

    theTree->SetBranchAddress("vertexInHotSpot"     , vertexInHotSpot);

    theTree->SetBranchAddress("firstGgHitLayer"     , firstGgHitLayer);
    theTree->SetBranchAddress("lastGgHitLayer"      , lastGgHitLayer);

    // number of unassociated prompt geiger hits
    theTree->SetBranchAddress("NAPromptGgHits"              , &NAPromptGgHits);
    theTree->SetBranchAddress("NAPromptGgHitsSide"          , NAPromptGgHitsSide);
    theTree->SetBranchAddress("NAPromptGgHitsDist2Vertex"   , NAPromptGgHitsDist2Vertex);
    theTree->SetBranchAddress("NAPromptGgHitsDist2Calo"     , NAPromptGgHitsDist2Calo);

    // not sure if these branches are in the tree?
    //theTree->SetBranchAddress("electronCaloZ"     ,     electronCaloZ);
    //theTree->SetBranchAddress("electronCaloR"     ,     electronCaloR);
    //theTree->SetBranchAddress("electronCaloSector",     electronCaloSector);
    //  theTree->SetBranchAddress( "electronPMT"       , electronPMT );

    theTree->SetBranchAddress("nGammaClusters"              , &nGammaClusters);
    theTree->SetBranchAddress("nInCluster"                 , &nInCluster);
    theTree->SetBranchAddress("clusterEnergy"               , clusterEnergy);
    theTree->SetBranchAddress("clusterTimeSpan"             , clusterTimeSpan);
    //double clusterLength : clusterLength[nGammaClusters]/D                    *
    //double eMingInternalVertexPull : eMingInternalVertexPull[nGammaClusters]/D*
    //double eMingInternalVertexProb : eMingInternalVertexProb[nGammaClusters]/D*
    //double coseMing  : coseMing[nGammaClusters]/D                             *
    //double eMaxgInternalVertexPull : eMaxgInternalVertexPull[nGammaClusters]/D*
    //double eMaxgInternalVertexProb : eMaxgInternalVertexProb[nGammaClusters]/D*
    //double coseMaxg  : coseMaxg[nGammaClusters]/D                             *
    //

    theTree->SetBranchAddress("nTotalClusterHits"           , &nTotalClusterHits);
    theTree->SetBranchAddress("clusterHitEnergy"            , clusterHitEnergy);
    theTree->SetBranchAddress("clusterHitPMT"               , clusterHitPMT);
    theTree->SetBranchAddress("clusterHitLDFlag"            , clusterHitLDFlag);
    theTree->SetBranchAddress("clusterHitLDCorr"            , clusterHitLDCorr);
    theTree->SetBranchAddress("clusterHitLDCorrErr"         , clusterHitLDCorrErr);
    theTree->SetBranchAddress("clusterHitSec"               , clusterHitSec);
    theTree->SetBranchAddress("clusterHitZ"                 , clusterHitZ);
    // NOTE: these were already removed by summer
    // theTree->SetBranchAddress( "lowEClusterR"      , lowEClusterR );
    // theTree->SetBranchAddress( "lowEMinDistPromptGg", lowEMinDistPromptGg );
    
    if((histName.CompareTo("nd150_rot_2b2n_m4") == 0)   ||
       (histName.CompareTo("nd150_rot_2n2b_m4") == 0))
    {
        theTree->SetBranchAddress("trueElectronEnergy"      , trueElectronEnergy);
    }
    else
    {
        trueElectronEnergy[0] = 0.;
        trueElectronEnergy[1] = 0.;
    }

    std::cout << "Processing: " << histName << std::endl;


    int nCuts = 20; // 26;
    Int_t cut_counter[nCuts];
    //Int_t cut_counter_index = 0;
    for ( int i = 0; i < nCuts; i++ )
    {
        cut_counter[i] = 0;
    }

    Long_t events = (Long_t)theTree->GetEntries();
    for(Long_t event_i = 0; event_i < events; event_i++)
    {
        bool cut = false;
        if(event_i % 1000 == 0)
        {
            std::cout << "\r Processing complete : " << 100 * event_i / events << "%" << std::flush;
        }
        theTree->GetEvent(event_i);


        // theTree->SetBranchAddress( "electronCaloZ"     , electronCaloZ );
        // theTree->SetBranchAddress( "electronCaloR"     , electronCaloR );
        //  theTree->SetBranchAddress( "electronCaloSector", electronCaloSector );
        //electronCaloZ[0] = 0;
        //electronCaloZ[1] = 0;
        //electronCaloR[0] = 0;
        //electronCaloR[1] = 0;
        //electronCaloSector[0] = 0;
        //electronCaloSector[1] = 0;
        
        //for(Int_t i = 0; i < 10; ++i)
        //{
            //lowEClusterR[i] = 0;
            //lowEMinDistPromptGg[i] = 0;
        //}

        double weight = 1.;
        if(histName.CompareTo("bi214_swire") == 0)
        {
            weight = radonWeight;
        }
        // TODO: due to secular equlibrium, the following MAY be required
        // (unsure at present time)
        else if(histName.CompareTo("pb214_swire") == 0)
        {
            weight = radonWeight;
        }
        else if(histName.CompareTo("bi210_swire") == 0)
        {
            weight = bi210Weight;
        }
        
        // TODO: re-implement this
        // http://nile.hep.utexas.edu/DocDB/ut-nemo/docs/0037/003716/005/position_paper.pdf
        // Page 20. Section on short-lived isotopes.
        // Top of page 31 has Bi events being weighted by half-life of Pb
        // (Doing a ctrl-f on weight works quite well, but those are the main two.)
        // Table 9.1
        // Some isotopes contributing to the background rates have relatively short half-lives
        // on the time scale of the NEMO-3 experiment. The events from these isotopes are
        // therefore weighted according to their known half-lives and the time of the event.
        // I assume that in order to do this, one must have a measurement of the activity of
        // these isotopes at a particular time, and then the weight is given by
        // A * exp( -k (t - t0) ), where A is the initial weight measured at time t0 (A=1), k is a decay rate
        // On Page 31 I find this "Events from 210Bi are weighted by the half life of 210Pb as described in Section 3.4.1.
        // Therefore, the activities reported here correspond to the activity as of 01/02/03"
        // I guess 01/02/03 is 1st Feb 2003, and is t0 for isotopes 210Bi and 210Pb
        // Similar info from Table 4.8 "Events from 207Bi and 152,154Eu have been weighted by their respective
        // half-lives as described in Section 3.4.1, and therefore the activities reported here correspond to their activities as of 01/02/03."
        //
        /*    else if ( histName.CompareTo("bi210_swire") == 0 ) {weight = bi210Weight * pow(0.5,eventTime/433123.);}
        else if ( (histName.CompareTo("bi210_swire") == 0) || (histName.CompareTo("bi210_sscin") == 0) || (histName.CompareTo("bi210_sfoil") == 0)) {weight = bi210Weight * pow(0.5,eventTime/433123.);}
        else if ((histName.CompareTo("co60_cuTower")==0) || (histName.CompareTo("co60_steel")==0) || (histName.CompareTo("co60_muMetal")==0) || (histName.CompareTo("co60_cuPetals")==0) )weight = pow(0.5,(eventTime/166344192.));
        else if ( (histName.CompareTo("eu152_sscin")==0) || (histName.CompareTo("eu152_int")==0)) {weight = pow(0.5,(eventTime/4.265e8));}
        else if ((histName.CompareTo("eu154_int")==0)) {weight = pow(0.5,(eventTime/2.71169e8));}
        */


        // count number of events before cuts
        // count the initial number of events
        int cc = 0;
        cut_counter[cc ++] ++;

        // 1: Accepted Run
        // (1) Accepted Run

        // is that the same as Phase cut?
        

        ///////////////////////////////////////////////////////////////////////
        // TODO: probably do NOT want to fill this histo with "weight" variable
        // do not expect to see these used in the fit
        // TODO: this is filled before all cuts, including phase cut!
        // in Summers code it was filled at the end!
        //hRun->Fill(run, weight);
        ///////////////////////////////////////////////////////////////////////


        // 1869 - 3395 is P1
        // 3396 - 9186 is P2
        cut = true;
        if(thePhase == 0)
        {
            if((1869 <= run) && (run <= 3395)) cut = false;
        }
        else if(thePhase == 1)
        {
            if((1396 <= run) && (run <= 9186)) cut = false;
        }
        else
        {
            cut = true;
        }
        // accept P1 and P2 for 
        //cut = false;
        if(cut == true) continue;
        cut_counter[cc ++]++;


        ///////////////////////////////////////////////////////////////////////
        // TODO: probably do NOT want to fill this histo with "weight" variable
        // do not expect to see these used in the fit
        //hNElectrons->Fill(nElectrons, weight);
        //moved to fill after cuts
        ///////////////////////////////////////////////////////////////////////

    /*

        // 2: Two electrons from 150Nd foil
        // (not 1, 2) - this is number not 1 -> 2 in position paper
        // (2) Two Electrons from 150Nd foil
        // 1. Exactly two electrons with reconstructed vertices inside of the 150Nd source foil boundaries.
        if(nElectrons != 2) continue;
        cut_counter[cc ++]++;

    */

        // Set the high and low energy index
        int highE_index, lowE_index;
        if(electronEnergy[0] > electronEnergy[1])
        {
            highE_index = 0;
            lowE_index = 1;
        }
        else
        {
            highE_index = 1;
            lowE_index = 0;
        }

    /*
     *
     * TODO: some of these were making cuts that removed events, but they should not have!

        ///////////////////////////////////////////////////////
        // here is a summary of all the cuts from references //
        ///////////////////////////////////////////////////////
        
        // *******************************
        // *** Position Paper, Page 61 ***
        // *******************************
        //
        // 1. Exactly two electrons with reconstructed vertices inside of the 150Nd source foil boundaries.
        // 2. The separation between each individually reconstructed vertex must be <= 4 cm radially
        // and <= 8 cm longitudinally to ensure that the electrons originate from a common vertex.
        // 3. Electron tracks with vertices located in hot spot regions are not accepted in order to
        // improve the signal to background ratio in this channel.
        // 4. Each of the electron tracks must be at least 30 cm long, have a hit in the first Gg cell layer,
        // and be associated to unique scintillator blocks.
        // 5. Tracks extrapolated to the petal blocks in the calorimeter layer closest to the source foil
        // are not accepted as their energy calibration is less precise than the other block types.
        // 6. The impact of the tracks must be on the front face of the scintillator block for proper light
        // collection.
        // 7. The energy of each electron is required to be at least 0.3 MeV, as single electron energies
        // below this threshold are not well described by the simulation.
        // 8. No alpha candidates are allowed within 15 cm of the electron vertices in order to reduce the
        // event rate coming from 214Bi decays.
        // 9. To further improve this rejection of 214Bi, no more than one prompt Geiger hit which is
        // not associated to a track is allowed within 15 cm of the track vertices. If the tracks on the
        // same side of the foil, no unassociated prompt Geiger hits are allowed on the other side of
        // the foil.
        // This allows for some level of Geiger cell noise which is not simulated, while rejecting a large
        // number of 214Bi decays.
        // 10. No unassociated Geiger hits are allowed within 15 cm of the calorimeter hits as this is
        // indicative of a reflected electron, which therefore has not deposited its full energy in the
        // calorimeter.
        // 11. Additional, isolated calorimeter hits which do not have a track associated to them are
        // allowed if the single hit energies are less than 0.2 MeV (indicative of noise which is not
        // simulated).
        // 12. The internal probability between the two electrons must be at least 1% and the external
        // probability must be no more than 1%

        // *****************************************
        // *** Position Paper, Table 6.1 Page 63 ***
        // *****************************************
        //
        // (1) Accepted Run
        // (2) Two electrons from 150Nd foil
        // (3) No other tracks in event
        // (4) Tracks intersect foil
        // (5) No alpha candidates near vertex
        // (6) Ee >= 0.2 MeV
        // (7) <= 20 total free scintillators and <= 10 gamma clusters (data management)
        // (8) All PMTs have good status
        // (9) Isolated electron calorimeter hits
        // (10) Electron hits in different scintillator blocks
        // (11) Impacts on block faces
        // (12) Electron LD flags <= 2
        // (13) Electrons have LD corrections
        // (14) Electron HS flags = 0
        // (15) Gamma hits have LD flags <= 2
        // (16) Gamma hits have LD corrections
        // (17) No gamma hits with HS = 1
        // (18) Tracks start in Gg layer 0 or 1
        // (19) Tracks miss <= 1 gg layer
        // (20) Track lengths each > 20 cm
        // (21) <= 5 NA prompt Gg hits
        // (22) Electron block types != L2/L3
        // (23) No gamma with E_gamma > 200 keV
        // (24) Tracks have negative curvature
        // (25) Track lengths > 30 cm
        // (26) Tracks not from hot spots
        // (27) No NA prompt Gg hits within 15 cm of scint hits
        // (28) <= 1 NA prompt Gg hit within 15 cm of vertex
        // (29) No NA prompt Gg hits on opposite side of foil if tracks are on same side
        // (30) Both tracks have hit in Gg L0
        // (31) Pint > 0.01 and Pext < 0.01
        // (32) Delta R <= 4 cm and delta Z <= 8 cm
        // (33) Each E_e >= 300 keV

        // *********************
        // *** Summer Thesis ***
        // *********************
        //
        // 1: Accepted Run
        // 2: Two electrons from 150Nd foil
        // 3: No other tracks in event
        // 4: No alpha candidates near vertex
        // 5: E_e >= 0.2 MeV
        // 6: <= 20 total free scintillators and <= 10 gamma clusters
        // 7: All PMTs have good status
        // 8: Isolated electron calorimeter hits
        // 9: Electron hits in different scintillator blocks
        // 10: Impacts on block faces
        // 11: All PMTs have good laser survey and counting rates identified by LD and HS flags
        // 12: Tracks start in Gg layer 0 or 1
        // 13: Tracks miss <= 1 Gg layer before impact on calorimeter
        // 14: Each L_e > 20 cm
        // 15: <= 5 unassociated prompt Gg hits
        // 16: Electron block types != L2/L3
        // 17: No gamma rays in event
        // 18: Tracks have negative curvature
        // 19: Each L_e ≥ 30 cm
        // 20: Tracks not from hot spots
        // 21: <= 1 unassociated, prompt Gg hit within 15 cm in XY coordinate of event vertex
        // 22: No unassociated, prompt Gg hits on opposite side of foil if tracks are on same side
        // 23: Both tracks have hit in Gg L0
        // 24: P_int >= 0.01 and P_ext <= 0.01
        // 25: Delta R <= 4 cm and delta Z <= 8 cm
        // 26: Each E_e >= 300 keV

        // Here is a map of how all the above points connect
        //
        // 1: (1) Accepted Run
        // 2: (2) Two electrons
        // 2: (2) from Nd150 foil
        //
        // 25: (32) 2. Delta R <= 4cm and Delta Z <= 8 cm
        //
        //


        // Numbers followed by colon : are from Summers Thesis
        // Numbers in (parenthesys) are from Table 6.1 of Page 63 from Position Paper
        // http://nile.hep.utexas.edu/DocDB/ut-nemo/docs/0037/003716/005/position_paper.pdf
        // Numbers followed by a period xx. are from List on Page 61 of Position paper

        // 3: No other tracks in event
        // (3) No other tracks in event
        // (4) Tracks intersect foil
        cut_counter[cc ++]++;
        
        // 4: No alpha candidates near vertex
        // (5) No alpha candidates near vertex
        // TODO: what does near mean?
        cut_counter[cc ++]++;
        
        // 5: E_e >= 0.2 MeV
        // (6) E_e >= 0.2 MeV
        // (23) No gamma with E_gamma > 200 keV
        if( (electronEnergy[0] < 0.2) || (electronEnergy[1] < 0.2) ) continue;
        cut_counter[cc ++]++;

        // 6: <= 20 total free scintillators and <= 10 gamma clusters
        // (7) <= 20 total free scintillators and <= 10 gamma clusters (data management)
        // NOTE: nothing to be done here, this cut is done in data management
        // by file /home/blotsd/NEMO3/release_1.0.0/SDBanalysis/src/Nd150_2e.cpp
        // by file /home/blotsd/NEMO3/release_1.0.0/SDBanalysis/src/Nd150_2eNg.cpp
        if(nGammaClusters > 10) continue; 
        cut_counter[cc ++]++;

        // 7: All PMTs have good status
        // (8) All PMTs have good status
        // TODO : not a clue how to implement this
        cut_counter[cc ++]++;

        // 8: Isolated electron calorimeter hits
        // (9) Isolated electron calorimeter hits
        // TODO: again not a clue
        cut_counter[cc ++]++;

        // 9: Electron hits in different scintillator blocks
        // (10) Electron hits in different scintillator blocks
        // [3]
        // NOTE: This is the same as cut number 27 and [3] 
        // 27: Each electron track must be associated to a unique scinitllator block
        // [3] this index appears multiple times
        if(electronPMT[0] == electronPMT[1]) continue;
        cut_counter[cc ++]++;

        // 10: Impacts on block faces
        // (11) Impacts on block faces
        // TODO: again not a clue
        cut_counter[cc ++]++;
        
        // 11: All PMT have good laser survey and counting rates identified by LD and HS flags
        // (12) Electron LD flags <= 2
        // (13) Electrons have LD corrections
        // (14) Electron HS flags = 0
        // (15) Gamma hits have LD flags <= 2
        // (16) Gamma hits have LD corrections
        // (17) No gamma hits with HS = 1
        // TODO: not sure about this, this might only be 1 part of it
        if((electronLDFlag[0] > 0) || (electronLDFlag[1] > 0)) continue;
        cut_counter[cc ++]++;
        
        // 12: Tracks start in Gg layer 0 or 1
        // (18) Tracks start in Gg later 0 or 1
        if((firstGgHitLayer[0] > 1) || (firstGgHitLayer[1] > 1)) continue;
        cut_counter[cc ++]++;
        
        // 13: Tracks miss <= 1 Gg layer before impact on calorimeter
        // (19) Tracks miss <= 1 gg layer
        // TODO: not a clue
        cut_counter[cc ++]++;
        
        // 14: Each L_e > 20 cm
        // (20) Track lengths each > 20 cm
        if((eTrackLength[0] <= 20. ) || (eTrackLength[1] <= 20.)) continue;
        cut_counter[cc ++]++;
        
        // 15: <= 5 unassociated prompt Gg hits
        // (21) <= 5 NA prompt Gg hits
        if(NAPromptGgHits > 5) continue;
        cut_counter[cc ++]++;

    */  


        // Fill histograms with relevant variables which we want to plot pre-cuts

        ///////////////////////////////////////////////////////////////////////
        /*
        Double_t NAPromptGgHitsDist2VertexMin = -1.;
        for(int i{0}; i < NAPromptGgHits; ++ i)
        {
            if(NAPromptGgHitsDist2VertexMin == -1.)
            {
                NAPromptGgHitsDist2VertexMin = NAPromptGgHitsDist2Vertex[i];
            }
            else if(NAPromptGgHitsDist2Vertex[i] < NAPromptGgHitsDist2VertexMin)
            {
                NAPromptGgHitsDist2VertexMin = NAPromptGgHitsDist2Vertex[i];
            }

            // TODO: changed to fill all values
            hNAPromptGgHitsDist2VertexMin->Fill(NAPromptGgHitsDist2Vertex[i], weight);
        }
        if(NAPromptGgHitsDist2VertexMin != -1.)
        {
            //hNAPromptGgHitsDist2VertexMin->Fill(NAPromptGgHitsDist2VertexMin, weight);
        }
        */
        // TODO: I disabled this because I couldn't be bothered to move it
        // to the end of the code to fill after all cuts with the rest of the
        // variables. Check Summers code to see how this was filled. May need
        // to move the Fill() function call with the for loop (most likely)
        // or just the Fill() by itself? (Unlikely)
        ///////////////////////////////////////////////////////////////////////



    /* WORKING HERE */

        // file: /unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2eNg_29Sep2015/CutFlow.txt
        // number of events: 690220
        // last line in file, line 24
        // matches number 15 in table from summers thesis
        // this is the first cut which should have an effect...TODO check

        // 16: Elecron block types != L2 / L3
        // (22) Electron block types != L2 / L3
        // 5. Tracks extrapolated to the petal blocks in the calorimeter layer closest to the source foil
        // are not accepted as their energy calibration is less precise than the other block types.
        if( (electronBlockType[0] == 4) || (electronBlockType[0] == 5) ) continue; // not in petals
        if( (electronBlockType[1] == 4) || (electronBlockType[1] == 5) ) continue; // not in petals
        cut_counter[cc ++] ++;
        
        // 17: No gamma rays in event
        if(nGammaClusters > 0) continue;
        // Think this is the correct way to cut any gamma rays
        // Summer removed this for Position paper
        // replacing it with no gamma rays with energy > 200 keV
        // check this TODO
        cut_counter[cc ++] ++;
        // Check that this does not cut any events - it shouldn't
        // TODO
        // Justin says there may be gamma rays in addition to 2e
        // but then what does 2eNg contain?
        // In her position paper analysis, she doesn’t require the no-γ cut, so I think we’re OK to not use it.
        // Instead, see table 6.1 on page 63 of her position paper, she just gets rid of gamma with E > 200 keV.

        // 18: Tracks have negative curvature
        // (24) Track have negative curvature
        // This appears to have been removed from Summers Position Paper
        // TODO CHECK THIS
        // how many events does this cut?
        if((trackSign[0] >= 0) || (trackSign[1] >= 0 )) continue;
        cut_counter[cc ++] ++;
        
        // 19: Each L_e >= 30 cm
        // (25) Track lengths > 30 cm
        // 4. - position paper number/index, this index appears multiple times
        // 4. Each of the electron tracks must be at least 30 cm long, [have a hit in the first Gg cell layer,
        // and be associated to unique scintillator blocks.]
        if((eTrackLength[0] < 30. ) || (eTrackLength[1] < 30.)) continue;
        cut_counter[cc ++] ++;


        ///////////////////////////////////////////////////////////////////////
        /*
        hElectronFirstGgMax->Fill(firstGgHitLayer[highE_index], weight);
        hElectronLastGgMax->Fill(lastGgHitLayer[highE_index], weight);
        hElectronFirstGgMin->Fill(firstGgHitLayer[lowE_index], weight);
        hElectronLastGgMin->Fill(lastGgHitLayer[lowE_index], weight);
        */
        // moved to fill with others after cuts
        ///////////////////////////////////////////////////////////////////////

        // 4. [Each of the electron tracks must be at least 30 cm long,] have a hit in the first Gg cell layer,
        // and be associated to unique scintillator blocks.
        if((firstGgHitLayer[0] > 0) || (firstGgHitLayer[1] > 0)) continue;
        cut_counter[cc ++] ++;
        // TODO unique scintillator blocks
        // from below
        // 23: Both tracks have hit in Gg L0
        // (not 4, 18) - position paper number/index, this index appears multiple times
        // (30) Both tracks have hit in Gg L0
        if(electronPMT[0] == electronPMT[1]) continue;
        cut_counter[cc ++] ++;
        // 27: Each electron track must be associated to a unique scinitllator block
        // (3) - position paper number/index, this index appears multiple times

        
        // 20: Tracks not from hot spots
        // (26) Tracks not from hot spots
        // 3. Electron tracks with vertices located in hot spot regions are not accepted in order to
        // improve the signal to background ratio in this channel.
        if(vertexInHotSpot[0] || vertexInHotSpot[1]) continue;
        cut_counter[cc ++] ++;
        
        // 21: <= 1 unassociated, prompt Gg hit within 15 cm in XY coordinate of event vertex
        // (28) <= 1 NA prompt Gg hit within 15 cm of vertex
        // 9. [To further improve this rejection of 214Bi, no more than one prompt Geiger hit which is
        // not associated to a track is allowed within 15 cm of the track vertices.] If the tracks on the
        // same side of the foil, no unassociated prompt Geiger hits are allowed on the other side of
        // the foil.
        // This allows for some level of Geiger cell noise which is not simulated, while rejecting a large
        // number of 214Bi decays.
        // NAPromptGgHitsDist2Vertex is distance between the Geiger hits and the average
        // electron vertex in the X,Y plane
        // The average electron vertex is obtained by adding the 2 electron vertices and
        // dividing by 2
        int ssum = 0;
        for(int i{0}; i < NAPromptGgHits; ++ i)
        {
            // assuming units of cm TODO
            if(NAPromptGgHitsDist2Vertex[i] < 15.)
            {
                ++ ssum;
            }
        }
        if(ssum >= 1) continue;
        cut_counter[cc ++] ++;
        // 22: No unassociated, prompt Gg hits on opposite side of foil if tracks are on same side
        // 9. To further improve this rejection of 214Bi, no more than one prompt Geiger hit which is
        // not associated to a track is allowed within 15 cm of the track vertices. [If the tracks on the
        // same side of the foil, no unassociated prompt Geiger hits are allowed on the other side of
        // the foil.]
        // This allows for some level of Geiger cell noise which is not simulated, while rejecting a large
        // number of 214Bi decays.
        // (29) No NA prompt Gg hits on opposite side of foil if tracks are on same side
        // TODO: not sure if this is right
        // this is clearly not right
        ssum = 0;
        if(electronSide[0] == electronSide[1])
        {
            // tracks on same side
            //if((electronFirstGG[0] > 1) && (electronFirstGG[1] > 1)) cut = true;
            for(int i{0}; i < NAPromptGgHits; ++ i)
            {
                if(electronSide[0] != NAPromptGgHitsSide[i])
                {
                    ++ ssum;
                }
            }
        }
        else if(electronSide[0] != electronSide[1])
        {
            // tracks on opposite side
            // do nothing in this case
        }
        if(ssum != 0) continue;
        cut_counter[cc ++] ++;

        // 10. No unassociated Geiger hits are allowed within 15 cm of the calorimeter hits as this is
        // indicative of a reflected electron, which therefore has not deposited its full energy in the
        // calorimeter.
        // (27) No NA prompt Gg hits within 15 cm of scintillator hits
        // This does not appear to be included in the Table from Summers thesis
        // TODO
        ssum = 0;
        for(int i{0}; i < NAPromptGgHits; ++ i)
        {
            if(NAPromptGgHitsDist2Calo[i] < 15.)
            {
                ++ ssum;
            }
        }
        if(ssum != 0) continue;
        cut_counter[cc ++] ++;
        

        
        
        // 23: Both tracks have hit in Gg L0
        // TODO: copied from above, but this looks like something which should
        // go here - might need to check against layer 1 and 11 ?
        //if((electronFirstGG[0] > 1) || (electronFirstGG[1] > 1)) cut = true; TODO: changed from this?
        // (not 4, 18) - position paper number/index, this index appears multiple times
        // (30) Both tracks have hit in Gg L0
        //if((firstGgHitLayer[0] > 0) || (firstGgHitLayer[1] > 0)) continue;
        //cut_counter[cc ++]++;
        // moved above

        // more position paper stuff TODO
        // 6. The impact of the tracks must be on the front face of the scintillator block for proper light
        // collection.
        // 8. No alpha candidates are allowed within 15 cm of the electron vertices in order to reduce the
        // event rate coming from 214Bi decays.
        // 11. Additional, isolated calorimeter hits which do not have a track associated to them are
        // allowed if the single hit energies are less than 0.2 MeV (indicative of noise which is not
        // simulated).


        ///////////////////////////////////////////////////////////////////////
        // TODO: why are these being filled with different weight
        // is it because these are not used in the fit?
        /*
        hInternalPullee->Fill(eeInternalPull, 1.0);
        hInternalProbee->Fill(eeInternalProb, 1.0);
        hExternalPullee->Fill(eeExternalPull, 1.0);
        hExternalProbee->Fill(eeExternalProb, 1.0);
        */
        // moved to fill after all cuts
        ///////////////////////////////////////////////////////////////////////

        // 24: P_{int} >= 0.01 and P_{ext} <= 0.01
        // (31) P_int > 0.01 and P_ext < 0.01
        // 12. The internal probability between the two electrons must be at least 1% and the external
        // probability must be no more than 1%
        if(eeInternalProb < 0.01) continue;
        if(eeExternalProb > 0.01) continue;
        cut_counter[cc ++] ++;
        

        double rPhi0 = (vertexR[0] * vertexSec[0] * TMath::TwoPi() / 20.);
        double rPhi1 = (vertexR[1] * vertexSec[1] * TMath::TwoPi() / 20.);
        ///////////////////////////////////////////////////////////////////////
        /*
        hVertexDZ->Fill(vertexZ[highE_index] - vertexZ[lowE_index], weight);
        hVertexDR->Fill(vertexR[highE_index] - vertexR[lowE_index], weight);
        hVertexDRPhi->Fill(rPhi0 - rPhi1, weight);
        */
        // moved to fill after all cuts
        //TODO DSec
        ///////////////////////////////////////////////////////////////////////

        // 25: Delta R <= 4cm and Delta Z <= 8 cm
        // (32) Delta R <= 4 cm and Delta Z <= 8 cm
        // 2. The separation between each individually reconstructed vertex must be ≤ 4 cm radially
        // and ≤ 8 cm longitudinally to ensure that the electrons originate from a common vertex.
        if(fabs(vertexR[0] - vertexR[1]) > 4.) continue;
        if(fabs(vertexZ[0] - vertexZ[1]) > 8.) continue;
        cut_counter[cc ++] ++;
        
        // 26: Each E_e >= 300 keV
        // (33) Each E_e >= 300 keV
        // 7. The energy of each electron is required to be at least 0.3 MeV, as single electron energies
        // below this threshold are not well described by the simulation.
        if((electronEnergy[0] < 0.3) || (electronEnergy[1] < 0.3)) continue;
        cut_counter[cc ++] ++;

        // 27: Each electron track must be associated to a unique scinitllator block
        // (3) - position paper number/index, this index appears multiple times
        //if(electronPMT[0] == electronPMT[1]) continue;
        //cut_counter[cc ++]
        // moved above

        //std::cout << "Total number of cuts recorded in cut_counter: " << cc << std::endl;
/**/

        //for(int i = 0; i < nUnAssocPromptGg; i++)
        /*
        if(thePath.CompareTo("nd150/") == 0)
        {
            if(NAPromptGgHits < 1)
            {
                std::cout << "THERE ARE NO PROMPT GG HITS FOR THIS SAMPLE" << std::endl;
            }
        }
        */
        // for Nd150, there do not appear to be any unassociated prompt geiger
        // hits
        for(int i = 0; i < NAPromptGgHits; i++)
        {
/*
            if(highE_index == 0)
            {
                hVertexMinDistPromptGgMax->Fill(dist2VertexPromptGg0[i]);
                hVertexMinDistPromptGgMin->Fill(dist2VertexPromptGg1[i]);
            }
            else
            {
                hVertexMinDistPromptGgMax->Fill(dist2VertexPromptGg1[i]);
                hVertexMinDistPromptGgMin->Fill(dist2VertexPromptGg0[i]);
            }
*/

            /*
            if(thePath.CompareTo("nd150/") == 0)
            {
                std::cout << "fill, value = " << NAPromptGgHitsDist2Vertex[i] << std::endl;
            }
            hVertexMinDistPromptGgMax->Fill(NAPromptGgHitsDist2Vertex[i]);
            */

        }
/*


*/
        //hNAPromptGgHits->Fill(NAPromptGgHits);
        // moved below
/*
        double totalLowE = 0.;
        for(int i = 0; i < nLowEClusters; i++)
        {
            totalLowE += lowEClusterE[i];
            hLowEGammaEnergy->Fill(lowEClusterE[i], weight);
            hLowEMinDistPromptGg->Fill(lowEMinDistPromptGg[i], weight);
        }
        hNLowEGammas->Fill(nLowEClusters, weight);
        hSummedLowEGammaE->Fill(totalLowE, weight);

        //if ( nLowEClusters > 0 )continue;
*/

        // moved back here from above
        hRun->Fill(run, weight);
        hNElectrons->Fill(nElectrons, weight);


        // moved from above
        hElectronFirstGgMax->Fill(firstGgHitLayer[highE_index], weight);
        hElectronLastGgMax->Fill(lastGgHitLayer[highE_index], weight);
        hElectronFirstGgMin->Fill(firstGgHitLayer[lowE_index], weight);
        hElectronLastGgMin->Fill(lastGgHitLayer[lowE_index], weight);

        hInternalPullee->Fill(eeInternalPull, 1.0);
        hInternalProbee->Fill(eeInternalProb, 1.0);
        hExternalPullee->Fill(eeExternalPull, 1.0);
        hExternalProbee->Fill(eeExternalProb, 1.0);

        hVertexDZ->Fill(vertexZ[highE_index] - vertexZ[lowE_index], weight);
        hVertexDR->Fill(vertexR[highE_index] - vertexR[lowE_index], weight);
        hVertexDRPhi->Fill(rPhi0 - rPhi1, weight);

        hNAPromptGgHits->Fill(NAPromptGgHits);
        // these were moved from above


        hTotalE->Fill(electronEnergy[0] + electronEnergy[1], weight);
        hCosee->Fill(cosee, weight);
        hCoseeWeighted->Fill(cosee, weight * cosee_weight);

        hEeMax->Fill(electronEnergy[highE_index], weight);
        hElectronLengthMax->Fill(eTrackLength[highE_index], weight);
        hVertexSectorMax->Fill(vertexSec[highE_index], weight);
        hVertexZMax->Fill(vertexZ[highE_index], weight);
        hVertexRMax->Fill(vertexR[highE_index], weight);
        hVertexZSecMax->Fill(vertexSec[highE_index], vertexZ[highE_index], weight);
        hElectronLDCorrMax->Fill(electronLDCorr[highE_index], weight);

        hEeMin->Fill(electronEnergy[lowE_index], weight);
        hElectronLengthMin->Fill(eTrackLength[lowE_index], weight);
        hVertexSectorMin->Fill(vertexSec[lowE_index], weight);
        hVertexZMin->Fill(vertexZ[lowE_index], weight);
        hVertexRMin->Fill(vertexR[lowE_index], weight);
        hVertexZSecMin->Fill(vertexSec[lowE_index], vertexZ[lowE_index], weight);
        hElectronLDCorrMin->Fill(electronLDCorr[lowE_index], weight);

        hEeMaxVEeMin->Fill(electronEnergy[highE_index], electronEnergy[lowE_index], weight);


        // filling of hNAPromptGgHitsDist2VertexMin moved before cuts

/*
        for(int i = 0; i < nUnAssocPromptGg; i++)
        {
            if(highE_index == 0)
            {
                hVertexMinDistPromptGgMax->Fill(dist2VertexPromptGg0[i]);
                hVertexMinDistPromptGgMin->Fill(dist2VertexPromptGg1[i]);
    //new variable names:
    double NAPromptGgHitsDist2Vertex[20]; // TODO check number 20
    double NAPromptGgHitsDist2Calo[20];
    //old variable names:
    //double      dist2VertexPromptGg0[20];
    //double      dist2VertexPromptGg1[20];
            }
            else
            {
                hVertexMinDistPromptGgMax->Fill(dist2VertexPromptGg1[i]);
                hVertexMinDistPromptGgMin->Fill(dist2VertexPromptGg0[i]);
            }
        }

        hNUnassocPromptGgHits->Fill(nUnAssocPromptGg);

        double totalLowE = 0.;
        for(int i = 0; i < nLowEClusters; i++)
        {
            totalLowE += lowEClusterE[i];
            hLowEGammaEnergy->Fill(lowEClusterE[i], weight);
            hLowEMinDistPromptGg->Fill(lowEMinDistPromptGg[i], weight);
        }

        hNLowEGammas->Fill(nLowEClusters, weight);
        hSummedLowEGammaE->Fill(totalLowE, weight);
*/
        //hNLowEGammas->Fill();
        //hLowEGammaEneegy->Fill();
        //hSummedLowEGammaE->Fill();
        //hLowEMinDistPromptGg->Fill();

    } //~event_i
        
    std::cout << std::endl;
    
    aFile->Close();
    
    std::cout << "Here are the cut counts" << std::endl;
    for(int i = 0; i < nCuts; i++)
    {
        std::cout << "Passed cut " << i << ": " << cut_counter[i] << " events." << std::endl;
    }

    TFile *myFile = TFile::Open("Nd150_2e_P" + Phase + ".root", "UPDATE");
    TDirectory *dir = gDirectory;
    dir->cd("singleHistos/unscaled");

    //dir->pwd();
    //dir->print();
    //dir->Print();

    std::cout << "writing histograms into \"singleHistos/unscaled\'" << std::endl;
    std::cout << "writing: singleHistos/unscaled/" << hTotalE->GetName() << std::endl;

    hRun->Write();
    hNElectrons->Write();
    hTotalE->Write();
    hCosee->Write();
    hInternalPullee->Write();
    hInternalProbee->Write();
    hExternalPullee->Write();
    hExternalProbee->Write();
    hCosee->Write();
    hCoseeWeighted->Write();
    hVertexDZ->Write();
    hVertexDR->Write();
    hVertexDRPhi->Write();
    hNAPromptGgHits->Write();

    hEeMax->Write();
    hElectronLengthMax->Write();
    hVertexZMax->Write();
    hVertexSectorMax->Write();
    hVertexRMax->Write();
    hElectronFirstGgMax->Write();
    hElectronLastGgMax->Write();
    hVertexMinDistPromptGgMax->Write();
    hElectronLDCorrMax->Write();
    hVertexZSecMax->Write(); //

    hEeMin->Write();
    hElectronLengthMin->Write();
    hVertexZMin->Write();
    hVertexSectorMin->Write();
    hVertexRMin->Write();
    hElectronFirstGgMin->Write();
    hElectronLastGgMin->Write();
    hVertexMinDistPromptGgMin->Write();
    hElectronLDCorrMin->Write();
    hVertexZSecMin->Write(); //

    hNLowEGammas->Write();
    hLowEGammaEnergy->Write();
    hSummedLowEGammaE->Write();
    hLowEMinDistPromptGg->Write();

    hEeMaxVEeMin->Write();


    for(int i = 0; i < nCuts; i++ )
    {
        // NOTE: bug, should be i + 1, bin 0 is underflow, bin -1 does not exist
        //hNAfterCuts->SetBinContent(i - 1, cut_counter[i]);
        hNAfterCuts->SetBinContent(i + 1, cut_counter[i]);
    }

    hNAfterCuts->Write();
    // hEffAfterCuts->Write();
    
    
    hNAPromptGgHitsDist2VertexMin->Write();

    myFile->Close();

}


/*
static const int nInternalBkgs = 10;
TString InternalBkgFiles[nInternalBkgs] = {"bi212_int_rot","ac228_int_rot","bi207_int_rot","bi214_int_rot","eu152_int_rot","eu154_int_rot","k40_int_rot","pa234m_int_rot","pb214_int_rot","tl208_int_rot"};
TString InternalBkgNames[nInternalBkgs] = {"^{212}Bi int rot","^{228}Ac int rot","^{207}Bi int rot","^{214}Bi int rot","^{152}Eu int rot","^{154}Eu int rot","^{40}K int rot","^{234m}Pa int rot","^{214}Pb int rot","^{208}Tl int rot"};
double InternalBkgActivity[2][nInternalBkgs], InternalBkgEfficiency[2][nInternalBkgs];
double  InternalBkgNGenMC[nInternalBkgs];
Color_t InternalBkgColor[2] = {kPink+1, kPink}; //for grouping
Color_t InternalBkgColors[nInternalBkgs] = {kCyan+2, kCyan, kTeal+2, kRed, kRed+2, kGreen, kGreen+2, kYellow, kOrange+7,kViolet+1}; 
*/
void scale(TFile* myFile,                       // INPUT: unscaled histograms are loaded from this root file
           const std::string& activityfile,     // INPUT: activity values read from this file
           const std::string& typedir,          // INPUT: forms part of path string, eg: "nd150/", "internals/", "externals/", "neighbours/", all radon files are in "externals/"
           const int nSamples,                     // INPUT: number of objects in dataFiles/dataNames
           const TString *const sampleFiles,      // INPUT: list of samples, eg: bi214_sfoil_rot
           const TString *const sampleNames,      // INPUT: corresponding human readable names, Latex format, eg: ^{214}Bi Foil Surface
           const Color_t *const sampleColors,     // INPUT: sample color set using values in this array
           const double *const AcceptedTime,    // INPUT: accepted time for each phase is passed in here 
           double *const sampleNGenMC,            // OUTPUT: number of generated MC events are written to this location
           double *const sampleActivity,          // OUTPUT: activity values written to this location
           double *const sampleEfficiency,        // OUTPUT: efficiency values written to this location
           TObjArray** allSamples)              // OUTPUT: each histogram has colors set and is scaled (using scaling with activity), before being stored as an output here
{
    std::cout << "scale" << std::endl;

    std::ofstream oftemp("oftemp1.txt", std::ofstream::out | std::ofstream::app);

    TH1F *tmpHist; //we'll be using this later.

    //Let's get the internals/externals/data and add them to an array
    //TObjArray *allSamples[nSamples];
    for(int i = 0; i < nSamples; i++)
    {
        // create new object array
        //std::cout << sampleFiles[i] << std::endl;
        allSamples[i] = new TObjArray();

        // read in the information needed for calculating efficiencies
        // read in number of generated MC events
        std::ifstream inFile;
        inFile.open(filePath + typedir + sampleFiles[i] + "/JobSummary.txt");
        // cout << "JobSummary stuff is here: "
        // 	 << filePath << "internals/" << InternalBkgFiles[i] << "/JobSummary.txt"
        // 	 << endl;
        std::string dummy;
        inFile >> dummy >> sampleNGenMC[i];
        // cout << "Read in as dummy: " << dummy << endl;
        // cout << "Read in as InternalBkgNGenMC[i]: " << InternalBkgNGenMC[i] << endl;
        inFile.close();

        // if we want to use a pre-measured background activity (which we do for the time being) then we
        // have to scale the histogram appropriately.  So now we need to fetch the activity from the txt file
        std::ifstream inFile2;
        inFile2.open(activityfile.c_str());

        while(!inFile2.eof())
        {
            // phase 1 and phase 2 activity
            double activityPhase1, activityPhase2; 
            // name of MC sample
            TString mcname;
            
            std::stringstream ss;
            std::string s;
            std::getline(inFile2, s);

            if(s.size() > 0)
            {
                ss << s;
                ss >> mcname >> activityPhase1 >> activityPhase2;
                //std::cout << mcname << ", " << activityPhase1 << ", " << activityPhase2 << std::endl;
                
                // cout << "Read in from activity file: " << mcname
                // 	   << ", " << activityPhase1
                // 	   << ", " << activityPhase2 << endl;
                // cout << "Trying to compare to " << InternalBkgFiles[i] << endl;
                
                if(mcname.CompareTo(sampleFiles[i]) == 0)
                {
                    //sampleActivity[0][i] = activityPhase1;
                    //sampleActivity[1][i] = activityPhase2;
                    sampleActivity[0 * nSamples + i] = activityPhase1;
                    sampleActivity[1 * nSamples + i] = activityPhase2;

                    if(mcname.CompareTo("mo100_99_rot_2n2b_m14") == 0)
                    {
                        std::cout << "tracking m0 bb" << std::endl;
                        std::cout << activityPhase1 << ", " << activityPhase2 << std::endl;
                        //std::cin.get();
                    }
                }

            }
        }

        inFile2.close();

        for(int j = 0; j < numHistograms; j++)
        {
            //std::cout << "histogramNames[" << j << "]=" << histogramNames[j] << " sampleFiles[" << i <<"]=" << sampleFiles[i] << std::endl;

            std::string name("singleHistos/unscaled/" + histogramNames[j] + sampleFiles[i]);
            std::cout << "Get() : " << name << " from file, Clone() : " << histogramNames[j] + sampleFiles[i] + "_fit" << std::endl;
            tmpHist = (TH1F*)myFile->Get(name.c_str())->Clone(histogramNames[j] + sampleFiles[i] + "_fit");
            //tmpHist = (TH1F*)myFile->Get("singleHistos/unscaled/" + histogramNames[j] + sampleFiles[i])->Clone(histogramNames[j] + sampleFiles[i] + "_fit");
            // cout << "The histogram I have got is "
            // 	   << "singleHistos/unscaled/"
            // 	   << histogramNames[j]
            // 	   << InternalBkgFiles[i]
            // 	   << endl;
            // cout << "At this point, the integral is " << tmpHist->Integral() << endl;

            //calculate efficiency
            // TODO: this is calculated using the 0th histogram, should perhaps
            // change to the hTotalE ?
            if(j == 0)
            {
                double NPass = tmpHist->Integral();
                // NPass is phase specific
                // sampleNGenMC[i] is for both phases
                std::cout << "NPass: " << NPass << std::endl;
                std::cout << "sampleNGenMC[i=" << i << "]=" << sampleNGenMC[i] << std::endl;
                
                // TODO: don't fully understand this
                // efficiency = (number of events which pass cuts * total time) / 
                // (number of generated events * accepted time)
                // accepted time is like total time - dead time?
                //double eff = NPass / (double)(sampleNGenMC[i] * (AcceptedTime[thePhase]/(double)TotalTime) );
                // this efficiency is phase specific
                double eff = (double)NPass / (double)sampleNGenMC[i];
                //double dead_time_scale_factor = (double)TotalTime / (double)AcceptedTime[thePhase];
                //double efficiency = eff * dead_time_scale_factor;
                double phase_scale_factor = (double)TotalTime / (double)AcceptedTime[thePhase];
                // this efficiency is rescaled for to be phase independent
                double efficiency = eff * phase_scale_factor;
                std::cout << "efficiency: " << efficiency << std::endl;

                //sampleEfficiency[thePhase][i] = eff;
                sampleEfficiency[thePhase * nSamples + i] = efficiency;
            }
            tmpHist->Sumw2();
            tmpHist->SetFillColor(sampleColors[i]);
            tmpHist->SetLineColor(sampleColors[i]);
            tmpHist->SetTitle(sampleNames[i]);
            // cout << "Integral before scaling: " << tmpHist->Integral() << endl;
            // cout << "InternalBkgActivity[thePhase][i]: " << InternalBkgActivity[thePhase][i] << endl
            // 	   << "AcceptedTime[thePhase]: " << AcceptedTime[thePhase] << endl
            // 	   << "InternalBkgNGenMC[i]: " << InternalBkgNGenMC[i] << endl
            // 	   << "AcceptedTime[thePhase]: " << AcceptedTime[thePhase] << endl
            // 	   << "TotalTime: " << TotalTime << endl;
            //tmpHist->Scale(sampleActivity[thePhase][i] * AcceptedTime[thePhase] * 0.001 / (sampleNGenMC[i] * (AcceptedTime[thePhase]/(double)TotalTime)) );
            //tmpHist->Scale(sampleActivity[thePhase * nSamples + i] * 0.001 * (double)TotalTime / sampleNGenMC[i] );
            //std::cout << "scaling, index=" << thePhase + nSamples + i << std::endl;
            
            // In Position Paper
            // Section 3.4.1, Equation 3.7
            // A = N / (t * eps)
            // A: activity, N: number of events, t: time, eps: efficiency
            // Can calculate number of events backwards from activity, efficiency and time

            //Double_t scale_factor = sampleActivity[thePhase * nSamples + i] * AcceptedTime[thePhase] * additional_scaling_factor / (sampleNGenMC[i] * (AcceptedTime[thePhase]/(double)TotalTime));
            Double_t the_activity = sampleActivity[thePhase * nSamples + i];
            //std::cout << sampleFiles[i] << " -> " << sampleNames[i] << std::endl;
            //Double_t n_expected_events = the_activity * (double)TotalTime;
            //Double_t n_expected_events = the_activity * (double)AcceptedTime[thePhase];
            Double_t n_expected_events = the_activity * (double)TotalTime;
            Double_t scale_factor = n_expected_events / sampleNGenMC[i];

            //Double_t new_scale_factor = (1.0 / sampleNGenMC[i]) * (the_activity * (double)AcceptedTime[thePhase]);

            if(sampleFiles[i].CompareTo("mo100_99_rot_2n2b_m14") == 0 || sampleFiles[i].CompareTo("nd150_rot_2n2b_m4") == 0)
            {
                std::cout << sampleFiles[i] << ", " << the_activity << ", " << AcceptedTime[thePhase] << ", " << n_expected_events << ", " << sampleNGenMC[i] << ", " << scale_factor << std::endl;
                //std::cin.get();
            }
            // using the above information from Position paper, this scale_factor should be N/(efficiency*Ngen)
            // (N is the number that pass cuts, Ngen is the number generated)
            // N/Ngen = efficiency
            // therefore scale_factor = 1.
            // this cannot possibly be true
            //Double_t effi = tmpHist->Integral() / sampleNGenMC[i];
            //scale_factor /= effi;
            
            
            
            if(std::string(tmpHist->GetName()).find("hTotalE") != std::string::npos ||
               std::string(tmpHist->GetName()).find("hRun") != std::string::npos)
            {

                // TODO: working here, the activity can change between phases
                // not sure what the best way is to work out the number
                // of expected events on each phase and compare this to
                // the histogram integrals
                // also the number of events appears to be different depending
                // on whether hTotalE is used or whether hRun is used
                // does hRun contain events for both phases?
                // what are the limits on hRun? how is it filled?
                // copied from filling code:
                //hRun->Fill(run, weight);
                //hTotalE->Fill(electronEnergy[0] + electronEnergy[1], weight);

                double NPass = tmpHist->Integral();
                //std::cout << tmpHist->GetName() << "->Integral()=" << tmpHist->Integral() << std::endl;
                std::cout << tmpHist->GetName() << "->Integral()=" << NPass << std::endl;
                std::cout << "thePhase=" << thePhase << ", nSamples=" << nSamples << ", i=" << i << std::endl;
                std::cout << "sampleEfficiency[" << thePhase * nSamples + i << "]=" << sampleEfficiency[thePhase * nSamples + i] << std::endl;
                double eff = (double)NPass / (double)sampleNGenMC[i];
                //double dead_time_scale_factor = (double)TotalTime / (double)AcceptedTime[thePhase];
                double phase_scale_factor = (double)TotalTime / (double)AcceptedTime[thePhase];
                //double efficiency = eff * dead_time_scale_factor;
                // this scales the efficiency back up to compensate for the
                // events that were cut due to the phase cut
                double efficiency = eff * phase_scale_factor;
                std::cout << "efficiency calculated using " << tmpHist->GetName() << ": eff=" << eff << ", efficiency=" << efficiency << ", NPass=" << NPass << ", sampleNGenMC[" << i << "]=" << sampleNGenMC[i] << std::endl;
                std::cout << "check activity read in: " << sampleActivity[thePhase * nSamples + i] << std::endl;
                // TODO: must check all these calculations, particularly what do the times actually correspond to
                // this value should match the integral
                double my_n_expected = sampleActivity[thePhase * nSamples + i] * TotalTime * efficiency;
                std::cout << "number of expected events is: " << my_n_expected << " ~ should match NPass=" << NPass << std::endl;
                if(std::string(tmpHist->GetName()).find("hTotalE") != std::string::npos)
                {
                    std::cin.get();
                }
            }




            
            //std::cout << "*** the_activity=" << the_activity << ", n_expected_events=" << n_expected_events << ", sampleNGenMC[i]=" << sampleNGenMC[i] << ", scale_factor=" << scale_factor << ", integral()-> " << tmpHist->Integral() << std::endl;

            oftemp << sampleFiles[i] << " -> " << sampleNames[i] << std::endl;
            oftemp << "*** the_activity=" << the_activity << ", n_expected_events=" << n_expected_events << ", sampleNGenMC[i]=" << sampleNGenMC[i] << ", scale_factor=" << scale_factor << ", integral()-> " << tmpHist->Integral() << std::endl;

            //tmpHist->Scale(scale_factor); 
            //tmpHist->Scale(new_scale_factor); 
            
            // moved from stackFunction()
            // 2020-03-30: 
            //(TH1F*)allSamples[j]->At(i)->Write();
            // saving the histogram without scaling it by the activity will
            // cause the fit output (and input) to be in units of activty
            // is this actually the case, what about the factor of the detector
            // livetime which I have left in above?
            Double_t scale_factor_without_activity = (double)TotalTime / sampleNGenMC[i]; //scale_factor / the_activity;
            tmpHist->Scale(scale_factor_without_activity);
            std::cout << "writing: " << tmpHist->GetName() << std::endl;
            tmpHist->Write();

            tmpHist->Scale(scale_factor / scale_factor_without_activity);

            //std::cout << " : " << histogramNames[j] << " -> sf=" << scale_factor << ", sampleActivity[]=" << sampleActivity[thePhase * nSamples + i] << " integral()=> " << tmpHist->Integral() << std::endl;

            //std::cout << "sf=" << sampleActivity[thePhase * nSamples + i] * AcceptedTime[thePhase] * additional_scaling_factor / (sampleNGenMC[i] * (AcceptedTime[thePhase]/(double)TotalTime)) << std::endl;
            oftemp << " : " << histogramNames[j] << " -> sf=" << scale_factor << ", sampleActivity[]=" << sampleActivity[thePhase * nSamples + i] << " integral()=> " << tmpHist->Integral() << std::endl << std::endl;
            //tmpHist->Scale(sampleActivity[thePhase * nSamples + i] * AcceptedTime[thePhase] / (sampleNGenMC[i] * (AcceptedTime[thePhase]/(double)TotalTime)) );
            // cout << "Integral after scaling: " << tmpHist->Integral() << endl;

            //tmpHist->Scale(InternalBkgActivity[thePhase][i] * AcceptedTime[thePhase] * 0.001    / (InternalBkgNGenMC[i] * (AcceptedTime[thePhase]/(double)TotalTime)) );
            //tmpHist->Scale(NeighbourActivity[thePhase][i]   * AcceptedTime[thePhase] * 0.001    / (NeighbourNGenMC[i]   * (AcceptedTime[thePhase]/(double)TotalTime)) );
            //tmpHist->Scale(Nd150Activity[thePhase][i]       * AcceptedTime[thePhase] * 0.001    / (Nd150NGenMC[i]       * (AcceptedTime[thePhase]/(double)TotalTime)) );
            //tmpHist->Scale(ExternalBkgActivity[thePhase][i] * AcceptedTime[thePhase]            / (ExternalBkgNGenMC[i] * (AcceptedTime[thePhase]/(double)TotalTime)) );
            //tmpHist->Scale(ExternalBkgActivity[thePhase][i] * AcceptedTime[thePhase]            / (ExternalBkgNGenMC[i] * (AcceptedTime[thePhase]/(double)TotalTime)) );
            //tmpHist->Scale(Rn222BkgActivity[thePhase][i]    * AcceptedTime[thePhase]            / (Rn222BkgNGenMC[i]    * (AcceptedTime[thePhase]/(double)TotalTime)) );
            //tmpHist->Scale(Rn220BkgActivity[thePhase][i]    * AcceptedTime[thePhase]            / (Rn220BkgNGenMC[i]    * (AcceptedTime[thePhase]/(double)TotalTime)) );


            // NOTE: done AFTER call to Write()
            // so sample is scaled with activity value if used again from allSamples
            allSamples[i]->Add(tmpHist);

        }//all histos
    }//internals

    oftemp.close();
}

// allSamples are individual histograms, index order: histogram type, background type
// index order: array index ([]), followed by ->At() function call index for TObjArray
// hAllSamples is a single histogram, index order: histogram type
// index order: array index ([])
// hAllMc is a single histogram, index order: histogram type
// index order: array index ([])
// same as hAllSamples, but contains all the different types of background/mc/data in a single
// histogram rather than splitting between several different objects (arrays of pointers) one
// for each type of background/mc/data
void stackfunction(int i,                           // INPUT: histogram index, should possibly be moved inside function (can be done by moving calls to new for hMinorStacks inside, and calling when allMC_first == true?)
                   int nSamples,                    // INPUT: number of samples in list (number of objects in lists below)
                   THStack **hMinorStacks,          // OUTPUT: ?
                   TObjArray **allSamples,          // INPUT: all samples
                   THStack **hMajorStacks,          // OUTPUT: ?
                   const TString& sampleHumanName,  // INPUT:
                   const TString& sampleName,       // INPUT: "Nd150", "Externals", "Internals", "Radon", used to when calling Clone() to change histogram name
                   TH1F **hAllSamples,              // OUTPUT: single histogram stack composed of samples components added together by function call to TH1::Add()
                   const TString *const sampleNames,
                   const TString *const sampleFiles,
                   TH1F **hAllMC,
                   bool allMC_first,
                   bool hAllSamples_addtoexisting,  // INPUT: flag, if true objects in allSamples are added to hAllSamples without call to Clone to change name
                   Color_t majorStackColor,
                   bool hMajorStacks_add)
{

    //    stackfunction(i,
    //                  nNd150Samples,
    //                  hMinorStacks,
    //                  allNd150Samples,
    //                  "Nd150",
    //                  hAllNd150,
    //                  Nd150Names,
    //                  Nd150Files,
    //                  hAllMC,
    //                  false,
    //                  false
    //                  );


    // i = histogram index (type, eg; electron energy min/max, total energy, track length...)
    // j = background/data/mc index

    TH1F *tmpHist; //we'll be using this later.

    /*
    if(genericName.CompareTo("Nd150") == 0)
    {
        std::cout << "This is Nd150" << std::endl;
    }
    */
    
    // this block relating to hAllSamples
    // which is one of hAllExternals, hAllInternals, etc...
    // add all samples of a type into a histogram for the group hAllSamples
    for(int j = 0; j < nSamples; j++)
    {
        if(hAllSamples_addtoexisting == false)
        {
            if(j == 0)
            {
                // NOTE: allSamples is either allExternals, allInternals, etc...
                hAllSamples[i] = (TH1F*)allSamples[j]->At(i)->Clone(sampleName);
            }
            else
            {
                // NOTE: allSamples is either allExternals, allInternals, etc...
                hAllSamples[i]->Add((TH1F*)allSamples[j]->At(i));
            }
        }
        else
        {
            // NOTE: allSamples is either allExternals, allInternals, etc...
            hAllSamples[i]->Add((TH1F*)allSamples[j]->At(i));
        }

        // 2020-03-31: moved from below block to this location
        tmpHist = (TH1F*)allSamples[j]->At(i);
        if(tmpHist->Integral() > 0)
        {
            //std::cout << ExternalBkgFiles[j] << " " << tmpHist->Integral() << std::endl;
            tmpHist->SetTitle(sampleNames[j]); // TODO: want this above?
            hMinorStacks[i]->Add((TH1F*)allSamples[j]->At(i)->Clone(sampleFiles[j]));
        }
        else
        {
            #if 0
            std::cout << "MC: " << sampleNames[j] << ", HIST: " << tmpHist->GetName() << ", Integral() == 0" << std::endl;
            #endif
        }


        // below block all related to hAllMC

        // allMC goes here
        if(allMC_first == true)
        {
            if(j == 0)
            {
                // NOTE: allSamples is either allExternals, allInternals, etc...
                hAllMC[i] = (TH1F*)allSamples[j]->At(i)->Clone("Total MC");
            }
            else
            {
                // NOTE: allSamples is either allExternals, allInternals, etc...
                hAllMC[i]->Add((TH1F*)allSamples[j]->At(i));
            }
        }
        else
        {
            // NOTE: allSamples is either allExternals, allInternals, etc...
            hAllMC[i]->Add((TH1F*)allSamples[j]->At(i));
        }

   
        //(TH1F*)allSamples[j]->At(i)->Write();
        // 2020-03-30: removed
        

    } // generic


    TString events;
    //std::cout << "150Nd: Add()" << std::endl;
    //events.Form("%i",(int)hAllNd150[i]->Integral());
    //hAllNd150[i]->SetTitle("^{150}Nd 2#nu2#beta ("+events+")");
    //hMajorStacks[i]->Add((TH1F*)hAllNd150[i]);
    

    // TODO: there are still bugs in this section

    
    if(hMajorStacks_add == true)
    {
        /*
        if(i == 2)
        {
            // hTotalE
            std::cout << "hAllSamples[" << i << "]->GetName() -> " << hAllSamples[i]->GetName() << std::endl;
            if(sampleName.CompareTo("Radon") == 0)
            {
                std::cout << " Integral() -> " << hAllSamples[i]->Integral() << std::endl;
            }
        }
        */

        std::cout << "generic sample name: Add()" << std::endl;
        events.Form("%i", (int)hAllSamples[i]->Integral());
        TString title_string = sampleHumanName + " (" + events + ")";
        //tmpHist->SetFillColor(dataColors[i]);
        //tmpHist->SetLineColor(dataColors[i]);
        hAllSamples[i]->SetFillColor(majorStackColor);
        hAllSamples[i]->SetLineColor(majorStackColor);
        hAllSamples[i]->SetTitle(title_string);
        hMajorStacks[i]->Add((TH1F*)hAllSamples[i]);
        
        //std::cout << "Externals: Add()" << std::endl;
        //events.Form("%i",(int)hAllExternals[i]->Integral());
        //hAllExternals[i]->SetTitle("Externals ("+events+")");
        //hMajorStacks[i]->Add((TH1F*)hAllExternals[i]);
    }
    else
    {
        // do nothing
    
        /*
        if(i == 2)
        {
            std::cout << "hAllSamples[" << i << "]->GetName() -> " << hAllSamples[i]->GetName() << std::endl;
            if(sampleName.CompareTo("Radon") == 0)
            {
                std::cout << "need to add this integral: " << hAllSamples[i]->Integral() << std::endl;
            }
        }
        std::cout << "generic sample name: Add()" << std::endl;
        int previous_count = ((TH1F*)hMajorStacks[i]->GetHists()->Last())->Integral();
        events.Form("%i", (int)hAllSamples[i]->Integral() + previous_count);
        TString title_string = sampleHumanName + " (" + events + ")";
        hAllSamples[i]->SetTitle(title_string);
        //hMajorStacks[i]->Add((TH1F*)hAllSamples[i]);
        ((TH1F*)hMajorStacks[i]->GetHists()->Last())->Add((TH1F*)hAllSamples[i]);
        */
    }
    

}


void fitHistograms()
{

    std::cout << "fitHistograms()" << std::endl;

    TFile *myFile = TFile::Open("Nd150_2e_P"+Phase+".root", "UPDATE");
    TDirectory* dir = gDirectory; 
    std::cout << "dir->ls():" << std::endl;
    dir->ls();
    TObjArray *allDataHistograms = new TObjArray();    
    
    //Get the data histograms
    for ( int i = 0; i < numHistograms; i++ ) {
        std::string name("singleHistos/unscaled/" + histogramNames[i] + DataFile);
        //std::cout << "name=" << name << std::endl;
        std::cout << "Get() : " << name << " from file, Clone() : " << histogramNames[i] + "data" << std::endl;
        allDataHistograms->Add((TH1F*)myFile->Get(name.c_str())->Clone(histogramNames[i] + "data"));
        // 2020-04-02: changed "_data" to "data"
    }


    // TODO: TotalTime is 33859178. + 133770114.
    // is this correct? perhaps this is wrong?
    // considering dead_time_scale_factor
    // unless dead_time_scale_factor is actually a re-scaling for the entire
    // detector P1 + P2 time
    AcceptedTime[0] = 33859178.;
    AcceptedTime[1] = 133770114.;
    AcceptedTime[2] = TotalTime;
                              
    std::ofstream oftemp("oftemp1.txt", std::ofstream::out | std::ofstream::trunc);
    oftemp.close();

/*
void scale(TFile* myFile,                       // INPUT: unscaled histograms are loaded from this root file
           const std::string& activityfile,     // INPUT: activity values read from this file
           const std::string& typedir,          // INPUT: forms part of path string, eg: "nd150/", "internals/", "externals/", "neighbours/", all radon files are in "externals/"
           const int nData,                     // INPUT: number of objects in dataFiles/dataNames
           const TString *const dataFiles,      // INPUT: list of samples, eg: bi214_sfoil_rot
           const TString *const dataNames,      // INPUT: corresponding human readable names, Latex format, eg: ^{214}Bi Foil Surface
           const Color_t *const dataColor,      // INPUT: (unused)
           const Color_t *const dataColors,     // INPUT: data color set using values in this array
           const double *const AcceptedTime,    // INPUT: accepted time for each phase is passed in here 
           double *const dataNGenMC,            // OUTPUT: number of generated MC events are written to this location
           double *const dataActivity,          // OUTPUT: activity values written to this location
           double *const dataEfficiency,        // OUTPUT: efficiency values written to this location
           TObjArray** allSamples)              // OUTPUT: each histogram has colors set and is scaled (using scaling with activity), before being stored as an output here
*/
    // documenation for scale() function:
    //
    // myFile is a TFile containing all raw histograms (unscaled)
    // "../include/activities.txt" -> activityfile is a path and filename for activity file which contains all sample activties (pre-fit)
    // "internals/" -> typedir forms part of path string to MC files
    // nInternalBkgs -> nData contains number of MC histograms to process
    // InternalBkgFiles -> dataFiles contains list of histogram names, eg: bi214_int_rot
    // InternalBkgNames -> dataNames contains corresponding human readable names, eg: ^{214}Bi Internal
    // InternalBkgColor -> dataColor is unused
    // InternalBkgColors -> dataColors is used to set the line and fill color for each MC sample
    // AcceptedTune Phase 1 and Phase 2 accepted time
    // InternalBkgNGenMC -> dataNGenMC is an output containing number of generated MC
    // InternalBkgActivity -> dataActivity is an output containing the activity for Phase 1 and Phase 2
    // InternalBkgEfficiency -> dataEfficiency is an output containing the efficiency for Phase 1 and Phase 2
    // allInternals -> allSamples is an output, format TObjArray, contains each MC histogram (as a clone) in a TObjArray

    // after scale function, these contain
    // TH1F objects
    // the named object below is an array of TObjArray pointers
    // the TObjArray is allocated with new in scale
    // there is one for each background type
    // a single TH1F* object is added to the TObjArray
    // so it is an array of length 1
    TObjArray *allExternals[nExternalBkgs];
    TObjArray *allRn220Bkgs[nRn220Bkgs];
    TObjArray *allRn222Bkgs[nRn222Bkgs];
    TObjArray *allInternals[nInternalBkgs];
    TObjArray *allNeighbours[nNeighbours];
    TObjArray *allNd150Samples[nNd150Samples];

#if INTERNALS_ON
    std::cout << "scale: internals" << std::endl;
    scale(myFile,
          "../include/activities.txt",
          "internals/",
          nInternalBkgs,
          InternalBkgFiles,
          InternalBkgNames,
          InternalBkgColors,
          AcceptedTime,
          InternalBkgNGenMC,
          &InternalBkgActivity[0][0],
          &InternalBkgEfficiency[0][0],
          allInternals);
#endif
#if NEIGHBOURS_ON
    std::cout << "scale: neighbours" << std::endl;
    scale(myFile,
          "../include/activities.txt",
          "neighbours/",
          nNeighbours,
          NeighbourFiles,
          NeighbourNames,
          NeighbourColors,
          AcceptedTime,
          NeighbourNGenMC,
          &NeighbourActivity[0][0],
          &NeighbourEfficiency[0][0],
          allNeighbours);
    
#endif
#if ND_ON
    std::cout << "scale: 150Nd" << std::endl;
    scale(myFile,
          "../include/activities.txt",
          "nd150/",
          nNd150Samples,
          Nd150Files,
          Nd150Names,
          Nd150Colors,
          AcceptedTime,
          Nd150NGenMC,
          &Nd150Activity[0][0],
          &Nd150Efficiency[0][0],
          allNd150Samples);
#endif
#if EXTERNALS_ON
    std::cout << "scale: externals" << std::endl;
    scale(myFile,
          "../include/activities.txt",
          "externals/",
          nExternalBkgs,
          ExternalBkgFiles,
          ExternalBkgNames,
          ExternalBkgColors,
          AcceptedTime,
          ExternalBkgNGenMC,
          &ExternalBkgActivity[0][0],
          &ExternalBkgEfficiency[0][0],
          allExternals);

#endif
#if RADON_ON
    std::cout << "scale: radon" << std::endl;
    // TODO: for Rn222BkgNames, 3 values are missing
    scale(myFile,
          "../include/activities.txt",
          "externals/",
          nRn222Bkgs,
          Rn222BkgFiles,
          Rn222BkgNames,
          Rn222BkgColors,
          AcceptedTime,
          Rn222BkgNGenMC,
          &Rn222BkgActivity[0][0],
          &Rn222BkgEfficiency[0][0],
          allRn222Bkgs);

    scale(myFile,
          "../include/activities.txt",
          "externals/",
          nRn220Bkgs,
          Rn220BkgFiles,
          Rn220BkgNames,
          Rn220BkgColors,
          AcceptedTime,
          Rn220BkgNGenMC,
          &Rn220BkgActivity[0][0],
          &Rn220BkgEfficiency[0][0],
          allRn220Bkgs);
#endif

    // myFile->Close();

    std::cout << "Stack" << std::endl;
    // hAllExternals, hAllRadon, hAllInternals, hAllNeighbours, hAllNd150 are
    // passed into stackfunction() as argument "hAllSamples"
    // stackfunction calls Add() on object hAllSamples to add each histogram
    // component from each "category" (externals/radon/internals/neighbours
    // /Nd150) to hAllSamples
    // each histogram is pulled from TObjArray allSamples, which is an argument
    // passed in from one of allExternals, allRadon, allInternals,
    // AllNeighbours, allNd150
    // 
    
    // single histogram containing sum of all sub-samples, created using
    // TH1::Add()
    TH1F *hAllExternals[numHistograms];
    TH1F *hAllRadon[numHistograms];
    TH1F *hAllInternals[numHistograms];
    TH1F *hAllNeighbours[numHistograms];
    TH1F *hAllNd150[numHistograms];

    TH1F *hAllMC[numHistograms];
    // for all MC samples to calculate chi2
    // also used to draw black line around stack of all MC samples

    // stacks for drawing: hMajorStacks: samples are grouped by type,
    // eg: radon, external, internal, neighbour, 150Nd
    // hMinorStacks: samples are not grouped, each individual sample
    // is displayed seperatly
    THStack *hMajorStacks[numHistograms];
    THStack *hMinorStacks[numHistograms];

    // I don't know why this is here
    //dir->cd("singleHistos/");
    
    // -1 because we don't want the 2D hist
    // TODO: think this -1 is wrong, because there are several 2d hists
    //for(int i = 0; i < numHistograms-1; i++)
    for(int i = 0; i < numHistograms; i++)
    {
        TString i_str;
        i_str.Form("%i", i);

        hMajorStacks[i] = new THStack("majorstack" + i_str, histogramNames[i]);
        hMinorStacks[i] = new THStack("minorstack" + i_str, histogramNames[i]);



#if EXTERNALS_ON
        if(i == 0)
        {
            std::cout << "stack: externals" << std::endl;
        }
        stackfunction(i, nExternalBkgs,
                      hMinorStacks, allExternals, // TObjArray *allExternals[nExternalBkgs]
                      hMajorStacks, "External Backgrounds",
                      "Externals", hAllExternals, // TH1F *hAllExt[numHistograms]
                      ExternalBkgNames, ExternalBkgFiles, hAllMC,
                      true, false,
                      ExternalBkgColor,
                      true);
#endif
#if RADON_ON
        if(i == 0)
        {
            std::cout << "stack: radon" << std::endl;
        }
        stackfunction(i, nRn222Bkgs,
                      hMinorStacks, allRn222Bkgs,
                      hMajorStacks, "Radon Backgrounds",
                      "Radon", hAllRadon,
                      Rn222BkgNames, Rn222BkgFiles, hAllMC,
                      false, false,
                      RadonBkgColor,
                      false);

        stackfunction(i, nRn220Bkgs,
                      hMinorStacks, allRn220Bkgs,
                      hMajorStacks, "Radon Backgrounds",
                      "Radon", hAllRadon,
                      Rn220BkgNames, Rn220BkgFiles, hAllMC,
                      false, true,
                      RadonBkgColor,
                      true);

#endif
#if INTERNALS_ON
        if(i == 0)
        {
            std::cout << "stack: internals" << std::endl;
        }
        stackfunction(i, nInternalBkgs,
                      hMinorStacks, allInternals,
                      hMajorStacks, "Internal Backgrounds",
                      "Internals", hAllInternals,
                      InternalBkgNames, InternalBkgFiles, hAllMC,
                      false, false,
                      InternalBkgColor,
                      true);
#endif
#if NEIGHBOURS_ON
        if(i == 0)
        {
            std::cout << "stack: neighbours" << std::endl;
        }
        stackfunction(i, nNeighbours,
                      hMinorStacks, allNeighbours,
                      hMajorStacks, "Neighbouring Foils",
                      "Neighbours", hAllNeighbours,
                      NeighbourNames, NeighbourFiles, hAllMC,
                      false, false,
                      NeighbourColor,
                      true);
#endif
#if ND_ON
        if(i == 0)
        {
            std::cout << "stack: 150Nd" << std::endl;
        }
        stackfunction(i, nNd150Samples,
                      hMinorStacks, allNd150Samples,
                      hMajorStacks, "^{150}Nd 2#nu2#beta",
                      "Nd150", hAllNd150,
                      Nd150Names, Nd150Files, hAllMC,
                      false, false,
                      Nd150Color,
                      true);
#endif

        //TString events;

        // TODO: move these into the stackfunction() function
        // TODO: cannot be done, because not on a per-histogram basis
        // yes it is
        #if EXTERNALS_ON
            //std::cout << "Externals: Add()" << std::endl;
            //events.Form("%i",(int)hAllExternals[i]->Integral());
            //hAllExternals[i]->SetTitle("Externals ("+events+")");
            //hMajorStacks[i]->Add((TH1F*)hAllExternals[i]);
        #endif

        #if RADON_ON
            //std::cout << "Radon: Add()" << std::endl;
            //events.Form("%i",(int)hAllRadon[i]->Integral());
            //hAllRadon[i]->SetTitle("Radon ("+events+")");
            //hMajorStacks[i]->Add((TH1F*)hAllRadon[i]);
        #endif

        #if INTERNALS_ON
            //std::cout << "Internals: Add()" << std::endl;
            //events.Form("%i",(int)hAllInternals[i]->Integral());
            //hAllInternals[i]->SetTitle("Internals ("+events+")");
            //hMajorStacks[i]->Add((TH1F*)hAllInternals[i]);
        #endif

        #if NEIGHBOURS_ON
            //std::cout << "Neighbours: Add()" << std::endl;
            //events.Form("%i",(int)hAllNeighbours[i]->Integral());
            //hAllNeighbours[i]->SetTitle("Neighbouring foils ("+events+")");
            //hMajorStacks[i]->Add((TH1F*)hAllNeighbours[i]);
        #endif

        #if ND_ON
            //std::cout << "150Nd: Add()" << std::endl;
            //events.Form("%i",(int)hAllNd150[i]->Integral());
            //hAllNd150[i]->SetTitle("^{150}Nd 2#nu2#beta ("+events+")");
            //hMajorStacks[i]->Add((TH1F*)hAllNd150[i]);
        #endif

        
    }

    std::cout << "Completed the numHistograms loop" << std::endl;

    // stack functions for 2D histograms
/*
    TH2F *hAllMC2D[numHistograms]; //for all MC samples to calculate chi2
    TH2F *tmp2D;
    for(int i = numHistograms-1; i < numHistograms; i++)
    {
        //stack doesnt support 3D hist
        TString i_str;
        i_str.Form("%i",i);

        #if EXTERNALS_ON
        for( int j = 0; j < nExternalBkgs; j++)
        {
            (TH2F*)allExternals[j]->At(i)->Write();
            if(j == 0)
            {
                // TODO: this one different because of clone option
                hAllMC2D[i] = (TH2F*)allExternals[j]->At(i)->Clone("allMC"+histogramNames[i]);
            }
            else
            {
                hAllMC2D[i]->Add((TH2F*)allExternals[j]->At(i)->Clone());
            }
        }//externals
        #endif
    
        #if RADON_ON
        for(int j = 0; j < nRn222Bkgs; j++)
        {
            (TH2F*)allRn222Bkgs[j]->At(i)->Write();
            hAllMC2D[i]->Add((TH2F*)allRn222Bkgs[j]->At(i)->Clone());
        }//radon222

        for(int j = 0; j < nRn220Bkgs; j++)
        {
            (TH2F*)allRn220Bkgs[j]->At(i)->Write();
            hAllMC2D[i]->Add((TH2F*)allRn220Bkgs[j]->At(i)->Clone());
        }//radon220
        #endif

        #if INTERNALS_ON
        for(int j = 0; j < nInternalBkgs; j++)
        {
            (TH2F*)allInternals[j]->At(i)->Write();
            hAllMC2D[i]->Add((TH2F*)allInternals[j]->At(i)->Clone());
        }//internals
        #endif

        #if NEIGHBOURS_ON
        for(int j = 0; j < nNeighbours; j++)
        {
            (TH2F*)allNeighbours[j]->At(i)->Write();
            hAllMC2D[i]->Add((TH2F*)allNeighbours[j]->At(i)->Clone());
        }//neighbours
        #endif

        #if ND_ON
        for(int j = 0; j < nNd150Samples; j++)
        {
            (TH2F*)allNd150Samples[j]->At(i)->Write();
            hAllMC2D[i]->Add((TH2F*)allNd150Samples[j]->At(i)->Clone());
        }//150Nd
        #endif

        hAllMC2D[i]->Write();

        tmp2D = (TH2F*)allDataHistograms->At(i);
        tmp2D->SetName(histogramNames[i]+DataName);
        tmp2D->Write();
    }
    std::cout << "Finished the cloning loop" << std::endl;
*/



    // loop to write out histograms before running canvas loop
    // also sets hist structure for use in following loop
    TH1F* hist[numHistograms];
    for(int i = 0; i < numHistograms; i++)
    {

        hist[i] = (TH1F*)allDataHistograms->At(i);
        hist[i]->SetName(histogramNames[i] + DataName);

        // inserted from below loop (2)
        TString Ndata_str;
        Ndata_str.Form("%i", (int)hist[i]->Integral());
        hist[i]->SetTitle("Data (" + Ndata_str + ")");
        hist[i]->SetLineWidth(2);
        hist[i]->SetMarkerStyle(20);

        // 2020-04-02: removed
        // 2020-04-02: re-enabled, this is required to write the data
        // histogram to file in the root directory
        // 2020-04-02: disabled again, will read data histogram from
        // /singleHistos/unscaled directory in fitting code
        // 2020-04-02: requires fiddling around with GetListOfKeys() and
        // cd() to different TDirectory, can't be bothered with that
        // so will re-enable
        std::cout << "writing histogram with name " << hist[i]->GetName() << " at end of code which was removed" << std::endl;
        hist[i]->Write();
        // TODO: don't like that other functions write out the equivalent
        // MC histograms to file in a different place from this Write() call
        // expecially since this call to Write() is inside a block of code
        // which creates a TCanvas - move this and SetName() elsewhere

    }


    TCanvas *c[numHistograms];
    TLegend *leg[numHistograms];
    TPad *c0;
    TPad *c1;
    TPad *c2;
    std::cout << "About to start the canvas loop" << std::endl;
    // for ( int i = 2; i < 3; i++) {
    //for( int i = 2; i < 4; i++)
    for(int j = 0; j < numHistograms; j++)
    {
        
        if(histogramDrawFlag[j] == 0) continue;

        int i = j;

        std::cout << "At start of canvas loop" << endl;

        TString i_str;
        TString Nmc_str;
        
        i_str.Form("%i", i);
        c[i] = new TCanvas("c" + i_str + "_");
        c[i]->SetFillColor(kWhite);
        // c[i]->SetGrid(1,1);
        //c[i]->SetTicks(1,1);
        //c[i]->Divide(2,2);

        //c[i]->cd(1);
        //c0 = (TPad*)c[i]->GetPad(1);
        // c1_1->SetBottomMargin(0.1);
        //c0->SetGrid(1,1);
        //c0->SetTicks(1,1);
        c[i]->SetGrid(1, 1);
        c[i]->SetTicks(1, 1);

        // moved to above (1)

        Nmc_str.Form("%i", (int)hAllMC[i]->Integral());
        // moved to above (2)
        hist[i]->Draw("PE");

        //std::cout << " data: " << hist[i]->Integral() << std::endl;
        //hSimpleStacks[i]->Draw("hist same");
        hMinorStacks[i]->Draw("hist same");
        //hMajorStacks[i]->Draw("hist same");

        hAllMC[i]->SetLineWidth(2);
        hAllMC[i]->SetLineColor(kBlack);
        hAllMC[i]->SetFillColor(kWhite);
        hAllMC[i]->SetFillStyle(0);
        hAllMC[i]->Sumw2();
        hAllMC[i]->Draw("hist sames");

        double chi2;
        int ndf, igood;
        TString chi2_str, ndf_str;
        double prob = hist[i]->Chi2TestX(hAllMC[i], chi2, ndf, igood, "UW");
        chi2_str.Form("%4.3f", chi2);
        ndf_str.Form("%i", ndf);
        hAllMC[i]->SetTitle("total MC (" + Nmc_str + ")");
        // TODO: moved Ndata_str out to above loop, consider moving out Nmc_str

        //hAllMC[i]->Write();
        // 2020-03-30: removed

        // TODO: forget difference between hAllMC and hMinorStacks
        // why is only 1 written to file? does this write all the
        // histograms to file? cannot be the case, as some
        // histos would be skipped by the continue statement above

        //leg[i] = c0->BuildLegend();
        leg[i] = c[i]->BuildLegend();
        leg[i]->SetName("leg" + i_str + "_");
        leg[i]->SetFillColor(kWhite);
        leg[i]->AddEntry((TObject*)0, "#chi^{2}/ndf = " + chi2_str + "/" + ndf_str, "");
        //c[i]->cd(2);
        TLegend *tmpLeg = (TLegend*)leg[i]->Clone("leg_copy");
        tmpLeg->Draw();
        //c0->cd();
        leg[i]->Delete();
        hist[i]->SetTitle(histogramNames[i]); //switch to titles once you've written them
        hist[i]->Draw("PEsame");

        //c[i]->cd(3);
        //c2 = (TPad*)c[i]->GetPad(3);
        // c1_1->SetBottomMargin(0.1);
        //c2->SetGrid(1,1);
        //c2->SetTicks(1,1);

        //make ratio plot
        //TH1F *ratio = (TH1F*)hist[i]->Clone("ratio_" + histogramNames[i]);
        //ratio->Sumw2();
        //ratio->Divide((TH1F*)allMC[i]);
        //ratio->SetLineWidth(2);
        //ratio->SetMarkerStyle(7);
        //ratio->SetTitle("");
        //ratio->GetYaxis()->SetTitle("data / MC");
        //ratio->GetYaxis()->CenterTitle();

        // ratio->GetYaxis()->SetTitleSize(.1);
        // ratio->GetYaxis()->SetLabelSize(.11);
        // ratio->GetYaxis()->SetNdivisions(9+(100*5));
        // ratio->GetXaxis()->SetTitleSize(.15);

        //  ratio->GetXaxis()->SetTitleOffset(.8);
        // ratio->GetXaxis()->SetLabelOffset(.01);
        //    ratio->GetXaxis()->SetLabelSize(.11);
        //ratio->Draw("PE");
        std::cout << "Iteration of canvasy loop finished" << std::endl;
   

/*
    TCanvas *c[numHistograms];
    TH1F* hist[numHistograms];
    TLegend *leg[numHistograms];
    TPad *c0;
    TPad *c1;
    TPad *c2;
    std::cout << "About to start the canvasy loop" << std::endl;
    // for ( int i = 2; i < 3; i++) {
    for( int i = 2; i < 4; i++)
    {
        std::cout << "At start of canvasy loop" << endl;

        TString i_str, Ndata_str, Nmc_str;
        i_str.Form("%i",i);
        c[i] = new TCanvas("c"+i_str+"_");
        c[i]->SetFillColor(kWhite);
        // c[i]->SetGrid(1,1);
        //c[i]->SetTicks(1,1);
        c[i]->Divide(2,2);

        c[i]->cd(1);
        c0 = (TPad*)c[i]->GetPad(1);
        // c1_1->SetBottomMargin(0.1);
        c0->SetGrid(1,1);
        c0->SetTicks(1,1);

        hist[i] = (TH1F*)allDataHistograms->At(i);
        hist[i]->SetName(histogramNames[i]+DataName);
        hist[i]->Write();

        Ndata_str.Form("%i",(int)hist[i]->Integral());
        Nmc_str.Form("%i",(int)allMC[i]->Integral());
        hist[i]->SetTitle("Data ("+Ndata_str+")");
        hist[i]->SetLineWidth(2);
        hist[i]->SetMarkerStyle(20);
        hist[i]->Draw("PE");

        //Justin printing out integrals:
        std::cout << "Justin printing out integrals:" << std::endl;
        std::cout << "Data integral: " << hist[i]->Integral() << endl;
        std::cout << "MC integral: " << allMC[i]->Integral() << endl;
        
        std::cout << " data: " << hist[i]->Integral() << std::endl;
        //hSimpleStacks[i]->Draw("hist same");
        hMinorStacks[i]->Draw("hist same");

        allMC[i]->SetLineWidth(2);
        allMC[i]->SetLineColor(kBlack);
        allMC[i]->SetFillColor(kWhite);
        allMC[i]->SetFillStyle(0);
        allMC[i]->Sumw2();
        allMC[i]->Draw("hist sames");

        double chi2;
        int ndf, igood;
        TString chi2_str, ndf_str;
        double prob = hist[i]->Chi2TestX(allMC[i], chi2, ndf, igood, "UW");
        chi2_str.Form("%4.3f", chi2);
        ndf_str.Form("%i", ndf);
        allMC[i]->SetTitle("total MC (" + Nmc_str + ")");

        allMC[i]->Write();

        leg[i] = c0->BuildLegend();
        leg[i]->SetName("leg" + i_str + "_");
        leg[i]->SetFillColor(kWhite);
        leg[i]->AddEntry((TObject*)0, "#chi^{2}/ndf = " + chi2_str + "/" + ndf_str, "");
        c[i]->cd(2);
        TLegend *tmpLeg = (TLegend*)leg[i]->Clone("leg_copy");
        tmpLeg->Draw();
        c0->cd();
        leg[i]->Delete();
        hist[i]->SetTitle(histogramNames[i]); //switch to titles once you've written them
        hist[i]->Draw("PEsame");

        c[i]->cd(3);
        c2 = (TPad*)c[i]->GetPad(3);
        // c1_1->SetBottomMargin(0.1);
        c2->SetGrid(1,1);
        c2->SetTicks(1,1);

        //make ratio plot
        TH1F *ratio = (TH1F*)hist[i]->Clone("ratio_" + histogramNames[i]);
        ratio->Sumw2();
        ratio->Divide((TH1F*)allMC[i]);
        ratio->SetLineWidth(2);
        ratio->SetMarkerStyle(7);
        ratio->SetTitle("");
        ratio->GetYaxis()->SetTitle("data / MC");
        ratio->GetYaxis()->CenterTitle();

        // ratio->GetYaxis()->SetTitleSize(.1);
        // ratio->GetYaxis()->SetLabelSize(.11);
        // ratio->GetYaxis()->SetNdivisions(9+(100*5));
        // ratio->GetXaxis()->SetTitleSize(.15);

        //  ratio->GetXaxis()->SetTitleOffset(.8);
        // ratio->GetXaxis()->SetLabelOffset(.01);
        //    ratio->GetXaxis()->SetLabelSize(.11);
        ratio->Draw("PE");
        std::cout << "Iteration of canvasy loop finished" << std::endl;

    */
    }
    std::cout << "Finished with canvas" << std::endl;
    // drawPlots(hEe_data, hEe_totalMC, hs_hEe);


    // myFile->Close();
}


void drawPlots(TH1F *data,TH1F *mc, THStack *hs)
{
    TString whichHist = mc->GetName();
    gStyle->SetTitleFillColor(kWhite);
    gROOT->SetStyle("Plain");

    TFile *myFile = TFile::Open("final"+whichHist+"_2e_P"+Phase+".root", "RECREATE");  

    std::cout << whichHist << std::endl;
    gStyle->SetOptStat(0);
    TCanvas *c1 = new TCanvas("c1_"+whichHist,"c1_"+whichHist,700,700);

    c1->SetFillColor(kWhite);
    c1->SetGrid(1,1);
    c1->Divide(2,2);

    c1->cd(1);
    TPad *c1_1;
    c1_1 = (TPad*)c1->GetPad(1);
    // c1_1->SetBottomMargin(0.1);
    c1_1->SetGrid(1,1);
    c1_1->SetTicks(1,1);

    data->SetMarkerStyle(20);
    data->SetLineWidth(2);
    data->SetTitle("");
    data->GetXaxis()->SetLabelOffset(2.);
    data->GetYaxis()->SetLabelSize(0.05);
    data->Draw("PE");

    hs->Draw("hist same");
    mc->SetLineWidth(2);
    mc->Draw("same");

    c1_1->Update();
    // TLegend *leg;
    // leg = c1_1->BuildLegend();
    // leg->SetName("leg_"+whichHist);
    // leg->SetFillColor(kWhite);
    data->Draw("PEsame");

    c1->cd(2);
    TPad *c1_2;
    c1_2 = (TPad*)c1->GetPad(2);
    TLegend *leg;
    leg = c1_1->BuildLegend();
    leg->SetName("leg_"+whichHist);
    leg->SetFillColor(kWhite);
    leg->Draw();

    c1->cd(3);
    TPad *c1_3;
    c1_3 = (TPad*)c1->GetPad(3);
    //  c1_2->SetTopMargin(0.001);
    // c1_2->SetBottomMargin(2.);
    c1_3->SetGrid(1,1);
    c1_3->SetTicks(1,1);
    double ymin;// = c1_1->GetY1();
    double ymax;// = c1_1->GetY2();
    double xmin;// = c1_1->GetX1();
    double xmax;// = c1_1->GetX2();

    c1_3->GetPadPar(xmin, ymin, xmax, ymax);
    c1_3->SetPad(xmin, ymin + ((ymax - ymin) / 4.), xmax, ymax - ((ymax - ymin) / 4.));
    // c1_2->SetBottomMargin(3.);

    TH1F *ratio = (TH1F*)data->Clone("ratio_" + whichHist); 
    ratio->Sumw2();
    ratio->Divide(mc);
    ratio->SetLineWidth(1);
    ratio->SetMarkerStyle(7);
    ratio->SetTitle("");
    ratio->GetYaxis()->SetTitle("data / MC");
    ratio->GetYaxis()->CenterTitle();
    ratio->GetYaxis()->SetTitleOffset(.2);
    ratio->GetYaxis()->SetTitleSize(.1);
    ratio->GetYaxis()->SetLabelSize(.11);
    ratio->GetYaxis()->SetNdivisions(9 + (100 * 5));
    ratio->GetXaxis()->SetTitleSize(.15);
    ratio->GetXaxis()->SetTitleOffset(.8);
    ratio->GetXaxis()->SetLabelOffset(.01);
    ratio->GetXaxis()->SetLabelSize(.11);
    ratio->Draw("PE");

    c1_3->Update();

    //  TDirectory* dir = gDirectory;
    //  dir->cd();
    //TODO this is wrong
    data->SetTitle("Phase 1OR2");
    data->SetName("data_" + whichHist);
    data->Write();
    //  ratio->Write();
    hs->SetName("hs_" + whichHist);
    hs->Write();
    leg->SetName("leg_" + whichHist);
    leg->Write();
    mc->SetName("mcSum_" + whichHist);
    mc->Write();

    nHistograms++;

}


double getActErr(double Npass, double Ngen, double Ndata, double sf_err)
{
    double err = TMath::Sqrt( (1.0 / Npass) + (1.0 / Ngen) + (1.0 / Ndata) + (sf_err * sf_err) );
    return err;
}

double getChi2(TH1F *data, TH1F* mc, Int_t &ndof)
{

    size_t N_bins =  data->GetXaxis()->GetNbins();
    double chi2 = 0.;
    ndof = 0;

    for(size_t i=0; i<N_bins ; i++)
    {

        double nData = data->GetBinContent(i + 1) ;

        double nMC = mc->GetBinContent(i + 1);

        if(nData + nMC == 0) continue;

        ndof ++;

        double diff = (nMC - nData)/sqrt(nMC + nData); // gaussian prob

        chi2 += diff*diff;
    }

    //  cout << " chi2 " << chi2 << " ndof " << ndof << endl;
    //fflush(stdout);
    std::cout.flush();

    double P = TMath::Prob(chi2,ndof);

    //  cout << " Prob " << P  << endl;
    //fflush(stdout);
    std::cout.flush();

    return chi2;
}

