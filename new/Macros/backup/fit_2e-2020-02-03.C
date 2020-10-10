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
    TFile *myFile = TFile::Open("Nd150_2e_P"+Phase+".root","RECREATE");
    TDirectory *dir = gDirectory;
    TDirectory *dir_sh = dir->mkdir("singleHistos"); 
    dir_sh->mkdir("unscaled");
    myFile->Close();

    // Read in all files and make the cuts you want.
    //First we will read in the data
    std::cout << "Reading in data file...";
    makeHistograms("betabeta/", DataFile);
    std::cout << "... DONE." << std::endl;

    //Now we'll read in all the backgrounds: externals, Rn222, Rn220, internals and neighbouring foils
    std::cout << "Reading in external backgrounds...";
    for(int i = 0; i < nExternalBkgs; i++)
    {
        makeHistograms("externals/", ExternalBkgFiles[i]);
    }
    std::cout << "... DONE." << std::endl;

    std::cout << "Reading in Rn220 backgrounds...";
    for(int i = 0; i < nRn220Bkgs; i++)
    {
        makeHistograms("externals/", Rn220BkgFiles[i]);
    }
    std::cout << "... DONE." << std::endl;

    std::cout << "Reading in Rn222 backgrounds...";
    for(int i = 0; i < nRn222Bkgs; i++)
    {
        makeHistograms("externals/", Rn222BkgFiles[i]);
    }
    std::cout << "... DONE." << std::endl;

    std::cout << "Reading in internal backgrounds...";
    for(int i = 0; i < nInternalBkgs; i++)
    {
        makeHistograms("internals/", InternalBkgFiles[i]);
    }
    std::cout << "... DONE." << std::endl;

    std::cout << "Reading in Nd150 backgrounds...";
    for(int i = 0; i < nNd150Samples; i++)
    {
        makeHistograms("nd150/", Nd150Files[i]);
    }
    std::cout << "... DONE." << std::endl;

    std::cout << "Reading in neighbouring foil backgrounds...";
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
                                       8000, 1000, 8000);

    hNElectrons             = new TH1F("hNElectrons_" + histName,
                                       "Phase " + Phase + " " + histName + " N tracks; N tracks/event",
                                       20, -0.5, 19.5);

    hNUnassocPromptGgHits   = new TH1F("hNUnassocPromptGgHits_" + histName,
                                       "Phase " + Phase + " " + histName + " N unassoc. prompt gg hits; N unassoc prompt gg hits",
                                       20, -0.5, 19.5);

    hTotalE                 = new TH1F("hTotalE_" + histName,
                                       "Phase " + Phase + " " + histName + " total energy; #SigmaE_{e} (MeV)",
                                       50, 0.0, 2.2);

    hInternalPullee         = new TH1F("hInternalPullee_" + histName,
                                       "Phase " + Phase + " " + histName + " internal hypothesis; ee Pull",
                                       50, -40., 40.);

    hInternalProbee         = new TH1F("hInternalProbee_" + histName,
                                       "Phase " + Phase + " " + histName + " internal hypothesis; ee Probability",
                                       50, 0.0, 1.);

    hExternalPullee         = new TH1F("hExternalPullee_" + histName,
                                       "Phase " + Phase + " " + histName + " external hypothesis; ee Pull",
                                       50, -40.,40.);

    hExternalProbee         = new TH1F("hExternalProbee_" + histName,
                                       "Phase " + Phase + " " + histName + " external hypothesis; ee Probability",
                                       50, 0.0, 1.);

    hCosee                  = new TH1F("hCosee_" + histName,
                                       "Phase " + Phase + " " + histName + " angular correlation of electron tracks; cos(#theta)_{ee}",
                                       50, -1., 1.);

    hCoseeWeighted          = new TH1F("hCoseeWeighted_" + histName,
                                       "Phase " + Phase + " " + histName + "corrected angular correlation of electron tracks; cos(#theta)_{ee}",
                                       50, -1., 1.);

    hVertexDZ               = new TH1F("hVertexDZ_" + histName,
                                       "Phase " + Phase + " " + histName + "distance between recon vertices; #DeltaZ (cm)",
                                       50, -30., 30.);

    hVertexDRPhi            = new TH1F("hVertexDRPhi_" + histName,
                                       "Phase " + Phase + " " + histName + "distance between recon vertices; #DeltaR#phi (cm*rad)",
                                       50,-30.,30.);//

    hEeMax                  = new TH1F("hEeMax_" + histName,
                                       "Phase " + Phase + " " + histName + " higher energy Ee; E_{e} (MeV)"              ,
                                       50, 0.0, 4);

    hElectronLengthMax      = new TH1F("hElectronLengthMax_" + histName,
                                       "Phase " + Phase + " " + histName + " higher energy e^{-} track length;  track length (cm)",
                                       50, 0, 600);

    hVertexZMax             = new TH1F("hVertexZMax_" + histName,
                                       "Phase " + Phase + " " + histName + " higher energy vertex position;  Z (cm)"                  ,
                                       50, -120, 120);

    hVertexSectorMax        = new TH1F("hVertexSectorMax_" + histName,
                                       "Phase " + Phase + " " + histName + " higher energy vertex position;  sector"                 ,
                                       50, 5.7, 5.9);

    hVertexRMax             = new TH1F("hVertexRMax_" + histName,
                                       "Phase " + Phase + " " + histName + " higher energy vertex position;  R (cm)"                  ,
                                       50, 154.7, 155.);

    hElectronFirstGgMax     = new TH1F("hElectronFirstGgMax_" + histName,
                                       "Phase " + Phase + " " + histName + " higher energy e^{-} first gg hit location;first gg layer",
                                       9, -0.5, 8.5);

    hElectronLastGgMax      = new TH1F("hElectronLastGgMax_" + histName,
                                       "Phase " + Phase + " " + histName + " higher energy e^{-} last gg hit location;last gg layer",
                                       9, -0.5, 8.5);

    hVertexMinDistPromptGgMax  = new TH1F("hVertexMinDistPromptGgMax_" + histName,
                                       "Phase " + Phase + " " + histName + " distance from higher energy vertex to closest, prompt unassoc. gg hit;(v_{x}, v_{y}, v_{z}) - (gg_{x},gg_{y},gg_{z}) (cm)",
                                       50, 0, 600);

    hElectronLDCorrMax      = new TH1F("hElectronLDCorrMax_" + histName,
                                       "Phase "+Phase+" "+histName+" higher energy e^{-} LD correction;LD corr", 
                                       50, 0.9,1.1);

    hElectronCaloZMax       = new TH1F("hElectronCaloZMax_" + histName,
                                       "Phase "+Phase+" "+histName+" higher energy e^{-} calo hit position; Z(cm)",
                                       50, -150, 150);

    hElectronCaloSectorMax  = new TH1F("hElectronCaloSectorMax_" + histName,
                                       "Phase "+Phase+" "+histName+" higher energy e^{-} calo hit position; sector",
                                       50, 5.7, 5.9);

    hElectronCaloRMax       = new TH1F("hElectronCaloRMax_" + histName,
                                       "Phase "+Phase+" "+histName+" higher energy e^{-} calo hit position; R (cm)",
                                       50, 90, 220);

    hVertexZSecMax          = new TH2F("hVertexZSecMax_" + histName,
                                       "Phase "+Phase+" "+histName+" higher energy vertex location; sector; Z (cm)",
                                       100, 5.7, 5.9,
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
                                       50, 154.7, 155.);

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

    hElectronCaloZMin       = new TH1F("hElectronCaloZMin_" + histName,
                                       "Phase "+Phase+" "+histName+" lower energy e^{-} calo hit position; Z(cm)"          ,
                                       50, -150, 150);

    hElectronCaloSectorMin  = new TH1F("hElectronCaloSectorMin_" + histName,
                                       "Phase "+Phase+" "+histName+" lower energy e^{-} calo hit position; sector"  ,
                                       50, 5.7, 5.9);

    hElectronCaloRMin       = new TH1F("hElectronCaloRMin_" + histName,
                                       "Phase "+Phase+" "+histName+" lower energy e^{-} calo hit position; R (cm)"  ,
                                       50, 90, 220);

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



    hNAfterCuts      = new TH1F("hNAfterCuts_" + histName,
                                       "Phase " + Phase + " " + histName + " Analysis Cut Flow; ; N events",
                                       50, 0, 50); //never know how many we'll make... :)

    hEffAfterCuts      = new TH1F("hEffAfterCuts_" + histName,
                                       "Phase " + Phase + " " + histName + " Analysis Cut Flow; ; efficiency",
                                       50, 0, 50); //never know how many we'll make... :)



    TFile *aFile = TFile::Open(filePath+thePath+histName+"/Nd150_2eNg_output.root");
    TTree *theTree = (TTree*)aFile->Get("Nd150_2eNg/Nd150_2eNg");

    // general
    int         event, run;
    double      eventTime;
    bool        mcFlag;
    //radon related stuff
    double      radonWeight, bi210Weight;

    int         electronSide[2], foilSide, trueVertexLayer;
    //counters 
    int         nGammas, nGammaClusters, nElectrons, nParticles;
    double      electronEnergy[2], eTrackLength[2];
    int         electronFirstGG[2], electronLastGG[2];
    double      trackSign[2];
    double      vertexZ[2], vertexSector[2], vertexR[2];
    double      electronCaloZ[2], electronCaloSector[2], electronCaloR[2];
    int         electronPMT[2], electronBlockType[2];
    int         electronLDFlag[2];
    double      electronLDCorr[2], electronLDCorrErr[2];

    int         nUnAssocPromptGg;
    double      dist2VertexPromptGg0[20];
    double      dist2VertexPromptGg1[20];
    double      eeInternalPull;
    double      eeInternalProb;
    double      eeExternalPull;
    double      eeExternalProb;
    double      cosee;
    double      cosee_weight;

    int         nLowEClusters;
    double      lowEClusterE[10];
    double      lowEClusterSector[10];
    double      lowEClusterZ[10];
    double      lowEClusterR[10];
    double      lowEMinDistPromptGg[10];
    int         lowEClusterPMTs[10];

    theTree->SetBranchAddress( "Run"               , &run );
    theTree->SetBranchAddress( "Event"             , &event );
    theTree->SetBranchAddress( "eventTime"         , &eventTime );
    theTree->SetBranchAddress( "nElectrons"        , &nElectrons );  
    theTree->SetBranchAddress( "nGammaClusters"    , &nGammaClusters );  
    theTree->SetBranchAddress( "nGammaClusters"    , &nGammas );  
    //  theTree->SetBranchAddress( "nParticles"        , &nParticles );  

    theTree->SetBranchAddress( "radonWeight"       , &radonWeight);
    theTree->SetBranchAddress( "bi210Weight"       , &bi210Weight );
    theTree->SetBranchAddress( "trueVertexLayer"   , &trueVertexLayer );
    theTree->SetBranchAddress( "electronSide"      , electronSide );
    theTree->SetBranchAddress( "electronEnergy"    , electronEnergy );
    theTree->SetBranchAddress( "eTrackLength"      , eTrackLength );
    theTree->SetBranchAddress( "firstGgHitLayer"   , electronFirstGG );
    theTree->SetBranchAddress( "lastGgHitLayer"    , electronLastGG);
    theTree->SetBranchAddress( "NAPromptGgHits"    , &nUnAssocPromptGg);
    theTree->SetBranchAddress( "NAPromptGgHitsDist2Vertex"  , dist2VertexPromptGg0);
    theTree->SetBranchAddress( "NAPromptGgHitsDist2Calo"    , dist2VertexPromptGg1);

    theTree->SetBranchAddress( "trackSign"          , trackSign);
    theTree->SetBranchAddress( "vertexZ"            , vertexZ );
    theTree->SetBranchAddress( "vertexSec"          , vertexSector );
    theTree->SetBranchAddress( "vertexR"            , vertexR );
    // theTree->SetBranchAddress( "electronCaloZ"     , electronCaloZ );
    // theTree->SetBranchAddress( "electronCaloR"     , electronCaloR );
    //  theTree->SetBranchAddress( "electronCaloSector", electronCaloSector );
    //  theTree->SetBranchAddress( "electronPMT"       , electronPMT );
    theTree->SetBranchAddress( "electronBlockType"  , electronBlockType );
    theTree->SetBranchAddress( "electronLDFlag"     , electronLDFlag );
    theTree->SetBranchAddress( "electronLDCorr"     , electronLDCorr );
    theTree->SetBranchAddress( "electronLDCorrErr"  , electronLDCorrErr );

    theTree->SetBranchAddress( "internalPull"       , &eeInternalPull );
    theTree->SetBranchAddress( "internalProb"       , &eeInternalProb );
    theTree->SetBranchAddress( "externalPull"       , &eeExternalPull );
    theTree->SetBranchAddress( "externalProb"       , &eeExternalProb );
    theTree->SetBranchAddress( "cosee"              , &cosee );
    theTree->SetBranchAddress( "cosee_weight"       , &cosee_weight );

    theTree->SetBranchAddress( "nTotalClusterHits"  , &nLowEClusters );
    theTree->SetBranchAddress( "clusterHitEnergy"   , lowEClusterE );
    theTree->SetBranchAddress( "clusterHitPMT"      , lowEClusterPMTs );
    theTree->SetBranchAddress( "clusterHitSec"      , lowEClusterSector );
    theTree->SetBranchAddress( "clusterHitZ"        , lowEClusterZ );
    // theTree->SetBranchAddress( "lowEClusterR"      , lowEClusterR );
    // theTree->SetBranchAddress( "lowEMinDistPromptGg", lowEMinDistPromptGg );

    std::cout <<"Processing: " << histName << std::endl;


    int nCuts = 26;
    Int_t cut_counter[nCuts];
    Int_t cut_counter_index = 0;
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
        // count number of events before cuts
        cut_counter_index = 0;
        cut_counter[cut_counter_index]++;


        // theTree->SetBranchAddress( "electronCaloZ"     , electronCaloZ );
        // theTree->SetBranchAddress( "electronCaloR"     , electronCaloR );
        //  theTree->SetBranchAddress( "electronCaloSector", electronCaloSector );
        electronCaloZ[0] = 0;
        electronCaloZ[1] = 0;
        electronCaloR[0] = 0;
        electronCaloR[1] = 0;
        electronCaloSector[0] = 0;
        electronCaloSector[1] = 0;
        
        for(Int_t i = 0; i < 10; ++i)
        {
            lowEClusterR[i] = 0;
            lowEMinDistPromptGg[i] = 0;
        }

        double weight = 1.;
        if(histName.CompareTo("bi214_swire") == 0)
        {
            weight = radonWeight;
        }
        else if(histName.CompareTo("bi210_swire") == 0)
        {
            weight = bi210Weight;
        }
        
        /*    else if ( histName.CompareTo("bi210_swire") == 0 ) {weight = bi210Weight * pow(0.5,eventTime/433123.);}
        else if ( (histName.CompareTo("bi210_swire") == 0) || (histName.CompareTo("bi210_sscin") == 0) || (histName.CompareTo("bi210_sfoil") == 0)) {weight = bi210Weight * pow(0.5,eventTime/433123.);}
        else if ((histName.CompareTo("co60_cuTower")==0) || (histName.CompareTo("co60_steel")==0) || (histName.CompareTo("co60_muMetal")==0) || (histName.CompareTo("co60_cuPetals")==0) )weight = pow(0.5,(eventTime/166344192.));
        else if ( (histName.CompareTo("eu152_sscin")==0) || (histName.CompareTo("eu152_int")==0)) {weight = pow(0.5,(eventTime/4.265e8));}
        else if ((histName.CompareTo("eu154_int")==0)) {weight = pow(0.5,(eventTime/2.71169e8));}
        */

        int cc = 0;
        cut_counter[cc ++] ++;

        // Two electrons from 150Nd foil
        if(nElectrons != 2) continue;
        cut_counter[cc ++]++;

        int highE_index, lowE_index;
        if ( electronEnergy[0] > electronEnergy[1] ) {
          highE_index = 0;
          lowE_index = 1;
        } else {
          highE_index = 1;
          lowE_index = 0;
        }

        // No other tracks in event
        
        // No alpha candidates near vertex
        
        // E_e >= 0.2 MeV

        // <= 20 total free scintillators and <= 10 gamma clusters

        // All PMTs have good status

        // Isolated electron calorimeter hits

        // Electron hits in different scintillator blocks

        // Impacts on block faces
        
        // All PMT have good laser survey and counting rates identified by LD and HS flags
        
        // Tracks start in Gg layer 0 or 1
        
        // Tracks miss <= 1 Gg layer before impact on calorimeter
        
        // Each L_e > 20 cm
        
        // <= 5 unassociated prompt Gg hits
        
        // Elecron block types != L2 / L3
        
        // No gamma rays in event
        
        // Tracks have negative curvature
        
        // Each L_e >= 30 cm
        
        // Tracks not from hot spots
        
        // <= 1 unassociated, prompt Gg hit within 15 cm in XY coordinate of event vertex
        
        // No unassociated, prompt Gg hits on opposite side of foil if tracks are on same side
        
        // Both tracks have hit in Gg L0
        
        // P_{int} >= 0.01 and P_{ext} <= 0.01
        
        // Delta R <= 4cm and Delta Z <= 8 cm
        
        // Each E_e >= 300 keV




        cut_counter[cc ++]++;
        cut_counter[cc ++]++;
        cut_counter[cc ++]++;
        cut_counter[cc ++]++;
        cut_counter[cc ++]++;
        cut_counter[cc ++]++;
        cut_counter[cc ++]++;
        cut_counter[cc ++]++;
        cut_counter[cc ++]++;
        cut_counter[cc ++]++;


        double rPhi0 = (vertexR[0]*vertexSector[0]*TMath::TwoPi()/20.);
        double rPhi1 = (vertexR[1]*vertexSector[1]*TMath::TwoPi()/20.);


/****
        //if (nClusters > 0 ) continue;
        if(nGammaClusters > 0) continue;
        cut_counter[1]++;
      
        //if (nParticles > 2 ) continue; //this should be redundant.
        if (false) continue; //this should be redundant.
        cut_counter[2]++;
       
        int highE_index, lowE_index;
        if ( electronEnergy[0] > electronEnergy[1] ) {
          highE_index = 0;
          lowE_index = 1;
        } else {
          highE_index = 1;
          lowE_index = 0;
        }

        //Make sure calo hits aren't in the petals
        // if ( (electronBlockType[0] > 2) || (electronBlockType[1] > 2) ) continue;
        cut_counter[3]++;
        
        if (fabs(vertexZ[0] >117.) || fabs(vertexZ[1] >117.)) continue;
        cut_counter[4]++;

        cut = false;
        //If electrons are both on the same foil side, require only one to have hits in the first few gg layers. 
        //If they are on opposite sides, require both.
        // if ( electronSide[0] == electronSide[1] ) {
          //     if ( (electronFirstGG[0] > 1) && (electronFirstGG[1] > 1) ) cut = true;
        // } else if ( electronSide[0] != electronSide[1] ) {
        if ( (electronFirstGG[0] > 1) || (electronFirstGG[1] > 1) ) cut = true;
          // }

        if ( cut ) continue;
        cut_counter[5]++;

        if ( (electronLDFlag[0] > 0 ) || (electronLDFlag[1] > 0) ) continue;
        cut_counter[6]++;

        if ( (electronEnergy[0] < 0.2 ) || (electronEnergy[1] < 0.2) ) continue;
        cut_counter[7]++;

        if ( (eTrackLength[0] < 30. ) || (eTrackLength[1] < 30.) ) continue;
        cut_counter[8]++;

    // if ( (fabs(vertexZ[0] - vertexZ[1] ) > 4 )   ) continue;
    //   cut_counter[8]++;

        double rPhi0 = (vertexR[0]*vertexSector[0]*TMath::TwoPi()/20.);
        double rPhi1 = (vertexR[1]*vertexSector[1]*TMath::TwoPi()/20.);

    // if ( (fabs(vertexR[0] - vertexR[1] ) > 2 )   ) continue;

        if ((trackSign[0] > 0) || (trackSign[1] > 0 )) continue;

        if ( eeInternalProb < 0.04 ) continue;
        cut_counter[9]++;
      
        if ( eeExternalProb > 0.01 ) continue;
        cut_counter[10]++;

        // this cut currently isn't passing any events
        // according to SBlott thesis
        // 1869 - 3395 is P1
        // 3396 - 9186 is P2
        cut = false;
        if ( thePhase == 0 ) {
          if ( run > 3395) cut = true;
        } else if (thePhase == 1 ) {
          if ( run <= 3395) cut = true;
        } else cut = false;

        if ( cut ) continue;
        cut_counter[11]++;
****/
/*
        //if ( nTracks > 2 ) continue;
        if ( nElectrons > 2 ) continue;
        cut_counter[0]++;

        //if (nClusters > 0 ) continue;
        if(nGammaClusters > 0) continue;
        cut_counter[1]++;
      
        //if (nParticles > 2 ) continue; //this should be redundant.
        if (false) continue; //this should be redundant.
        cut_counter[2]++;
       
        int highE_index, lowE_index;
        if ( electronEnergy[0] > electronEnergy[1] ) {
          highE_index = 0;
          lowE_index = 1;
        } else {
          highE_index = 1;
          lowE_index = 0;
        }

        //Make sure calo hits aren't in the petals
        // if ( (electronBlockType[0] > 2) || (electronBlockType[1] > 2) ) continue;
        cut_counter[3]++;
        
        if (fabs(vertexZ[0] >117.) || fabs(vertexZ[1] >117.)) continue;
        cut_counter[4]++;

        cut = false;
        //If electrons are both on the same foil side, require only one to have hits in the first few gg layers. 
        //If they are on opposite sides, require both.
        // if ( electronSide[0] == electronSide[1] ) {
          //     if ( (electronFirstGG[0] > 1) && (electronFirstGG[1] > 1) ) cut = true;
        // } else if ( electronSide[0] != electronSide[1] ) {
        if ( (electronFirstGG[0] > 1) || (electronFirstGG[1] > 1) ) cut = true;
          // }

        if ( cut ) continue;
        cut_counter[5]++;

        if ( (electronLDFlag[0] > 0 ) || (electronLDFlag[1] > 0) ) continue;
        cut_counter[6]++;

        if ( (electronEnergy[0] < 0.2 ) || (electronEnergy[1] < 0.2) ) continue;
        cut_counter[7]++;

        if ( (eTrackLength[0] < 30. ) || (eTrackLength[1] < 30.) ) continue;
        cut_counter[8]++;

    // if ( (fabs(vertexZ[0] - vertexZ[1] ) > 4 )   ) continue;
    //   cut_counter[8]++;

        double rPhi0 = (vertexR[0]*vertexSector[0]*TMath::TwoPi()/20.);
        double rPhi1 = (vertexR[1]*vertexSector[1]*TMath::TwoPi()/20.);

    // if ( (fabs(vertexR[0] - vertexR[1] ) > 2 )   ) continue;

        if ((trackSign[0] > 0) || (trackSign[1] > 0 )) continue;

        if ( eeInternalProb < 0.04 ) continue;
        cut_counter[9]++;
      
        if ( eeExternalProb > 0.01 ) continue;
        cut_counter[10]++;

        // this cut currently isn't passing any events
        // according to SBlott thesis
        // 1869 - 3395 is P1
        // 3396 - 9186 is P2
        cut = false;
        if ( thePhase == 0 ) {
          if ( run > 3395) cut = true;
        } else if (thePhase == 1 ) {
          if ( run <= 3395) cut = true;
        } else cut = false;

        if ( cut ) continue;
        cut_counter[11]++;
*/
/*        
        if(nElectrons != 2) continue;

        // leave this here so that intial number of events before cuts
        // is printed
        // note that this will not cut any events, since file contains
        // only 2e events
        cut_counter[cut_counter_index++] ++;

        if(false) continue; // number of other tracks = 0
        cut_counter[cut_counter_index ++] ++;

        if(false) continue; // no alpha candidates near vertex
        cut_counter[cut_counter_index ++] ++;

        if((electronEnergy[0] < 0.2 ) || (electronEnergy[1] < 0.2)) continue;
        cut_counter[cut_counter_index ++] ++;
        //cut_counter[7]++;

        if(false) continue; // <= 20 total free scintillators and <= 10 gamma clusters
        if(nLowEClusters > 10) continue;
        //if(nTotalClusterHits > 10) continue;
        cut_counter[cut_counter_index ++] ++;

        if(false) continue; // all PMTs have good status
        cut_counter[cut_counter_index ++] ++;

        if(false) continue; // isolated electron calorimeter hits
        cut_counter[cut_counter_index ++] ++;

        if(false) continue; // electron hits in different scintillator blocks
        cut_counter[cut_counter_index ++] ++;

        if(false) continue; // impacts on block faces
        cut_counter[cut_counter_index ++] ++;

        if(false) continue; // all PMTs have good laser survery and counter rates identified by LD and HS flags
        cut_counter[cut_counter_index ++] ++;
           
        if(!( (electronFirstGG[0] == 0) || (electronFirstGG[0] == 1) )) continue; // tracks start in Gg layer 0 or 1
        if(!( (electronFirstGG[1] == 0) || (electronFirstGG[1] == 1) )) continue; // tracks start in Gg layer 0 or 1
        cut_counter[cut_counter_index ++] ++;

        if(false) continue; // tracks miss <= 1 gg layer before impact on calorimeter
        cut_counter[cut_counter_index ++] ++;

        //if ( (eTrackLength[0] < 30. ) || (eTrackLength[1] < 30.) ) continue;
        if((eTrackLength[0] <= 20. ) || (eTrackLength[1] <= 20.)) continue;
        cut_counter[cut_counter_index ++] ++;
        //cut_counter[8]++;     

        if(nUnAssocPromptGg > 5) continue;
        //if (nGammaClusters > 0 ) continue;
        cut_counter[cut_counter_index ++] ++;

        //Make sure calo hits aren't in the petals
        // if ( (electronBlockType[0] > 2) || (electronBlockType[1] > 2) ) continue;
        if((electronBlockType[0] == 4) || (electronBlockType[0] == 5)) continue;
        if((electronBlockType[1] == 4) || (electronBlockType[1] == 5)) continue;
        cut_counter[cut_counter_index ++] ++;

        //nParticles = nElectrons+nGammaClusters;
        //if (nParticles > 2 ) continue; //this should be redundant.
        if (nGammas > 0 ) continue; //this should be redundant.
        cut_counter[cut_counter_index ++] ++;
        //cut_counter[2]++;

        // tracks have negative curvature
        if ((trackSign[0] > 0) || (trackSign[1] > 0 )) continue;
        cut_counter[cut_counter_index ++] ++;

        // another length cut, this time longer
        if ( (eTrackLength[0] < 30. ) || (eTrackLength[1] < 30.) ) continue;
        //if ( (eTrackLength[0] <= 20. ) || (eTrackLength[1] <= 20.) ) continue;
        cut_counter[cut_counter_index ++] ++;
        //cut_counter[8]++;     

        // tracks not from hot spots
        if(false) continue;
        cut_counter[cut_counter_index ++] ++;

        // <= 1 unassociated, prompt Gg hit within 15cm in XY coordinate of event vertex
        if(false) continue;
        cut_counter[cut_counter_index ++] ++;

        cut = false;
        //If electrons are both on the same foil side, require only one to have hits in the first few gg layers. 
        //If they are on opposite sides, require both.
        // if ( electronSide[0] == electronSide[1] ) {
          //     if ( (electronFirstGG[0] > 1) && (electronFirstGG[1] > 1) ) cut = true;
        // } else if ( electronSide[0] != electronSide[1] ) {
        //if ( (electronFirstGG[0] > 1) || (electronFirstGG[1] > 1) ) cut = true;
          // }
        if( electronSide[0] == electronSide[1] )
        {
            if(nUnAssocPromptGg > 0) continue;
        }
        // TODO: this will not work because need to check nUnAssocPromptGg hits are on the opposite side
        cut_counter[cut_counter_index ++] ++;
        //if ( cut ) continue;
        //cut_counter[cut_counter_index++]++;
        //cut_counter[5]++;

        // both tracks have hits in Gg L0
        if(false) continue;
        cut_counter[cut_counter_index ++] ++;


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

//        cut_counter[cut_counter_index++]++;
//        //cut_counter[3]++;

//        if (fabs(vertexZ[0] >117.) || fabs(vertexZ[1] >117.)) continue;
//        cut_counter[cut_counter_index++]++;
//        //cut_counter[4]++;


//        if ( (electronLDFlag[0] > 0 ) || (electronLDFlag[1] > 0) ) continue;
//        cut_counter[cut_counter_index++]++;
//        //cut_counter[6]++;


        // if ( (fabs(vertexZ[0] - vertexZ[1] ) > 4 )   ) continue;
        //   cut_counter[8]++;

        double rPhi0 = (vertexR[0] * vertexSector[0] * TMath::TwoPi() / 20.);
        double rPhi1 = (vertexR[1] * vertexSector[1] * TMath::TwoPi() / 20.);

        // if ( (fabs(vertexR[0] - vertexR[1] ) > 2 )   ) continue;


        //if ( eeInternalProb < 0.04 ) continue;
        if(eeInternalProb < 0.01) continue;
        //cut_counter[cut_counter_index++]++;
        //cut_counter[9]++;
        if(eeExternalProb > 0.01) continue;
        //cut_counter[10]++;
        cut_counter[cut_counter_index ++] ++;

        if(std::abs(vertexR[0] - vertexR[1]) > 4.) continue;
        if(std::abs(vertexZ[0] - vertexZ[1]) > 8.) continue;
        cut_counter[cut_counter_index ++] ++;

        //
        if( (electronEnergy[0] < 0.3) || (electronEnergy[1] < 0.3) ) continue;
        cut_counter[cut_counter_index ++] ++;

//        cut = false;
//        if ( thePhase == 0 ) {
//          if ( run > 3395) cut = true;
//        } else if (thePhase == 1 ) {
//          if ( run <= 3395) cut = true;
//        } else cut = false;
//
//        if ( cut ) continue;
//        //cut_counter[11]++;
//        cut_counter[cut_counter_index++]++;




        for(int i = 0; i < nUnAssocPromptGg; i++)
        {
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

        //if ( nLowEClusters > 0 )continue;
*/


        hRun->Fill(run, weight);
        hNElectrons->Fill(nElectrons, weight);
        hTotalE->Fill(electronEnergy[0] + electronEnergy[1], weight);
        // TODO: why are these being filled with different weight
        hInternalPullee->Fill(eeInternalPull, 1.0);
        hInternalProbee->Fill(eeInternalProb, 1.0);
        hExternalPullee->Fill(eeExternalPull, 1.0);
        hExternalProbee->Fill(eeExternalProb, 1.0);
        hCosee->Fill(cosee, weight);
        hCoseeWeighted->Fill(cosee, weight * cosee_weight);
        hVertexDZ->Fill(vertexZ[highE_index] - vertexZ[lowE_index], weight);
        hVertexDRPhi->Fill(rPhi0 - rPhi1, weight);

        hEeMax->Fill(electronEnergy[highE_index], weight);
        hElectronLengthMax->Fill(eTrackLength[highE_index], weight);
        hVertexSectorMax->Fill(vertexSector[highE_index], weight);
        hVertexZMax->Fill(vertexZ[highE_index], weight);
        hVertexRMax->Fill(vertexR[highE_index], weight);
        hVertexZSecMax->Fill(vertexSector[highE_index], vertexZ[highE_index], weight);
        hElectronFirstGgMax->Fill(electronFirstGG[highE_index], weight);
        hElectronLastGgMax->Fill(electronLastGG[highE_index], weight);
        hElectronLDCorrMax->Fill(electronLDCorr[highE_index], weight);
        hElectronCaloSectorMax->Fill(electronCaloSector[highE_index], weight);
        hElectronCaloZMax->Fill(electronCaloZ[highE_index], weight);
        hElectronCaloRMax->Fill(electronCaloR[highE_index], weight);

        hEeMin->Fill(electronEnergy[lowE_index], weight);
        hElectronLengthMin->Fill(eTrackLength[lowE_index], weight);
        hVertexSectorMin->Fill(vertexSector[lowE_index], weight);
        hVertexZMin->Fill(vertexZ[lowE_index], weight);
        hVertexRMin->Fill(vertexR[lowE_index], weight);
        hVertexZSecMin->Fill(vertexSector[lowE_index], vertexZ[lowE_index], weight);
        hElectronFirstGgMin->Fill(electronFirstGG[lowE_index], weight);
        hElectronLastGgMin->Fill(electronLastGG[lowE_index], weight);
        hElectronLDCorrMin->Fill(electronLDCorr[lowE_index], weight);
        hElectronCaloSectorMin->Fill(electronCaloSector[lowE_index], weight);
        hElectronCaloZMin->Fill(electronCaloZ[lowE_index], weight);
        hElectronCaloRMin->Fill(electronCaloR[lowE_index], weight);

        hEeMaxVEeMin->Fill(electronEnergy[highE_index], electronEnergy[lowE_index], weight);

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
    hVertexDRPhi->Write();
    hNUnassocPromptGgHits->Write();

    hEeMax->Write();
    hElectronLengthMax->Write();
    hVertexZMax->Write();
    hVertexSectorMax->Write();
    hVertexRMax->Write();
    hElectronFirstGgMax->Write();
    hElectronLastGgMax->Write();
    hVertexMinDistPromptGgMax->Write();
    hElectronLDCorrMax->Write();
    hElectronCaloSectorMax->Write();
    hElectronCaloZMax->Write();
    hElectronCaloRMax->Write();
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
    hElectronCaloSectorMin->Write();
    hElectronCaloZMin->Write();
    hElectronCaloRMin->Write();
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
void scale(TFile* myFile, const int nData, TString *dataFiles, TString *dataNames, std::string typedir, double *dataNGenMC, std::string activityfile, double *dataActivity, double *dataEfficiency, Color_t *dataColor, Color_t *dataColors, TObjArray** allData, double *AcceptedTime, double additional_scaling_factor)
{
    std::cout << "scale" << std::endl;

    std::ofstream oftemp("oftemp1.txt", std::ofstream::out | std::ofstream::app);

    TH1F *tmpHist; //we'll be using this later.

    //Let's get the internals/externals/data and add them to an array
    //TObjArray *allData[nData];
    for(int i = 0; i < nData; i++)
    {
        // create new object array
        std::cout << dataFiles[i] << std::endl;
        allData[i] = new TObjArray();

        // read in the information needed for calculating efficiencies
        // read in number of generated MC events
        std::ifstream inFile;
        inFile.open(filePath+typedir+dataFiles[i]+"/JobSummary.txt");
        // cout << "JobSummary stuff is here: "
        // 	 << filePath << "internals/" << InternalBkgFiles[i] << "/JobSummary.txt"
        // 	 << endl;
        std::string dummy;
        inFile >> dummy >> dataNGenMC[i];
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
            
            inFile2 >> mcname >> activityPhase1 >> activityPhase2;
            std::cout << mcname << ", " << activityPhase1 << ", " << activityPhase2 << std::endl;
            // cout << "Read in from activity file: " << mcname
            // 	   << ", " << activityPhase1
            // 	   << ", " << activityPhase2 << endl;
            // cout << "Trying to compare to " << InternalBkgFiles[i] << endl;
            
            if(mcname.CompareTo(dataFiles[i]) == 0)
            {
                //dataActivity[0][i] = activityPhase1;
                //dataActivity[1][i] = activityPhase2;
                dataActivity[0 * nData + i] = activityPhase1;
                dataActivity[1 * nData + i] = activityPhase2;
            }
        }

        inFile2.close();

        for(int j = 0; j < numHistograms; j++)
        {
            tmpHist = (TH1F*)myFile->Get("singleHistos/unscaled/"+histogramNames[j]+dataFiles[i])->Clone(histogramNames[j]+dataFiles[i]+"_fit");
            // cout << "The histogram I have got is "
            // 	   << "singleHistos/unscaled/"
            // 	   << histogramNames[j]
            // 	   << InternalBkgFiles[i]
            // 	   << endl;
            // cout << "At this point, the integral is " << tmpHist->Integral() << endl;

            //calculate efficiency
            if(j == 0)
            {
                double NPass = tmpHist->Integral();
                std::cout << "NPass: " << NPass << std::endl;
                
                // TODO: don't fully understand this
                // efficiency = (number of events which pass cuts * total time) / 
                // (number of generated events * accepted time)
                // accepted time is like total time - dead time?
                //double eff = NPass / (double)(dataNGenMC[i] * (AcceptedTime[thePhase]/(double)TotalTime) );
                double eff = (double)NPass / (double)dataNGenMC[i];
                double dead_time_scale_factor = (double)TotalTime / (double)AcceptedTime[thePhase];
                double efficiency = eff * dead_time_scale_factor;
                std::cout << "efficiency: " << efficiency << std::endl;

                //dataEfficiency[thePhase][i] = eff;
                dataEfficiency[thePhase * nData + i] = efficiency;
            }
            tmpHist->Sumw2();
            tmpHist->SetFillColor(dataColors[i]);
            tmpHist->SetLineColor(dataColors[i]);
            tmpHist->SetTitle(dataNames[i]);
            // cout << "Integral before scaling: " << tmpHist->Integral() << endl;
            // cout << "InternalBkgActivity[thePhase][i]: " << InternalBkgActivity[thePhase][i] << endl
            // 	   << "AcceptedTime[thePhase]: " << AcceptedTime[thePhase] << endl
            // 	   << "InternalBkgNGenMC[i]: " << InternalBkgNGenMC[i] << endl
            // 	   << "AcceptedTime[thePhase]: " << AcceptedTime[thePhase] << endl
            // 	   << "TotalTime: " << TotalTime << endl;
            //tmpHist->Scale(dataActivity[thePhase][i] * AcceptedTime[thePhase] * 0.001 / (dataNGenMC[i] * (AcceptedTime[thePhase]/(double)TotalTime)) );
            //tmpHist->Scale(dataActivity[thePhase * nData + i] * 0.001 * (double)TotalTime / dataNGenMC[i] );
            //std::cout << "scaling, index=" << thePhase + nData + i << std::endl;
            
            //Double_t scale_factor = dataActivity[thePhase * nData + i] * AcceptedTime[thePhase] * additional_scaling_factor / (dataNGenMC[i] * (AcceptedTime[thePhase]/(double)TotalTime));
            Double_t the_activity = dataActivity[thePhase * nData + i] * additional_scaling_factor;
            Double_t n_expected_events = the_activity * (double)TotalTime;
            Double_t scale_factor = n_expected_events / dataNGenMC[i];
            tmpHist->Scale(scale_factor); 

            //std::cout << "sf=" << dataActivity[thePhase * nData + i] * AcceptedTime[thePhase] * additional_scaling_factor / (dataNGenMC[i] * (AcceptedTime[thePhase]/(double)TotalTime)) << std::endl;
            oftemp << dataFiles[i] << " -> " << dataNames[i] << " : " << histogramNames[j] << " -> sf=" << scale_factor << ", " << dataActivity[thePhase * nData + i] << " integral()=> " << tmpHist->Integral() << std::endl;
            //tmpHist->Scale(dataActivity[thePhase * nData + i] * AcceptedTime[thePhase] / (dataNGenMC[i] * (AcceptedTime[thePhase]/(double)TotalTime)) );
            // cout << "Integral after scaling: " << tmpHist->Integral() << endl;


            //tmpHist->Scale(InternalBkgActivity[thePhase][i] * AcceptedTime[thePhase] * 0.001    / (InternalBkgNGenMC[i] * (AcceptedTime[thePhase]/(double)TotalTime)) );
            //tmpHist->Scale(NeighbourActivity[thePhase][i]   * AcceptedTime[thePhase] * 0.001    / (NeighbourNGenMC[i]   * (AcceptedTime[thePhase]/(double)TotalTime)) );
            //tmpHist->Scale(Nd150Activity[thePhase][i]       * AcceptedTime[thePhase] * 0.001    / (Nd150NGenMC[i]       * (AcceptedTime[thePhase]/(double)TotalTime)) );
            //tmpHist->Scale(ExternalBkgActivity[thePhase][i] * AcceptedTime[thePhase]            / (ExternalBkgNGenMC[i] * (AcceptedTime[thePhase]/(double)TotalTime)) );
            //tmpHist->Scale(ExternalBkgActivity[thePhase][i] * AcceptedTime[thePhase]            / (ExternalBkgNGenMC[i] * (AcceptedTime[thePhase]/(double)TotalTime)) );
            //tmpHist->Scale(Rn222BkgActivity[thePhase][i]    * AcceptedTime[thePhase]            / (Rn222BkgNGenMC[i]    * (AcceptedTime[thePhase]/(double)TotalTime)) );
            //tmpHist->Scale(Rn220BkgActivity[thePhase][i]    * AcceptedTime[thePhase]            / (Rn220BkgNGenMC[i]    * (AcceptedTime[thePhase]/(double)TotalTime)) );



            allData[i]->Add(tmpHist);

            // cout << "Justin is here. numHistograms is " << numHistograms  << endl;
            // cout << "Integral is " << tmpHist->Integral() << endl;

        }//all histos
    }//internals

    oftemp.close();
}


void stackfunction(int i, int nGeneric, THStack **hAllStacks, TObjArray **allGeneric, TString genericName, TH1F **hAllGeneric, TString *GenericNames, TString *GenericFiles, TH1F **hAllMC, bool allMC_first, bool hAllGeneric_addtoexisting)
{
    // i = histogram index (type, eg; electron energy min/max, total energy, track length...)
    // j = background/data/mc index

    TH1F *tmpHist; //we'll be using this later.

    // generic
    for(int j = 0; j < nGeneric; j++)
    {
        if(hAllGeneric_addtoexisting == false)
        {
            if(j == 0)
            {
                hAllGeneric[i] = (TH1F*)allGeneric[j]->At(i)->Clone(genericName);
                //allMC[i] = (TH1F*)hAllExt[i]->Clone("Total MC");
            }
            else
            {
                hAllGeneric[i]->Add((TH1F*)allGeneric[j]->At(i));
                //allMC[i]->Add((TH1F*)allExternals[j]->At(i));
            }
        }
        else
        {
            hAllGeneric[i]->Add((TH1F*)allGeneric[j]->At(i));
        }

        // allMC goes here
        if(allMC_first == true)
        {
            if(j == 0)
            {
                //allMC[i] = (TH1F*)hAllGeneric[i]->Clone("Total MC");
                // TODO: this is hAllExt - bug?
                hAllMC[i] = (TH1F*)hAllGeneric[j]->Clone("Total MC");
            }
            else
            {
                // TODO: this is allExternals - bug?
                hAllMC[i]->Add((TH1F*)hAllGeneric[j]);
            }
        }
        else
        {
            // TODO: this is allExternals/allRn222Bkgs - bug?
            hAllMC[i]->Add((TH1F*)hAllGeneric[j]);
        }
      
        tmpHist = (TH1F*)allGeneric[j]->At(i);
        if(tmpHist->Integral() > 0)
        {
            //std::cout << ExternalBkgFiles[j] << " " << tmpHist->Integral() << std::endl;
            tmpHist->SetTitle(GenericNames[j]);
            hAllStacks[i]->Add((TH1F*)allGeneric[j]->At(i)->Clone(GenericFiles[j]));
        }
      
        (TH1F*)allGeneric[j]->At(i)->Write();

    } // generic

}


void fitHistograms()
{

    std::cout << "Time to get fitting!" << std::endl;

    TFile *myFile = TFile::Open("Nd150_2e_P"+Phase+".root", "UPDATE");
    TDirectory* dir = gDirectory; 
    dir->ls();
    TObjArray *allDataHistograms = new TObjArray();    
    
    //Get the data histograms
    for ( int i = 0; i < numHistograms; i++ ) {
        std::string name("singleHistos/unscaled/"+histogramNames[i]+DataFile);
        //std::cout << "name=" << name << std::endl;
        allDataHistograms->Add( (TH1F*)myFile->Get(name.c_str())->Clone(histogramNames[i]+"_data") );
    }


    AcceptedTime[0] = 33859178.;
    AcceptedTime[1] = 133770114.;
    AcceptedTime[2] = TotalTime;
                              
    TH1F *tmpHist; //we'll be using this later.



    std::ofstream oftemp("oftemp1.txt", std::ofstream::out | std::ofstream::trunc);
    oftemp.close();

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

    scale(myFile, nInternalBkgs, InternalBkgFiles, InternalBkgNames, "internals/", InternalBkgNGenMC,
        "../include/internalActivities.txt", &InternalBkgActivity[0][0], &InternalBkgEfficiency[0][0],
        InternalBkgColor, InternalBkgColors, allInternals, AcceptedTime, 0.001);

    scale(myFile, nNeighbours, NeighbourFiles, NeighbourNames, "neighbours/", NeighbourNGenMC,
        "../include/internalActivities2.txt", &NeighbourActivity[0][0], &NeighbourEfficiency[0][0],
        NeighbourBkgColors, NeighbourBkgColors, allNeighbours, AcceptedTime, 0.001);

    scale(myFile, nNd150Samples, Nd150Files, Nd150Names, "nd150/", Nd150NGenMC,
        "../include/internalActivities2.txt", &Nd150Activity[0][0], &Nd150Efficiency[0][0],
        Nd150BkgColor, Nd150BkgColor, allNd150Samples, AcceptedTime, 0.001);

    scale(myFile, nExternalBkgs, ExternalBkgFiles, ExternalBkgNames, "externals/", ExternalBkgNGenMC,
        "../include/externalActivities.txt", &ExternalBkgActivity[0][0], &ExternalBkgEfficiency[0][0],
        ExternalBkgColor, ExternalBkgColor, allExternals, AcceptedTime, 1.0);

    // TODO: for Rn222BkgNames, 3 values are missing
    scale(myFile, nRn222Bkgs, Rn222BkgFiles, Rn222BkgNames, "externals/", Rn222BkgNGenMC,
        "../include/radonActivities.txt", &Rn222BkgActivity[0][0], &Rn222BkgEfficiency[0][0],
        Rn222BkgColor, Rn222BkgColor, allRn222Bkgs, AcceptedTime, 1.0);

    scale(myFile, nRn220Bkgs, Rn220BkgFiles, Rn220BkgNames, "externals/", Rn220BkgNGenMC,
        "../include/radonActivities.txt", &Rn220BkgActivity[0][0], &Rn220BkgEfficiency[0][0],
        Rn220BkgColor, Rn220BkgColor, allRn220Bkgs, AcceptedTime, 1.0);


    // myFile->Close();
    std::cout << "Now time to stack em" << std::endl;

    TH1F *hAllExternals[numHistograms];
    TH1F *hAllRadon[numHistograms];
    TH1F *hAllInternals[numHistograms];
    TH1F *hAllNeighbours[numHistograms];
    TH1F *hAllNd150[numHistograms];
    TH1F *hAllMC[numHistograms]; //for all MC samples to calculate chi2
    THStack *hStacks[numHistograms];
    THStack *hAllStacks[numHistograms];
    //std::cout << "Before the cd to singleHistos" << std::endl;
    dir->cd("singleHistos/");
    //std::cout << "After the cd to singleHistos" << std::endl;
    
    // -1 because we don't want the 2D hist
    // TODO: think this -1 is wrong, because there are several 2d hists
    for(int i = 0; i < numHistograms-1; i++)
    {
        TString i_str;
        i_str.Form("%i",i);

        hStacks[i] = new THStack("stack"+i_str,histogramNames[i]);
        hAllStacks[i] = new THStack("allstack"+i_str,histogramNames[i]);

  /*      
        //////////Externals////////////////
        for(int j = 0; j < nExternalBkgs; j++)
        {
        
            if(j == 0)
            {
                hAllExt[i] = (TH1F*)allExternals[j]->At(i)->Clone("Externals");
                allMC[i] = (TH1F*)hAllExt[i]->Clone("Total MC");
            }
            else
            {
                hAllExt[i]->Add((TH1F*)allExternals[j]->At(i));
                allMC[i]->Add((TH1F*)allExternals[j]->At(i));
            }
          
            tmpHist = (TH1F*)allExternals[j]->At(i);
            if(tmpHist->Integral() > 0)
            {
                //std::cout << ExternalBkgFiles[j] << " " << tmpHist->Integral() << std::endl;
                tmpHist->SetTitle(ExternalBkgNames[j]);
                hAllStacks[i]->Add((TH1F*)allExternals[j]->At(i)->Clone(ExternalBkgFiles[j]));
            }
          
            ((TH1F*)allExternals[j])->At(i)->Write();
        }//externals
*/

//        void stackfunction(int i,
//                   int nGeneric,
//                   THStack **hAllStacks,
//                   TObjArray **allGeneric,
//                   TString genericName,
//                   TH1F **hAllGnrc,
//                   TString *GenericNames,
//                   TString *GenericFiles)

        stackfunction(i,
                      nExternalBkgs,
                      hAllStacks,
                      allExternals, // TObjArray *allExternals[nExternalBkgs]
                      "Externals",
                      hAllExternals, // TH1F *hAllExt[numHistograms]
                      ExternalBkgNames,
                      ExternalBkgFiles,
                      hAllMC,
                      true,
                      false);

        /*
        //////////Externals////////////////
        for(int j = 0; j < nExternalBkgs; j++)
        {
        
            if(j == 0)
            {
                ////hAllExt[i] = (TH1F*)allExternals[j]->At(i)->Clone("Externals");
                // TODO: is this a bug? using Clone?
                //
                // this line is different from every other following code block
                allMC[i] = (TH1F*)hAllExt[i]->Clone("Total MC");
            }
            else
            {
                ////hAllExt[i]->Add((TH1F*)allExternals[j]->At(i));
                // TODO: is this a bug, was allExt above, now allExternals
                //
                // this line is different from every other following code block
                allMC[i]->Add((TH1F*)allExternals[j]->At(i));
            }
          
            ////tmpHist = (TH1F*)allExternals[j]->At(i);
            ////if(tmpHist->Integral() > 0)
            ////{
            ////    //std::cout << ExternalBkgFiles[j] << " " << tmpHist->Integral() << std::endl;
            ////    tmpHist->SetTitle(ExternalBkgNames[j]);
            ////    hAllStacks[i]->Add((TH1F*)allExternals[j]->At(i)->Clone(ExternalBkgFiles[j]));
            ////}
          
            //(TH1F*)allExternals[j]->At(i)->Write();
        }//externals
        */


        stackfunction(i,
                      nRn222Bkgs,
                      hAllStacks,
                      allRn222Bkgs,
                      "Radon",
                      hAllRadon,
                      Rn222BkgNames,
                      Rn222BkgFiles,
                      hAllMC,
                      false,
                      false);

        /*
        //////////Rn222////////////////
        for(int j = 0; j < nRn222Bkgs; j++)
        {
            if (j == 0 )
            {
                // done by stackfunction
                //hAllRn[i] = (TH1F*)allRn222Bkgs[j]->At(i)->Clone("Radon");
                // TODO: there is no allMC here?
            }
            else
            {
                // done by stackfunction
                //hAllRn[i]->Add((TH1F*)allRn222Bkgs[j]->At(i));
            }

            // not done by stackfunction
            allMC[i]->Add((TH1F*)allRn222Bkgs[j]->At(i));
            // TODO: but there is an allMC here

            // done by stackfunction
            //tmpHist = (TH1F*)allRn222Bkgs[j]->At(i);
            //if( tmpHist->Integral() > 0 )
            //{
            //    tmpHist->SetTitle(Rn222BkgNames[j]);
            //    hAllStacks[i]->Add((TH1F*)allRn222Bkgs[j]->At(i)->Clone(Rn222BkgFiles[j]));
            //}

            //(TH1F*)allRn222Bkgs[j]->At(i)->Write();
        }//radon222
        */

        stackfunction(i,
                      nRn220Bkgs,
                      hAllStacks,
                      allRn220Bkgs,
                      "Radon",
                      hAllRadon,
                      Rn220BkgNames,
                      Rn220BkgFiles,
                      hAllMC,
                      false,
                      true);

        // this one is different
        //////////Rn220////////////////
        //for(int j = 0; j < nRn220Bkgs; j++)
        //{
        //    hAllRn[i]->Add((TH1F*)allRn220Bkgs[j]->At(i));
        //
        //    allMC[i]->Add((TH1F*)allRn220Bkgs[j]->At(i));
        //
        //    tmpHist = (TH1F*)allRn220Bkgs[j]->At(i);
        //    if(tmpHist->Integral() > 0)
        //    {
        //        tmpHist->SetTitle(Rn220BkgNames[j]);
        //        hAllStacks[i]->Add((TH1F*)allRn220Bkgs[j]->At(i)->Clone(Rn220BkgFiles[j]));
        //    }
        //
        //    (TH1F*)allRn220Bkgs[j]->At(i)->Write();
        //
        //}//rn220

        stackfunction(i,
                      nInternalBkgs,
                      hAllStacks,
                      allInternals,
                      "Internals",
                      hAllInternals,
                      InternalBkgNames,
                      InternalBkgFiles,
                      hAllMC,
                      false,
                      false);

        //////////Internals////////////////
        //for(int j = 0; j < nInternalBkgs; j++)
        //{
        //    if( j == 0 )
        //    {
        //        // done by stackfunction
        //        //hAllInt[i] = (TH1F*)allInternals[j]->At(i)->Clone("Internals");
        //    }
        //    else
        //    {
        //        // done by stackfunction
        //        //hAllInt[i]->Add((TH1F*)allInternals[j]->At(i));
        //    }
        //
        //    allMC[i]->Add((TH1F*)allInternals[j]->At(i));
        //
        //    // done by stackfunction
        //    //tmpHist = (TH1F*)allInternals[j]->At(i);
        //    //if(tmpHist->Integral() > 0)
        //    //{
        //    //    tmpHist->SetTitle(InternalBkgNames[j]);
        //    //    hAllStacks[i]->Add((TH1F*)allInternals[j]->At(i)->Clone(InternalBkgFiles[j]));
        //    //}
        //
        //    //(TH1F*)allInternals[j]->At(i)->Write();
        //}//internals
        
        stackfunction(i,
                      nNeighbours,
                      hAllStacks,
                      allNeighbours,
                      "Neighbours",
                      hAllNeighbours,
                      NeighbourNames,
                      NeighbourFiles,
                      hAllMC,
                      false,
                      false);

        /*
        //////////Neighbours////////////////
        for(int j = 0; j < nNeighbours; j++)
        {
            if(j == 0)
            {
                hAllNeighbours[i] = (TH1F*)allNeighbours[j]->At(i)->Clone("Neighbours");
            }
            else
            {
                hAllNeighbours[i]->Add((TH1F*)allNeighbours[j]->At(i));
            }

            tmpHist = (TH1F*)allNeighbours[j]->At(i);
            if(tmpHist->Integral() > 0)
            {
                tmpHist->SetTitle(NeighbourNames[j]);
                hAllStacks[i]->Add((TH1F*)allNeighbours[j]->At(i)->Clone(NeighbourFiles[j]));
            }

            allMC[i]->Add((TH1F*)allNeighbours[j]->At(i));
            (TH1F*)allNeighbours[j]->At(i)->Write();
        }//neighbouring foils
        */

        stackfunction(i,
                      nNd150Samples,
                      hAllStacks,
                      allNd150Samples,
                      "Nd150",
                      hAllNd150,
                      Nd150Names,
                      Nd150Files,
                      hAllMC,
                      false,
                      false);

    /*
        for(int j = 0; j < nNd150Samples; j++)
        {
            if(j == 0)
            {
                hAllNd150[i] = (TH1F*)allNd150Samples[j]->At(i)->Clone("Nd150");
            }
            else
            {
                hAllNd150[i]->Add((TH1F*)allNd150Samples[j]->At(i));
            }

            tmpHist = (TH1F*)allNd150Samples[j]->At(i);
            if(tmpHist->Integral() > 0)
            {
                tmpHist->SetTitle(Nd150Names[j]);
                hAllStacks[i]->Add((TH1F*)allNd150Samples[j]->At(i)->Clone(Nd150Files[j]));
            }

            allMC[i]->Add((TH1F*)allNd150Samples[j]->At(i));
            (TH1F*)allNd150Samples[j]->At(i)->Write();
            // cout << "Justin:" << endl;
            // cout << "allMC[i] integral: " << allMC[i]->Integral() << endl;
        }//nd150
    */

        TString events;
        events.Form("%i",(int)hAllExternals[i]->Integral());
        hAllExternals[i]->SetTitle("Externals ("+events+")");

        events.Form("%i",(int)hAllRadon[i]->Integral());
        hAllRadon[i]->SetTitle("Radon ("+events+")");

        events.Form("%i",(int)hAllInternals[i]->Integral());
        hAllInternals[i]->SetTitle("Internals ("+events+")");

        events.Form("%i",(int)hAllNeighbours[i]->Integral());
        hAllNeighbours[i]->SetTitle("Neighbouring foils ("+events+")");

        events.Form("%i",(int)hAllNd150[i]->Integral());
        hAllNd150[i]->SetTitle("^{150}Nd 2#nu2#beta ("+events+")");

        hStacks[i]->Add((TH1F*)hAllExternals[i]);
        hStacks[i]->Add((TH1F*)hAllRadon[i]);
        hStacks[i]->Add((TH1F*)hAllInternals[i]);
        hStacks[i]->Add((TH1F*)hAllNeighbours[i]);
        hStacks[i]->Add((TH1F*)hAllNd150[i]);

    }

    std::cout << "Completed the numHistograms loop" << std::endl;
    TH2F *hAllMC2D[numHistograms]; //for all MC samples to calculate chi2
    TH2F *tmp2D;
    for(int i = numHistograms-1; i < numHistograms; i++)
    {
        //stack doesnt support 3D hist
        TString i_str;
        i_str.Form("%i",i);

        for( int j = 0; j < nExternalBkgs; j++)
        {
            (TH2F*)allExternals[j]->At(i)->Write();
            if(j == 0)
            {
                hAllMC2D[i] = (TH2F*)allExternals[j]->At(i)->Clone("allMC"+histogramNames[i]);
            }
            else
            {
                hAllMC2D[i]->Add((TH2F*)allExternals[j]->At(i)->Clone());
            }
        }//externals

        for(int j = 0; j < nRn222Bkgs; j++)
        {
            (TH2F*)allRn222Bkgs[j]->At(i)->Write();
            hAllMC2D[i]->Add((TH2F*)allRn222Bkgs[j]->At(i)->Clone());
        }//radon222

        for(int j = 0; j < nRn220Bkgs; j++)
        {
            (TH2F*)allRn220Bkgs[j]->At(i)->Write();
            hAllMC2D[i]->Add((TH2F*)allRn220Bkgs[j]->At(i)->Clone());
        }//rn220

        for(int j = 0; j < nInternalBkgs; j++)
        {
            (TH2F*)allInternals[j]->At(i)->Write();
            hAllMC2D[i]->Add((TH2F*)allInternals[j]->At(i)->Clone());
        }//internals

        for(int j = 0; j < nNeighbours; j++)
        {
            (TH2F*)allNeighbours[j]->At(i)->Write();
            hAllMC2D[i]->Add((TH2F*)allNeighbours[j]->At(i)->Clone());
        }//internals

        for(int j = 0; j < nNd150Samples; j++)
        {
            (TH2F*)allNd150Samples[j]->At(i)->Write();
            hAllMC2D[i]->Add((TH2F*)allNd150Samples[j]->At(i)->Clone());
        }//internals

        hAllMC2D[i]->Write();

        tmp2D = (TH2F*)allDataHistograms->At(i);
        tmp2D->SetName(histogramNames[i]+DataName);
        tmp2D->Write();
    }
    std::cout << "Finished the cloning loop" << std::endl;


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
        //c[i]->Divide(2,2);

        //c[i]->cd(1);
        //c0 = (TPad*)c[i]->GetPad(1);
        // c1_1->SetBottomMargin(0.1);
        //c0->SetGrid(1,1);
        //c0->SetTicks(1,1);
        c[i]->SetGrid(1,1);
        c[i]->SetTicks(1,1);

        hist[i] = (TH1F*)allDataHistograms->At(i);
        hist[i]->SetName(histogramNames[i]+DataName);
        hist[i]->Write();

        Ndata_str.Form("%i",(int)hist[i]->Integral());
        Nmc_str.Form("%i",(int)hAllMC[i]->Integral());
        hist[i]->SetTitle("Data ("+Ndata_str+")");
        hist[i]->SetLineWidth(2);
        hist[i]->SetMarkerStyle(20);
        hist[i]->Draw("PE");

        //Justin printing out integrals:
        std::cout << "Justin printing out integrals:" << std::endl;
        std::cout << "Data integral: " << hist[i]->Integral() << endl;
        std::cout << "MC integral: " << hAllMC[i]->Integral() << endl;
        
        //std::cout << " data: " << hist[i]->Integral() << std::endl;
        //hSimpleStacks[i]->Draw("hist same");
        hAllStacks[i]->Draw("hist same");

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

        hAllMC[i]->Write();

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
        hAllStacks[i]->Draw("hist same");

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
    std::cout << "Finished the canvasy bit" << std::endl;
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

