#include "TH1D.h"
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
#include "TH1D.h"
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
#include <exception>

#include "../include/fit_2e.h"


// TODO:
//
// - change to double precision floating point

//-----------------------------------------
//    Functions
//----------------------------------------
void loadFiles();

void makeHistograms(TString thePath, TString sampleName, std::ofstream &ofile_cutcount, const Int_t mode_flag, const uint64_t sample_split_flags = 0);
//void  makeHistograms(TTree *tree, TString sampleName, Int_t isotope);
void fitHistograms();
void drawPlots(TH1D *data,TH1D *mc, THStack *hs);
double getActErr(double Npass, double Ngen, double Ndata, double sf_err );
double getChi2(TH1D *data, TH1D* mc, Int_t &ndof);

void loadFiles()
{

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);

    //First, create a root file to hold all of the histograms
    TFile *myFile = new TFile("Nd150_2e_P" + Phase + ".root", "RECREATE");
    TDirectory *dir = gDirectory;
    //TDirectory *dir_sh = dir->mkdir("singleHistos"); 
    //dir_sh->mkdir("unscaled");
    
    // raw data, no cuts
    TDirectory *dir_rawdata = dir->mkdir("rawdata");
    // processed data with cuts
    TDirectory *dir_processeddata = dir->mkdir("processeddata");
    // processed data with cuts and scaling
    TDirectory *dir_scaled = dir->mkdir("scaled");

    
    myFile->Close();


    // map of flags, used for selecting events which are split by
    // L0, L>0, IN, OUT
    // 0 = default
    // 1 = select IN events
    // 2 = select OUT events
    // 4 = select L0 IN events
    // 8 = select L0 OUT events
    // 16 = select L>0 (IN+OUT) events
    /*
    std::map<std::string, uint64_t> sample_split_flags;
    sample_split_flags["bi214_sfoil"] = (uint64_t)1 | (uint64_t)2;
    sample_split_flags["pb214_sfoil"] = (uint64_t)1 | (uint64_t)2;
    sample_split_flags["bi214_swire"] = (uint64_t)4 | (uint64_t)8 | (uint64_t)16;
    sample_split_flags["pb214_swire"] = (uint64_t)4 | (uint64_t)8 | (uint64_t)16;
    */


    // open and reset file for storing the number of events that pass cuts
    std::ofstream ofile_cutcount("cutcount.txt", std::ofstream::out | std::ofstream::trunc);
    // TODO: need to count the number of events that pass cuts and store
    // into variable ExternalBkgPassCut[thePhase][i]
    // cannot pass that variable here because this would require
    // loadFiles() to be called before fitHistograms()
    // which is slow
    // instead, store results into a file

    // Read in all files and make the cuts you want.
    //First we will read in the data
    std::cout << ">>>>> Reading in data file...";
    makeHistograms("betabeta/", DataFile, ofile_cutcount, 0);
    makeHistograms("betabeta/", DataFile, ofile_cutcount, 1);
    std::cout << "... DONE." << std::endl;

    // Read in all the backgrounds: externals, Rn222, Rn220, internals and neighbouring foils
    std::cout << std::endl;
    std::cout << ">>>>> Reading in external backgrounds..." << std::endl;
    for(int i = 0; i < nExternalBkgs; i++)
    {
        makeHistograms("externals/", ExternalBkgFiles[i], ofile_cutcount, 0);
        makeHistograms("externals/", ExternalBkgFiles[i], ofile_cutcount, 1);
    }
    std::cout << "... DONE." << std::endl;


    std::cout << std::endl;
    std::cout << ">>>>> Reading in Rn220 backgrounds..." << std::endl;
    for(int i = 0; i < nRn220Bkgs; i++)
    {
        makeHistograms("externals/", Rn220BkgFiles[i], ofile_cutcount, 0);
        makeHistograms("externals/", Rn220BkgFiles[i], ofile_cutcount, 1);
    }
    std::cout << "... DONE." << std::endl;

    #if 0
    std::cout << std::endl;
    std::cout << ">>>>> Reading in Rn222 backgrounds..." << std::endl;
    for(int i = 0; i < nRn222BkgsInput; i++)
    {
        TString sample_name = Rn222BkgFilesInput[i];
        std::string std_string_sample_name = std::string(sample_name.Data());
        if(sample_split_flags.count(std_string_sample_name) > 0)
        {
            //std::cout << "catch: flags: " << std_string_sample_name << std::endl;
            //std::cin.get();

            uint64_t flags = sample_split_flags.at(std_string_sample_name);
            //std::cout << "flags=" << flags << std::endl;
            int counter = 0;
            //               |  ||  ||  ||  |
            while(flags & (0xFFFFFFFFFFFFFFFF << counter))
            {
                //std::cout << "while" << std::endl;

                //uint64_t flag_bit = flags & (uint64_t)0x01;
                //                             |  ||  ||  ||  |
                uint64_t flag_bit = flags & (0x0000000000000001 << counter);
                //std::cout << "flag_bit=" << flag_bit << std::endl;
                if(flag_bit)
                {
                    //uint64_t f = flag_bit >> counter;
                    //std::cout << "calling makeHistograms with flag_bit=" << flag_bit << std::endl;
                    //std::cin.get();
                    makeHistograms("externals/", Rn222BkgFilesInput[i], ofile_cutcount, 0, flag_bit);
                    makeHistograms("externals/", Rn222BkgFilesInput[i], ofile_cutcount, 1, flag_bit);
                }

                //flags_shift << 1;
                ++ counter;

                //std::cout << std::hex;
                //std::cout << "0xFFFFFFFFFFFFFFFF << counter = " << (0xFFFFFFFFFFFFFFFF << counter) << " flags=" << flags << std::endl;
                //std::cout << std::dec;
            }
        }
        else
        {
            makeHistograms("externals/", Rn222BkgFilesInput[i], ofile_cutcount, 0);
            makeHistograms("externals/", Rn222BkgFilesInput[i], ofile_cutcount, 1);
        }

    /*
        if(sample_name == "bi214_sfoil" ||
           sample_name == "pb214_sfoil")
        {
            makeHistograms("externals/", Rn222BkgFiles[i], ofile_cutcount, 0, sample_split_flags[sample_name]);
            makeHistograms("externals/", Rn222BkgFiles[i], ofile_cutcount, 1, sample_split_flags[sample_name]);
        }
        else if(sample_name == "bi214_swire" ||
                sample_name == "pb214_swire")
        {
            makeHistograms("externals/", Rn222BkgFiles[i], ofile_cutcount, 0, sample_split_flags[sample_name]);
            makeHistograms("externals/", Rn222BkgFiles[i], ofile_cutcount, 1, sample_split_flags[sample_name]);
        }
        else
        {
            makeHistograms("externals/", Rn222BkgFiles[i], ofile_cutcount, 0);
            makeHistograms("externals/", Rn222BkgFiles[i], ofile_cutcount, 1);
        }
    */
    }
    std::cout << "... DONE." << std::endl;
    #endif




    std::cout << std::endl;
    std::cout << ">>>>> Reading in Rn222 backgrounds..." << std::endl;
    for(int i = 0; i < nRn222Bkgs; i++)
    {
        makeHistograms("externals/", Rn222BkgFilesNew[i], ofile_cutcount, 0);
        makeHistograms("externals/", Rn222BkgFilesNew[i], ofile_cutcount, 1);
    }
    std::cout << "... DONE." << std::endl;




    std::cout << std::endl;
    std::cout << ">>>>> Reading in internal backgrounds..." << std::endl;
    for(int i = 0; i < nInternalBkgs; i++)
    {
        makeHistograms("internals/", InternalBkgFiles[i], ofile_cutcount, 0);
        makeHistograms("internals/", InternalBkgFiles[i], ofile_cutcount, 1);
    }
    std::cout << "... DONE." << std::endl;


    std::cout << std::endl;
    std::cout << ">>>>> Reading in neighbouring foil backgrounds..." << std::endl;
    for(int i = 0; i < nNeighbours; i++)
    {
        makeHistograms("neighbours/", NeighbourFiles[i], ofile_cutcount, 0);
        makeHistograms("neighbours/", NeighbourFiles[i], ofile_cutcount, 1);
    }
    std::cout << "... DONE." << std::endl;


    std::cout << std::endl;
    std::cout << ">>>>> Reading in Nd150 backgrounds..." << std::endl;
    for(int i = 0; i < nNd150Samples; i++)
    {
        makeHistograms("nd150/", Nd150Files[i], ofile_cutcount, 0);
        makeHistograms("nd150/", Nd150Files[i], ofile_cutcount, 1);
    }
    std::cout << "... DONE." << std::endl;

    std::cout << std::endl;
    std::cout << ">>>>> Reading in Nd150 systematics..." << std::endl;
    for(int i = 0; i < nNd150Systematics; i++)
    {
        makeHistograms("nd150/", Nd150Systematics[i], ofile_cutcount, 0);
        //makeHistograms("nd150/", Nd150Systematics[i], ofile_cutcount, 1);
    }
    std::cout << "... DONE." << std::endl;
    


    ofile_cutcount.close();

}


// mode_flag = 0: run function with cuts enabled (normal), write to (normal)
//                location (/rawdata/...)
// mode_flag = 1: run function with cuts disabled, write to location
//                ...
//
//                TODO / NOTE: currently using same NAMES and pointers
//                to hold histogram data, may wish to change behaviour
//                depending on value of the flag such that the NAMES
//                and pointers are always unique (surpress root warnings)
void makeHistograms(
    TString thePath,
    TString sampleName,
    std::ofstream& ofile_cutcount,
    const Int_t mode_flag = 0,
    const uint64_t sample_split_flags = 0)
{

    // enable / disable additional cuts
    // by default the basic cuts are enabled
    // this includes sector and phase cut
    const Int_t mode_flag_2 = 1;

    // swap between raw mode and processed mode
    // raw mode disables all cuts
    TString name_append;
    if(mode_flag == 0)
    {
        //name_append = "";
    }
    else if(mode_flag == 1)
    {
        name_append = "_raw";
    }


    //Declare the histograms.
    TH1D* hRun; // 0
    TH1D* hNElectrons;
    TH1D* hTotalE;
    TH1D* hInternalPullee;
    TH1D* hInternalProbee; 
    TH1D* hExternalPullee; // 5
    TH1D* hExternalProbee;
    TH1D* hCosee;
    TH1D* hCoseeWeighted;
    TH1D* hVertexDZ;
    TH1D* hVertexDR; // 10
    TH1D* hVertexDRPhi;
    TH1D* hNAPromptGgHits;

    TH1D* hEeMax;
    TH1D* hElectronLengthMax;
    TH1D* hVertexZMax; // 15
    TH1D* hVertexSectorMax;
    TH1D* hVertexRMax;
    TH1D* hElectronFirstGgMax;
    TH1D* hElectronLastGgMax;
    TH1D* hVertexMinDistPromptGgMax; // 20
    TH1D* hElectronLDCorrMax;
    TH2D* hVertexZSecMax; // usefull really only for data.

    TH1D* hEeMin;
    TH1D* hElectronLengthMin;
    TH1D* hVertexZMin; // 25
    TH1D* hVertexSectorMin;
    TH1D* hVertexRMin;
    TH1D* hElectronFirstGgMin;
    TH1D* hElectronLastGgMin;
    TH1D* hVertexMinDistPromptGgMin; // 30
    TH1D* hElectronLDCorrMin;
    TH2D* hVertexZSecMin; // usefull really only for data.

    // TODO: obsolete, delete these
    TH1D* hNLowEGammas;
    TH1D* hLowEGammaEnergy;
    TH1D* hSummedLowEGammaE; // 35
    TH1D *hLowEMinDistPromptGg;

    TH2D *hEeMaxVEeMin;

    //Histograms for keeping track of analysis cuts
    TH1D* hNAfterCuts;
    TH1D* hEffAfterCuts;

    // my new histograms
    //
    //
    // missing items?
    TH1D* hNAPromptGgHitsDist2VertexMin; // 40
    TH1D *hTrackSignMin;
    TH1D *hTrackSignMax;
    TH1D *hnGammaClusters;
    TH1D *hnInCluster;
    TH1D *hclusterHitEnergy; // 45
    TH1D *hclusterHitEnergyMin;
    TH1D *hclusterHitEnergyMax;
    TH1D *hnLowEnergyHits;
    TH1D *hclusterEnergy;

    TH1D *hFoilSide;
    TH1D *hTrueVertexR;
    TH2D *hFoilSideTrueVertexR;
    TH1D *hTrueVertexLayer;

    TH1D *hSingleEnergy;
    TH1D *hHighEnergy;
    TH1D *hLowEnergy;
    TH2D *hHighLowEnergy;

    TH1D *hTotalClusterEnergy;
    TH2D *hHotSpotLocation;

    TH1D *hEnergySum;
    TH1D *hEnergyDiff;

    // NOTE: this idea (multuplying histograms) will not work because
    // we lose the information of which true energy values match which
    // reconstructed energy values
    // processing must be done event by event
    //#if TRUTH_ENABLE
    //TH2D *hMCElectronEnergy;
    //TH2D *hMCTrueElectronEnergy;
    //#endif


    TString name_sample_split_additional = "";
    /*
    if(sample_split_flags == 0x01)
    {
        // IN
        std::cout << "sample_split_flags is 1" << std::endl;
        name_sample_split_additional = "_IN";
    }
    else if(sample_split_flags == 0x02)
    {
        // OUT
        std::cout << "sample_split_flags is 2" << std::endl;
        name_sample_split_additional = "_OUT";
    }
    else if(sample_split_flags == 0x04)
    {
        // L0 IN
        std::cout << "sample_split_flags is 4" << std::endl;
        name_sample_split_additional = "_INL0";
    }
    else if(sample_split_flags == 0x08)
    {
        // L0 OUT
        std::cout << "sample_split_flags is 8" << std::endl;
        name_sample_split_additional = "_OUTL0";
    }
    else if(sample_split_flags == 0x10)
    {
        // L>0 IN+OUT
        std::cout << "sample_split_flags is 16" << std::endl;
        name_sample_split_additional = "_Lg0";
    }
    std::cout << "name_sample_split_additional=" << name_sample_split_additional << std::endl;
    */

    hRun                    = new TH1D("hRun_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " runs; run numbers",
                                       8000, 0, 10000); // 1000, 8000

    hNElectrons             = new TH1D("hNElectrons_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " N tracks; N tracks/event",
                                       20, -0.5, 19.5); // limits? all events are 2 for nd150

    hNAPromptGgHits         = new TH1D("hNAPromptGgHits_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " N unassoc. prompt gg hits; N unassoc prompt gg hits",
                                       20, -0.5, 19.5); // limits? blank histogram for nd150?

    hTotalE                 = new TH1D("hTotalE_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " total energy; #SigmaE_{e} (MeV)",
                                       50, 0.0, 5.0); // 0.0, 2.2
                                       // TODO was 4.0

    hInternalPullee         = new TH1D("hInternalPullee_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " internal hypothesis; ee Pull",
                                       50, -40., 40.); // limits? all events within -10, 10 for nd150

    hInternalProbee         = new TH1D("hInternalProbee_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " internal hypothesis; ee Probability",
                                       50, 0.0, 1.); // limits?

    hExternalPullee         = new TH1D("hExternalPullee_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " external hypothesis; ee Pull",
                                       50, -20., 20.); // limits? all events within -30, 10 for nd150

    hExternalProbee         = new TH1D("hExternalProbee_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " external hypothesis; ee Probability",
                                       50, 0.0, 1.); // limits?

    hCosee                  = new TH1D("hCosee_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " angular correlation of electron tracks; cos(#theta)_{ee}",
                                       50, -1., 1.); // limits?

    hCoseeWeighted          = new TH1D("hCoseeWeighted_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + "corrected angular correlation of electron tracks; cos(#theta)_{ee}",
                                       50, -1., 1.); // limits?

    hVertexDZ               = new TH1D("hVertexDZ_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + "distance between recon vertices; #DeltaZ (cm)",
                                       50, -10., 10.); // limits? all events within -10, 10 for nd150

    hVertexDR               = new TH1D("hVertexDR_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + "distance between recon vertices; #DeltaR (cm)",
                                       50, -0.01, 0.01); // limits?

    hVertexDRPhi            = new TH1D("hVertexDRPhi_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + "distance between recon vertices; #DeltaR#phi (cm*rad)",
                                       50, -5., 5.); // limits? all events witnin -10, 10 for nd150

    hEeMax                  = new TH1D("hEeMax_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " higher energy Ee; E_{e} (MeV)"              ,
                                       50, 0.0, 5.0); // limits ok
                                       // TODO was 4.0

    hElectronLengthMax      = new TH1D("hElectronLengthMax_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " higher energy e^{-} track length;  track length (cm)",
                                       50, 0, 600); // limits?

    hVertexZMax             = new TH1D("hVertexZMax_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " higher energy vertex position;  Z (cm)"                  ,
                                       50, -120, 120); // limits?
                                       //50, -150, 150); // limits?

    hVertexSectorMax        = new TH1D("hVertexSectorMax_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " higher energy vertex position;  sector"                 ,
                                       50, 5.7, 5.9); // limits?
                                       //50, 5.5, 6.1); // limits?

    hVertexRMax             = new TH1D("hVertexRMax_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " higher energy vertex position;  R (cm)"                  ,
                                       //50, 154.7, 155.0); // limits? all events within 154.85, 154.9 for nd150
                                       50, 154.8904, 154.8908); // limits? all events within 154.85, 154.9 for nd150

    hElectronFirstGgMax     = new TH1D("hElectronFirstGgMax_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " higher energy e^{-} first gg hit location;first gg layer",
                                       9, -0.5, 8.5); // all in 0 for nd150

    hElectronLastGgMax      = new TH1D("hElectronLastGgMax_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " higher energy e^{-} last gg hit location;last gg layer",
                                       9, -0.5, 8.5); // all in 7/8 for nd150

    hVertexMinDistPromptGgMax  = new TH1D("hVertexMinDistPromptGgMax_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " distance from higher energy vertex to closest, prompt unassoc. gg hit;(v_{x}, v_{y}, v_{z}) - (gg_{x},gg_{y},gg_{z}) (cm)",
                                       50, 0, 600); // blank?

    hElectronLDCorrMax      = new TH1D("hElectronLDCorrMax_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" higher energy e^{-} LD correction;LD corr", 
                                       50, 0.9,1.1); // limits? don't know what this histogram is


    hVertexZSecMax          = new TH2D("hVertexZSecMax_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" higher energy vertex location; sector; Z (cm)",
                                       100, 5.7, 5.9,
                                       //100, 5.5, 6.1,
                                       100, -120, 120); 

    hEeMin                  = new TH1D("hEeMin_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" lower energy Ee; E_{e} (MeV)",
                                       50, 0.0, 5.0);

    hElectronLengthMin      = new TH1D("hElectronLengthMin_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" lower energy e^{-} track length;  track length (cm)",
                                       50, 0, 600);

    hVertexZMin             = new TH1D("hVertexZMin_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" lower energy vertex position;  Z (cm)"                  ,
                                       50, -120, 120);
                                       //50, -150, 150);

    hVertexSectorMin        = new TH1D("hVertexSectorMin_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" lower energy vertex position;  sector"                 ,
                                       50, 5.7, 5.9);
                                       //50, 5.5, 6.1);

    hVertexRMin             = new TH1D("hVertexRMin_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" lower energy vertex position;  R (cm)"                  ,
                                       //50, 154.7, 155.0);
                                       50, 154.8904, 154.8908);

    hElectronFirstGgMin     = new TH1D("hElectronFirstGgMin_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" lower energy e^{-} first gg hit location;first gg layer",
                                       9, -0.5, 8.5);

    hElectronLastGgMin      = new TH1D("hElectronLastGgMin_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" lower energy e^{-} last gg hit location;last gg layer",
                                       9, -0.5, 8.5);

    hVertexMinDistPromptGgMin   = new TH1D("hVertexMinDistPromptGgMin_" + sampleName + name_sample_split_additional + name_append,
                                           "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " distance from lower energy vertex to closest, prompt unassoc. gg hit;(v_{x},v_{y},v_{z}) - (gg_{x},gg_{y},gg_{z}) (cm)",
                                           50, 0, 600);

    hElectronLDCorrMin      = new TH1D("hElectronLDCorrMin_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" lower energy e^{-} LD correction;LD corr"            ,
                                       50, 0.9, 1.1);

    hVertexZSecMin          = new TH2D("hVertexZSecMin_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" lower energy vertex location; sector; Z (cm)"              ,
                                       150, 5.8904, 5.8907,
                                       100, -120, 120);


    hNLowEGammas            = new TH1D("hNLowEGammas_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" unassoc. scint hits E < 200 keV; N hits"     ,
                                       11, -0.5, 10.5);

    hLowEGammaEnergy        = new TH1D("hLowEGammaEnergy_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" unassoc. scint hits E < 200 keV; Energy (indiv. hits) (MeV)" ,
                                       50, 0, 0.2);

    hSummedLowEGammaE       = new TH1D("hSummedLowEGammaE_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" unassoc. scint hits E < 200 keV; Energy (indiv. hits) (MeV)" ,
                                       50, 0, 0.5);

    hLowEMinDistPromptGg    = new TH1D("hLowEMinDistPromptGg_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" distance from unassoc. calo hits to closest, prompt unassoc. gg hit;(v_{x},v_{y},v_{z}) - (gg_{x},gg_{y},gg_{z}) (cm)",
                                       100, 0, 600.);//all low E hits


    hEeMaxVEeMin = new TH2D("hEeMaxVEeMin_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " electron energies; Ee^{Max} (MeV); Ee^{min} (MeV)",
                                       50, 0, 4, 50, 0, 4);

    hNAPromptGgHitsDist2VertexMin = new TH1D("hNAPromptGgHitsDist2VertexMin_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " min vertex distance to unassociated prompt geiger hit; min vertex distance (cm); y label",
                                       50, 0., 100.);
                                        



    hNAfterCuts      = new TH1D("hNAfterCuts_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Analysis Cut Flow; ; N events",
                                       20, 0, 20); //never know how many we'll make... :)

    hEffAfterCuts      = new TH1D("hEffAfterCuts_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Analysis Cut Flow; ; efficiency",
                                       20, 0, 20); //never know how many we'll make... :)

    hTrackSignMax   = new TH1D("hTrackSignMax_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Track Sign (Max)",
                                3, -1, 1);

    hTrackSignMin   = new TH1D("hTrackSignMin_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Track Sign (Min)",
                                3, -1, 1);

    hnGammaClusters = new TH1D("hnGammaClusters_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Num Gamma Clusters",
                                11, 0, 10);

    hnInCluster     = new TH1D("hnInCluster_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Num in Cluster",
                                11, 0, 10);

    hclusterHitEnergy     = new TH1D("hclusterHitEnergy_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Cluster Hit Energy",
                                50, 0, 0.5);
    
    hclusterHitEnergyMin     = new TH1D("hclusterHitEnergyMin_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Cluster Hit Energy Min",
                                50, 0, 0.5);

    hclusterHitEnergyMax     = new TH1D("hclusterHitEnergyMax_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Cluster Hit Energy Max",
                                50, 0, 0.5);

    hnLowEnergyHits     = new TH1D("hnLowEnergyHits_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Num in Cluster",
                                11, 0, 10);

    hclusterEnergy     = new TH1D("hclusterEnergy_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Cluster Energy",
                                50, 0, 1.0);


    hFoilSide     = new TH1D("hFoilSide_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Foil Side",
                                21, -10, 10);

    hTrueVertexR     = new TH1D("hTrueVertexR_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " True Vertex R",
                                100, -10.0, 400.0);

    hFoilSideTrueVertexR  = new TH2D("hFoilSideTrueVertexR_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Foil Side True Vertex R",
                                21, -10, 10, 100, -10.0, 400.0);

    hTrueVertexLayer     = new TH1D("hTrueVertexLayer_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " True Vertex Layer",
                                20, -5, 14);

    hSingleEnergy     = new TH1D("hSingleEnergy_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Single Electron Energy;Energy (MeV)",
                                50, 0.0, 5.0);

    hHighEnergy     = new TH1D("hHighEnergy_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " High Energy Electron Energy; Energy (MeV)",
                                50, 0.0, 5.0);

    hLowEnergy     = new TH1D("hLowEnergy_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Low Energy Electron Energy;Energy (MeV)",
                                50, 0.0, 5.0);

    hHighLowEnergy     = new TH2D("hHighLowEnergy_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Electron Energy;Low Energy Electron Energy (MeV);High Energy Electron Energy (MeV)",
                                50, 0.0, 5.0, 50, 0.0, 5.0);


    hTotalClusterEnergy     = new TH1D("hTotalClusterEnergy_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Total Cluster Hit Energy",
                                50, 0, 1.0);

    hHotSpotLocation = new TH2D("hHotSpotLocation_" + sampleName + name_sample_split_additional + name_append,
                        "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Hot Spot Locations",
                        50, 5.7371, 5.8706, 50, -120.0, 120.0);

    hEnergySum     = new TH1D("hEnergySum_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Single Electron Energy;Energy (MeV)",
                                50, 0.0, 5.0);
                                
    hEnergyDiff     = new TH1D("hEnergyDiff_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Single Electron Energy;Energy (MeV)",
                                50, 0.0, 5.0);




    //#if TRUTH_ENABLE
    // limits used for both are same as hTotalE limits, since these
    // datasets will be used to reweight 150 Nd MC and create new hTotalE
    // histogram

    // 150 Nd MC electron energy, both electrons in 2D histogram
    //hMCElectronEnergy = new TH2D("hMCElectronEnergy_" + sampleName + name_append,
    //                                "Phase " + Phase + " " + sampleName + name_append + " MC Electron Energy",
    //                                50, 0.0, 4.0);

    // 150 Nd MC true electron energy
    //hMCTrueElectronEnergy = new TH2D("hMCTrueElectronEnergy_" + sampleName + name_append,
    //                                "Phase " + Phase + " " + sampleName + name_append + " MC True Electron Energy",
    //                                50, 0.0, 4.0, 50, 0.0, 4.0);
    //#endif



    ///////////////////////////////////////////////////////////////////////////


    // map for data structure for raw data
    std::map<TString, TH1*> histogramPointers_rawdata;
    // map for data structure for processed data
    std::map<TString, TH1*> histogramPointers;


    std::map<TString, TH1*> *hpmap_p = nullptr;
    if(mode_flag == 0)
    {
        hpmap_p = &histogramPointers;
    }
    else if(mode_flag == 1)
    {
        hpmap_p = &histogramPointers_rawdata;
    }
    std::map<TString, TH1*> &hpmap{*hpmap_p};
    hpmap["hRun_"] = hRun;
    hpmap["hNElectrons_"] = hNElectrons;
    hpmap["hTotalE_"] = hTotalE;
    hpmap["hInternalPullee_"] = hInternalPullee;
    hpmap["hInternalProbee_"] = hInternalProbee;
    hpmap["hExternalPullee_"] = hExternalPullee;
    hpmap["hExternalProbee_"] = hExternalProbee;
    hpmap["hCosee_"] = hCosee;
    hpmap["hCoseeWeighted_"] = hCoseeWeighted;
    hpmap["hVertexDZ_"] = hVertexDZ;
    hpmap["hVertexDR_"] = hVertexDR;
    hpmap["hVertexDRPhi_"] = hVertexDRPhi;
    hpmap["hNAPromptGgHits_"] = hNAPromptGgHits;
    hpmap["hEeMax_"] = hEeMax;
    hpmap["hElectronLengthMax_"] = hElectronLengthMax;
    hpmap["hVertexZMax_"] = hVertexZMax;
    hpmap["hVertexSectorMax_"] = hVertexSectorMax;
    hpmap["hVertexRMax_"] = hVertexRMax;
    hpmap["hElectronFirstGgMax_"] = hElectronFirstGgMax;
    hpmap["hElectronLastGgMax_"] = hElectronLastGgMax;
    hpmap["hVertexMinDistPromptGgMax_"] = hVertexMinDistPromptGgMax;
    hpmap["hElectronLDCorrMax_"] = hElectronLDCorrMax;
    hpmap["hVertexZSecMax_"] = hVertexZSecMax;
    hpmap["hEeMin_"] = hEeMin;
    hpmap["hElectronLengthMin_"] = hElectronLengthMin;
    hpmap["hVertexZMin_"] = hVertexZMin;
    hpmap["hVertexSectorMin_"] = hVertexSectorMin;
    hpmap["hVertexRMin_"] = hVertexRMin;
    hpmap["hElectronFirstGgMin_"] = hElectronFirstGgMin;
    hpmap["hElectronLastGgMin_"] = hElectronLastGgMin;
    hpmap["hVertexMinDistPromptGgMin_"] = hVertexMinDistPromptGgMin;
    hpmap["hElectronLDCorrMin_"] = hElectronLDCorrMin;
    hpmap["hVertexZSecMin_"] = hVertexZSecMin;
    hpmap["hNLowEGammas_"] = hNLowEGammas;
    hpmap["hLowEGammaEnergy_"] = hLowEGammaEnergy;
    hpmap["hSummedLowEGammaE_"] = hSummedLowEGammaE;
    hpmap["hLowEMinDistPromptGg_"] = hLowEMinDistPromptGg;
    hpmap["hEeMaxVEeMin_"] = hEeMaxVEeMin;
    hpmap["hNAfterCuts_"] = hNAfterCuts;
    hpmap["hEffAfterCuts_"] = hEffAfterCuts;
    hpmap["hNAPromptGgHitsDist2VertexMin_"] = hNAPromptGgHitsDist2VertexMin;
    hpmap["hTrackSignMax_"] = hTrackSignMax;
    hpmap["hTrackSignMin_"] = hTrackSignMin;

    hpmap["hnGammaClusters_"] = hnGammaClusters;
    hpmap["hnInCluster_"] = hnInCluster;
    hpmap["hclusterHitEnergy_"] = hclusterHitEnergy;
    hpmap["hclusterHitEnergyMin_"] = hclusterHitEnergyMin;
    hpmap["hclusterHitEnergyMax_"] = hclusterHitEnergyMax;
    hpmap["hnLowEnergyHits_"] = hnLowEnergyHits;
    hpmap["hclusterEnergy_"] = hclusterEnergy;

    hpmap["hFoilSide_"] = hFoilSide;
    hpmap["hTrueVertexR_"] = hTrueVertexR;
    hpmap["hFoilSideTrueVertexR_"] = hFoilSideTrueVertexR;
    hpmap["hTrueVertexLayer_"] = hTrueVertexLayer;

    hpmap["hSingleEnergy_"] = hSingleEnergy;
    hpmap["hHighEnergy_"] = hHighEnergy;
    hpmap["hLowEnergy_"] = hLowEnergy;
    hpmap["hHighLowEnergy_"] = hHighLowEnergy;

    hpmap["hTotalClusterEnergy_"] = hTotalClusterEnergy;

    hpmap["hHotSpotLocation_"] = hHotSpotLocation;

    hpmap["hEnergySum_"] = hEnergySum;
    hpmap["hEnergyDiff_"] = hEnergyDiff;


    for(std::map<TString, TH1*>::iterator it{hpmap.begin()}; it != hpmap.end(); ++ it)
    {
        //std::cout << it->second->GetName() << std::endl;
        it->second->Sumw2();
    }

    //#if TRUTH_ENABLE
    //    if((sampleName.CompareTo("nd150_rot_2n2b_m4") == 0) ||
    //       (sampleName.CompareTo("nd150_rot_2b2n_m4") == 0))
    //    {
    //        hpmap["hMCElectronEnergy_"] = hMCElectronEnergy;
    //        hpmap["hMCTrueElectronEnergy_"] = hMCTrueElectronEnergy;
    //    }
    //#endif

    ///////////////////////////////////////////////////////////////////////////


    Long64_t total_event_count = 0;
    Long64_t event_pass_count = 0;


    // this idea didn't work - designed to stop program crashing when file
    // and path is not found
    /*
    TFile *aFile = new TFile;
    for(;;)
    {
        TString filepathname;

        if((sampleName.CompareTo("nd150_rot_2b2n_m4") == 0)   ||
           (sampleName.CompareTo("nd150_rot_2n2b_m4") == 0))
        {
            filepathname = filePath + thePath + sampleName + "/Nd150_2eNg_output_truth.root";
        }
        else
        {
            filepathname = filePath + thePath + sampleName + "/Nd150_2eNg_output.root";
        }

        //aFile = TFile::Open(filepathname);
        aFile->Open(filepathname);

        if(aFile->IsOpen())
        {
            break;
        }
        else
        {
            std::cout << "Failed to open input file " << filepathname << std::endl;
            std::cin.get();
        }
    }
    */
    TFile *aFile;
    if((sampleName.CompareTo("nd150_rot_2b2n_m4") == 0)   ||
       (sampleName.CompareTo("nd150_rot_2n2b_m4") == 0))
    {
        #if TRUTH_ENABLE
        //aFile = TFile::Open(filePath + thePath + sampleName + "/Nd150_2eNg_output_truth.root");
        aFile = TFile::Open(filePath + thePath + sampleName + "/Nd150_2eNg_output_truth_NEW_all_3.root");
        #else
        aFile = TFile::Open(filePath + thePath + sampleName + "/Nd150_2eNg_output.root");
        #endif
    }
    else
    {
        aFile = TFile::Open(filePath + thePath + sampleName + "/Nd150_2eNg_output.root");
    }
    TTree *theTree = (TTree*)aFile->Get("Nd150_2eNg/Nd150_2eNg");

    #if TRUTH_ENABLE
    TFile *outputFile = nullptr;
    TTree *outputTree = nullptr;
    TFile *outputFile_small = nullptr;
    TTree *outputTree_small = nullptr;
    if(mode_flag == 0)
    {
        outputFile_small = new TFile(filePath + thePath + sampleName + "/Nd150_2eNg_output_postprocessed_small.root", "recreate");
        outputFile_small->mkdir("Nd150_2eNg");                                                                           
        outputFile_small->cd("Nd150_2eNg");                                                  
        outputTree_small = new TTree("Nd150_2eNg", "Nd150_2eNg");
    
        #if 0
        if((sampleName.CompareTo("nd150_rot_2b2n_m4") == 0)   ||
           (sampleName.CompareTo("nd150_rot_2n2b_m4") == 0))
        {
            // create output tree with information saved
            /*
            outputFile = new TFile("Nd150_2eNg_output_truth_postprocessed.root", "recreate");
            outputFile->mkdir("Nd150_2eNg");                                                                           
            outputFile->cd("Nd150_2eNg");                                                  
            outputTree = new TTree("Nd150_2eNg", "Nd150_2eNg");

            outputFile_small = new TFile("Nd150_2eNg_output_truth_postprocessed_small.root", "recreate");
            outputFile_small->mkdir("Nd150_2eNg");                                                                           
            outputFile_small->cd("Nd150_2eNg");                                                  
            outputTree_small = new TTree("Nd150_2eNg", "Nd150_2eNg");
            */
            //#if TRUTH_ENABLE
            /*
            outputFile = new TFile(filePath + thePath + sampleName + "/Nd150_2eNg_output_truth_postprocessed.root", "recreate");
            outputFile->mkdir("Nd150_2eNg");                                                                           
            outputFile->cd("Nd150_2eNg");                                                  
            outputTree = new TTree("Nd150_2eNg", "Nd150_2eNg");
            */

            /*
            outputFile_small = new TFile(filePath + thePath + sampleName + "/Nd150_2eNg_output_truth_postprocessed_small.root", "recreate");
            outputFile_small->mkdir("Nd150_2eNg");                                                                           
            outputFile_small->cd("Nd150_2eNg");                                                  
            outputTree_small = new TTree("Nd150_2eNg", "Nd150_2eNg");
            */
            //#else
            /*
            outputFile = new TFile(filePath + thePath + sampleName + "/Nd150_2eNg_output_postprocessed.root", "recreate");
            outputFile->mkdir("Nd150_2eNg");                                                                           
            outputFile->cd("Nd150_2eNg");                                                  
            outputTree = new TTree("Nd150_2eNg", "Nd150_2eNg");
            */
            // NOTE: changed 2020-08-28: all files uniform naming scheme
            outputFile_small = new TFile(filePath + thePath + sampleName + "/Nd150_2eNg_output_postprocessed_small.root", "recreate");
            outputFile_small->mkdir("Nd150_2eNg");                                                                           
            outputFile_small->cd("Nd150_2eNg");                                                  
            outputTree_small = new TTree("Nd150_2eNg", "Nd150_2eNg");
            //#endif
        }
        else
        {
            /*
            outputFile = new TFile(filePath + thePath + sampleName + "/Nd150_2eNg_output_postprocessed.root", "recreate");
            outputFile->mkdir("Nd150_2eNg");                                                                           
            outputFile->cd("Nd150_2eNg");                                                  
            outputTree = new TTree("Nd150_2eNg", "Nd150_2eNg");
            */
            outputFile_small = new TFile(filePath + thePath + sampleName + "/Nd150_2eNg_output_postprocessed_small.root", "recreate");
            outputFile_small->mkdir("Nd150_2eNg");                                                                           
            outputFile_small->cd("Nd150_2eNg");                                                  
            outputTree_small = new TTree("Nd150_2eNg", "Nd150_2eNg");
        }
        #endif
    }
    #endif

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

    #if TRUTH_ENABLE
    double      trueElectronEnergy[2]; // TODO
    #endif


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
   
    #if TRUTH_ENABLE
    if(mode_flag == 0)
    {
        if((sampleName.CompareTo("nd150_rot_2b2n_m4") == 0)   ||
           (sampleName.CompareTo("nd150_rot_2n2b_m4") == 0))
        {
            theTree->SetBranchAddress("trueElectronEnergy" , trueElectronEnergy);
        }
        else
        {
            trueElectronEnergy[0] = 0.;
            trueElectronEnergy[1] = 0.;
        }

        ///////////////////////////////////////////////////////////////////////
        // do this for all sample types
        ///////////////////////////////////////////////////////////////////////

        // big tree, all information

        /*
        outputTree->Branch("Event"                      , &event                        , "event/I");
        outputTree->Branch("Run"                        , &run                          , "run/I");
        outputTree->Branch("runStatus"                  , &runStatus                    , "runStatus/I");
        outputTree->Branch("nElectrons"                 , &nElectrons                   , "nElectrons/I");

        outputTree->Branch("radonWeight"                , &radonWeight                  , "radonWeight/D");
        outputTree->Branch("bi210Weight"                , &bi210Weight                  , "bi210Weight/D");
        outputTree->Branch("foilSide"                   , &foilSide                     , "foilSide/I");
        outputTree->Branch("eventTime"                  , &eventTime                    , "eventTime/D");

        outputTree->Branch("trueVertexR"                , &trueVertexR                  , "trueVertexR/D");
        outputTree->Branch("trueVertexZ"                , &trueVertexZ                  , "trueVertexZ/D");
        outputTree->Branch("trueVertexSector"           , &trueVertexSector             , "trueVertexSector/D");
        outputTree->Branch("trueVertexLayer"            , &trueVertexLayer              , "trueVertexLayer/I");

        outputTree->Branch("electronEnergy"             , electronEnergy                , "electronEnergy[2]/D");
        outputTree->Branch("eTrackLength"               , eTrackLength                  , "eTrackLength[2]/D");
        outputTree->Branch("electronSide"               , electronSide                  , "electronSide[2]/I");
        outputTree->Branch("trackSign"                  , trackSign                     , "trackSign[2]/D");
        outputTree->Branch("electronMeasTime"           , electronMeasTime              , "electronMeasTime[2]/D");
        outputTree->Branch("electronDMeasTime"          , electronDMeasTime             , "electronDMeasTime[2]/D");
        outputTree->Branch("electronBlockType"          , electronBlockType             , "electronBlockType[2]/I");

        outputTree->Branch("internalPull"               , &eeInternalPull               , "eeInternalPull/D");
        outputTree->Branch("internalProb"               , &eeInternalProb               , "eeInternalProb/D");
        outputTree->Branch("externalPull"               , &eeExternalPull               , "eeExternalPull/D");
        outputTree->Branch("externalProb"               , &eeExternalProb               , "eeExternalProb/D");
        outputTree->Branch("cosee"                      , &cosee                        , "cosee/D");
        outputTree->Branch("cosee_weight"               , &cosee_weight                 , "cosee_weight/D");

        outputTree->Branch("electronPMT"                , electronPMT                   , "electronPMT[2]/I");
        outputTree->Branch("electronLDFlag"             , electronLDFlag                , "electronLDFlag[2]/I");
        outputTree->Branch("electronLDCorr"             , electronLDCorr                , "electronLDCorr[2]/D");
        outputTree->Branch("electronLDCorrErr"          , electronLDCorrErr             , "electronLDCorrErr[2]/D");

        outputTree->Branch("vertexZ"                    , vertexZ                       , "vertexZ[2]/D");
        outputTree->Branch("vertexSec"                  , vertexSec                     , "vertexSec[2]/D");
        outputTree->Branch("vertexR"                    , vertexR                       , "vertexR[2]/D");

        outputTree->Branch("vertexInHotSpot"            , vertexInHotSpot               , "vertexInHotSpot[2]/O");

        outputTree->Branch("firstGgHitLayer"            , firstGgHitLayer               , "firstGgHitLayer[2]/I");
        outputTree->Branch("lastGgHitLayer"             , lastGgHitLayer                , "lastGgHitLayer[2]/I");

        outputTree->Branch("NAPromptGgHits"             , &NAPromptGgHits               , "NAPromptGgHits/I");
        outputTree->Branch("NAPromptGgHitsSide"         , NAPromptGgHitsSide            , "NAPromptGgHitsSide/I"); // TODO: should be array?
        outputTree->Branch("NAPromptGgHitsDist2Vertex"  , NAPromptGgHitsDist2Vertex     , "NAPromptGgHitsDist2Vertex[NAPromptGgHits]/D");
        outputTree->Branch("NAPromptGgHitsDist2Calo"    , NAPromptGgHitsDist2Calo       , "NAPromptGgHitsDist2Calo[NAPromptGgHits]/D");

        outputTree->Branch("nGammaClusters"             , &nGammaClusters               , "nGammaClusters/I");
        outputTree->Branch("nInCluster"                 , &nInCluster                   , "nInCluster[nGammaClusters]/I");
        outputTree->Branch("clusterEnergy"              , clusterEnergy                 , "clusterEnergy[nGammaClusters]/D");
        outputTree->Branch("clusterTimeSpan"            , clusterTimeSpan               , "clusterTimeSpan[nGammaClusters]/D");

        outputTree->Branch("nTotalClusterHits"          , &nTotalClusterHits            , "nTotalClusterHits/I");
        outputTree->Branch("clusterHitEnergy"           , clusterHitEnergy              , "clusterHitEnergy[nTotalClusterHits]/D");
        outputTree->Branch("clusterHitPMT"              , clusterHitPMT                 , "clusterHitPMT[nTotalClusterHits]/I");
        outputTree->Branch("clusterHitLDFlag"           , clusterHitLDFlag              , "clusterHitLDFlag[nTotalClusterHits]/I");
        outputTree->Branch("clusterHitLDCorr"           , clusterHitLDCorr              , "clusterHitLDCorr[nTotalClusterHits]/D");
        outputTree->Branch("clusterHitLDCorrErr"        , clusterHitLDCorrErr           , "clusterHitLDCorrErr[nTotalClusterHits]/D");
        outputTree->Branch("clusterHitSec"              , clusterHitSec                 , "clusterHitSec[nTotalClusterHits]/D");
        outputTree->Branch("clusterHitZ"                , clusterHitZ                   , "clusterHitZ[nTotalClusterHits]/D");

        outputTree->Branch("trueElectronEnergy"         , trueElectronEnergy            , "trueElectronEnergy[2]/D"); 
        */

        // small tree, essential information

        outputTree_small->Branch("Run"                          , &run                      , "run/I");
        outputTree_small->Branch("eventTime"                    , &eventTime                , "eventTime/D");
        //outputTree_small->Branch("nElectrons"                 , &nElectrons               , "nElectrons/I"); // could remove
        outputTree_small->Branch("radonWeight"                  , &radonWeight              , "radonWeight/D");
        outputTree_small->Branch("bi210Weight"                  , &bi210Weight              , "bi210Weight/D");
        outputTree_small->Branch("electronEnergy"               , electronEnergy            , "electronEnergy[2]/D");
        // truth information, save in all files to maintain homogenaity
        outputTree_small->Branch("trueElectronEnergy"           , trueElectronEnergy        , "trueElectronEnergy[2]/D");
    }
    #else
    if(mode_flag == 0)
    {
        // TODO:
        // don't really care about filling this in because never run code
        // without TRUTH_ENABLE now
    }
    #endif

    if(mode_flag == 0)
    {
        std::cout << "Processing: " << sampleName << std::endl;
    }

    Int_t electron_cut_A = 0;
    //Int_t electron_cut_B = 0;

    int nCuts = 17; // 26;
    Int_t cut_counter[nCuts];
    //Int_t cut_counter_index = 0;
    for(int i = 0; i < nCuts; ++ i)
    {
        cut_counter[i] = 0;
    }

    Long_t events = (Long_t)theTree->GetEntries();
    for(Long_t event_i = 0; event_i < events; event_i++)
    {

        bool cut = false;
        //if(event_i % 1000 == 0)
        //{
        //    std::cout << "\r Processing complete : " << 100 * event_i / events << "%" << std::flush;
        //}
        theTree->GetEvent(event_i);


        if((electronEnergy[0] < 0.3) || (electronEnergy[1] < 0.3))
        {
            ++ electron_cut_A;
        }
        
        ++ total_event_count;







        double weight = 1.;
        if(sampleName.CompareTo("bi214_swire") == 0)
        {
            weight = radonWeight;
        }
        // TODO: due to secular equlibrium, the following MAY be required
        // (unsure at present time)
        else if(sampleName.CompareTo("pb214_swire") == 0)
        {
            weight = radonWeight;
        }
        else if(sampleName.Contains("bi210"))
        {
            // halflife
            // https://periodictable.com/Isotopes/083.210/index.p.full.dm.html
            const double T21 = 86400. * 5.011574074074;
            // lambda
            const double lambda = std::log(2.0) / T21;
            double arg = -lambda * eventTime;
            // I assume event time is in SI units
            weight = bi210Weight * std::exp(arg);
        }
        else if(sampleName.Contains("co60"))
        {
            // halflife
            // https://periodictable.com/Isotopes/027.60/index.p.html
            const double T21 = 31557600. * 5.274923896499;
            // lambda
            const double lambda = std::log(2.0) / T21;
            double arg = -lambda * eventTime;
            // I assume event time is in SI units
            weight = std::exp(arg);
        }
        else if(sampleName.Contains("eu152"))
        {
            // halflife
            // https://periodictable.com/Isotopes/063.152/index.p.html
            const double T21 = 31557600. * 13.537;
            // lambda
            const double lambda = std::log(2.0) / T21;
            double arg = -lambda * eventTime;
            // I assume event time is in SI units
            weight = std::exp(arg);
        }
        else if(sampleName.Contains("eu154"))
        {
            // halflife
            // https://periodictable.com/Isotopes/063.154/index.html
            const double T21 = 31557600. * 8.593;
            // lambda
            const double lambda = std::log(2.0) / T21;
            double arg = -lambda * eventTime;
            // I assume event time is in SI units
            weight = std::exp(arg);
        }
        // TODO: can fill 208Tl with 208Tl weight here instead of changing weight arbitrarily later
        // TODO: I get completely different numbers
        /*    else if ( histName.CompareTo("bi210_swire") == 0 ) {weight = bi210Weight * pow(0.5,eventTime/433123.);}
        else if ( (histName.CompareTo("bi210_swire") == 0) || (histName.CompareTo("bi210_sscin") == 0) || (histName.CompareTo("bi210_sfoil") == 0)) {weight = bi210Weight * pow(0.5,eventTime/433123.);}
        else if ((histName.CompareTo("co60_cuTower")==0) || (histName.CompareTo("co60_steel")==0) || (histName.CompareTo("co60_muMetal")==0) || (histName.CompareTo("co60_cuPetals")==0) )weight = pow(0.5,(eventTime/166344192.));
        else if ( (histName.CompareTo("eu152_sscin")==0) || (histName.CompareTo("eu152_int")==0)) {weight = pow(0.5,(eventTime/4.265e8));}
        else if ((histName.CompareTo("eu154_int")==0)) {weight = pow(0.5,(eventTime/2.71169e8));}
        */
        
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


        // count number of events before cuts
        // count the initial number of events
        int cc = 0;
        ++ cut_counter[cc]; // cut 0 - input count
        ++ cc;






        // additional sample splitting
        // side = 1: OUT
        // side = 0: IN
        // side = -1: FAIL
        // layer: 0 - 8 inclusive = valid
        // -1 = fail
        #if 0
        if(mode_flag_2 == 1)
        {
            if(sample_split_flags == 0x01)
            {
                // IN
                if(foilSide != 0)
                {
                    continue;
                }
            }
            else if(sample_split_flags == 0x02)
            {
                // OUT
                if(foilSide != 1)
                {
                    continue;
                }
            }
            else if(sample_split_flags == 0x04)
            {
                // L0 IN
                name_sample_split_additional = "_INL0";
                if(foilSide != 0)
                {
                    continue;
                }
                if(trueVertexLayer != 0)
                {
                    continue;
                }
            }
            else if(sample_split_flags == 0x08)
            {
                // L0 OUT
                name_sample_split_additional = "_OUTL0";
                if(foilSide != 1)
                {
                    continue;
                }
                if(trueVertexLayer != 0)
                {
                    continue;
                }
            }
            else if(sample_split_flags == 0x10)
            {
                // L>0 IN+OUT
                name_sample_split_additional = "_Lg0";
                if(trueVertexLayer < 1 || trueVertexLayer > 8)
                {
                    continue;
                }
            }
        }
        #endif


        ++ cut_counter[cc]; // cut 0 - input count
        ++ cc;






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


        // phase cut moved to end


        // maybe 5.87
        //if(vertexSec[0] < 5.74 || vertexSec[0] > 5.88 || vertexSec[1] < 5.74 || vertexSec[1] > 5.88)
        // The following values taken from Position Paper Section 2,
        // seems to match data from 2e directory
        // No cut in 2eNg, do not know why
        if(mode_flag_2 == 1)
        {
            //if(vertexSec[0] < 5.7371 + 0.001 || vertexSec[0] > 5.8706 - 0.001 ||
            //   vertexSec[1] < 5.7571 + 0.001 || vertexSec[1] > 5.8706 - 0.001)
            //if(vertexSec[0] < 5.7371 + 0.000 || vertexSec[0] > 5.8706 - 0.001 ||
            //   vertexSec[1] < 5.7571 + 0.000 || vertexSec[1] > 5.8706 - 0.001)
            if(vertexSec[0] < 5.7371 || vertexSec[0] > 5.8706 || vertexSec[1] < 5.7371 || vertexSec[1] > 5.8706)
            // when working with 2e, this should not cut any events
            //if(vertexSec[0] < 5.7 || vertexSec[0] > 5.9 || vertexSec[1] < 5.7 || vertexSec[1] > 5.9)
            {
                continue;
            }
        }
        ++ cut_counter[cc]; // cut 2 - vertex sector
        ++ cc;


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
       
        ++ cut_counter[cc];
        ++ cc;

    */

        // Set the high and low energy index
        int highE_index = -1;
        int lowE_index = -1;
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



        if(vertexInHotSpot[0])
        {
            hHotSpotLocation->Fill(vertexSec[0], vertexZ[0], 1.0);
        }
        if(vertexInHotSpot[1])
        {
            hHotSpotLocation->Fill(vertexSec[1], vertexZ[1], 1.0);
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
        ++ cut_counter[cc];
        ++ cc;
        
        // 4: No alpha candidates near vertex
        // (5) No alpha candidates near vertex
        // TODO: what does near mean?
        ++ cut_counter[cc];
        ++ cc
        
        // 5: E_e >= 0.2 MeV
        // (6) E_e >= 0.2 MeV
        // (23) No gamma with E_gamma > 200 keV
        if( (electronEnergy[0] < 0.2) || (electronEnergy[1] < 0.2) ) continue;
        ++ cut_counter[cc];
        ++ cc;

        // 6: <= 20 total free scintillators and <= 10 gamma clusters
        // (7) <= 20 total free scintillators and <= 10 gamma clusters (data management)
        // NOTE: nothing to be done here, this cut is done in data management
        // by file /home/blotsd/NEMO3/release_1.0.0/SDBanalysis/src/Nd150_2e.cpp
        // by file /home/blotsd/NEMO3/release_1.0.0/SDBanalysis/src/Nd150_2eNg.cpp
        if(nGammaClusters > 10) continue; 
        ++ cut_counter[cc];
        ++ cc;

        // 7: All PMTs have good status
        // (8) All PMTs have good status
        // TODO : not a clue how to implement this
        ++ cut_counter[cc];
        ++ cc;

        // 8: Isolated electron calorimeter hits
        // (9) Isolated electron calorimeter hits
        // TODO: again not a clue
        ++ cut_counter[cc];
        ++ cc;

        // 9: Electron hits in different scintillator blocks
        // (10) Electron hits in different scintillator blocks
        // [3]
        // NOTE: This is the same as cut number 27 and [3] 
        // 27: Each electron track must be associated to a unique scinitllator block
        // [3] this index appears multiple times
        if(electronPMT[0] == electronPMT[1]) continue;
        ++ cut_counter[cc];
        ++ cc

        // 10: Impacts on block faces
        // (11) Impacts on block faces
        // TODO: again not a clue
        ++ cut_counter[cc];
        ++ cc;

        // 11: All PMT have good laser survey and counting rates identified by LD and HS flags
        // (12) Electron LD flags <= 2
        // (13) Electrons have LD corrections
        // (14) Electron HS flags = 0
        // (15) Gamma hits have LD flags <= 2
        // (16) Gamma hits have LD corrections
        // (17) No gamma hits with HS = 1
        // TODO: not sure about this, this might only be 1 part of it
        if((electronLDFlag[0] > 0) || (electronLDFlag[1] > 0)) continue;
        ++ cut_counter[cc];
        ++ cc
        
        // 12: Tracks start in Gg layer 0 or 1
        // (18) Tracks start in Gg later 0 or 1
        if((firstGgHitLayer[0] > 1) || (firstGgHitLayer[1] > 1)) continue;
        ++ cut_counter[cc];
        ++ cc;
        
        // 13: Tracks miss <= 1 Gg layer before impact on calorimeter
        // (19) Tracks miss <= 1 gg layer
        // TODO: not a clue
        ++ cut_counter[cc];
        ++ cc;
        
        // 14: Each L_e > 20 cm
        // (20) Track lengths each > 20 cm
        if((eTrackLength[0] <= 20. ) || (eTrackLength[1] <= 20.)) continue;
        ++ cut_counter[cc];
        ++ cc;
        
        // 15: <= 5 unassociated prompt Gg hits
        // (21) <= 5 NA prompt Gg hits
        if(NAPromptGgHits > 5) continue;
        ++ cut_counter[cc];
        ++ cc;

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
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1)
            {
                
                //if( (electronBlockType[0] == 4) || (electronBlockType[1] == 4) ) continue; // not in petals
                //if( (electronBlockType[0] == 5) || (electronBlockType[1] == 5) ) continue; // not in petals
                //if( (electronBlockType[0] < 0) || (electronBlockType[1] < 0) ) continue; // does not do anything (makes no difference to cut count)
                //if( (electronBlockType[0] > 5) || (electronBlockType[1] > 5) ) continue; // cuts way to many events
                                                                                         // disabling cut on 4 doesn't give enough events back (cuts too many events in total)
                                                                                         // disabling cut on 5 gives too many events back (doesn't cut enough events in total)

                if( (electronBlockType[0] == 4) || (electronBlockType[1] == 4) ) continue; // not in petals
                if( (electronBlockType[0] == 5) || (electronBlockType[1] == 5) ) continue; // not in petals

                //if( (electronBlockType[0] > 2) || (electronBlockType[1] > 2) ) continue;
                
                //if( (electronBlockType[0] == 3) || (electronBlockType[1] == 3) ) continue;
                //if( (electronBlockType[0] == 4) || (electronBlockType[1] == 4) ) continue; // nope.

            }
        }
        ++ cut_counter[cc]; // cut 3 - block type
        ++ cc;
        
        // 17: No gamma rays in event
        //if(nGammaClusters > 0) continue;
        // Think this is the correct way to cut any gamma rays
        // Summer removed this for Position paper
        // replacing it with no gamma rays with energy > 200 keV
        // check this TODO
        //if(mode_flag == 0)
        //{
        //    if(mode_flag_2 == 1)
        //    {
        //        int cluster_hit_index = 0;
        //        int cesum = 0;
        //        for(int i{0}; i < nGammaClusters; ++ i)
        //        {
        //            /*
        //            if(clusterEnergy[i] > 0.2)
        //            {
        //                ++ cesum;
        //            }
        //            */
        //            for(int j{0}; j < nInCluster[i]; ++ j)
        //            {
        //                if(clusterHitEnergy[cluster_hit_index] > 0.2)
        //                {
        //                    ++ cesum;
        //                }
        //                ++ cluster_hit_index;
        //            }
        //        }
        //        if(cesum > 0) continue;
        //    }
        //}
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1)
            {
                double cluster_total_energy = 0.;
                int cluster_hit_index = 0;
                int cesum = 0;

                //if(nGammaClusters > 0) continue;

                for(int i{0}; i < nGammaClusters; ++ i)
                {
                    /*
                    if(clusterEnergy[i] > 0.2)
                    {
                        ++ cesum;
                    }
                    */
                    /*
                    for(int j{0}; j < nInCluster[i]; ++ j)
                    {
                        if(clusterHitEnergy[cluster_hit_index] > 0.2)
                        {
                            ++ cesum;
                            cluster_total_energy += clusterHitEnergy[cluster_hit_index];
                        }
                        ++ cluster_hit_index;
                    }
                    */
                    //if(clusterEnergy[i] > 0.2) continue;
                    cluster_total_energy += clusterEnergy[i];
                }
                //if(cesum > 0) continue;
                //if(cluster_total_energy > 0.2) continue;
                if(cluster_total_energy > 0.0) continue;
                
            }
        }
        ++ cut_counter[cc]; // cut 4 - gamma energy
        ++ cc;
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
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1 || mode_flag_2 == 0)
            {
                if((trackSign[0] >= 0) || (trackSign[1] >= 0)) continue;
                // 2020-05-03 changed TODO change back?
                // result: makes NO DIFFERENCE ?!
                //if((trackSign[0] > 0) || (trackSign[1] > 0)) continue;
            }
        }
        ++ cut_counter[cc]; // cut 5 - track negative curvature
        ++ cc;
        
        // 19: Each L_e >= 30 cm
        // (25) Track lengths > 30 cm
        // 4. - position paper number/index, this index appears multiple times
        // 4. Each of the electron tracks must be at least 30 cm long, [have a hit in the first Gg cell layer,
        // and be associated to unique scintillator blocks.]
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1)
            {
                if((eTrackLength[0] < 30. ) || (eTrackLength[1] < 30.)) continue;
            }
        }
        ++ cut_counter[cc]; // cut 6 - track length
        ++ cc;

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
        //if((firstGgHitLayer[0] > 0) || (firstGgHitLayer[1] > 0)) continue;
        //++ cut_counter[cc]; // cut 7
        //++ cc;
        // TODO unique scintillator blocks
        // from below
        // 23: Both tracks have hit in Gg L0
        // (not 4, 18) - position paper number/index, this index appears multiple times
        // (30) Both tracks have hit in Gg L0
        //if(electronPMT[0] == electronPMT[1]) continue;
        //++ cut_counter[cc]; // cut 8
        //++ cc;
        // 27: Each electron track must be associated to a unique scinitllator block
        // (3) - position paper number/index, this index appears multiple times

        
        // 20: Tracks not from hot spots
        // (26) Tracks not from hot spots
        // 3. Electron tracks with vertices located in hot spot regions are not accepted in order to
        // improve the signal to background ratio in this channel.
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1)
            {
                //if(vertexInHotSpot[0] || vertexInHotSpot[1]) continue;

                double vertexZ_1 = vertexZ[0];
                double vertexZ_2 = vertexZ[1];
                double vertexSec_1 = vertexSec[0];
                double vertexSec_2 = vertexSec[1];

                bool flag_continue = false;

                for(int elix = 0; elix < N_HOTSPOT_ELIPSE; ++ elix)
                {

                    double x1 = vertexSec_1;
                    double y1 = vertexZ_1;
                    double h1 = hotspot_elipse[elix][0];
                    double k1 = hotspot_elipse[elix][1];
                    double a1 = hotspot_elipse[elix][2];
                    double b1 = hotspot_elipse[elix][3];

                    double p1 = (std::pow(x1 - h1, 2.0) / pow(a1, 2.0)) +
                               (std::pow(y1 - k1, 2.0) / pow(b1, 2.0));

                    if(p1 <= 1.0)
                    {
                        flag_continue = true;
                        break;
                    }

                    double x2 = vertexSec_2;
                    double y2 = vertexZ_2;
                    double h2 = hotspot_elipse[elix][0];
                    double k2 = hotspot_elipse[elix][1];
                    double a2 = hotspot_elipse[elix][2];
                    double b2 = hotspot_elipse[elix][3];

                    double p2 = (std::pow(x2 - h2, 2.0) / pow(a2, 2.0)) +
                               (std::pow(y2 - k2, 2.0) / pow(b2, 2.0));

                    if(p2 <= 1.0)
                    {
                        flag_continue = true;
                        break;
                    }
                               
                }

                if(flag_continue == true) continue;

            }
        }
        ++ cut_counter[cc]; // cut 7 (was cut 9) - track hot spot
        ++ cc;
        
        // 10. No unassociated Geiger hits are allowed within 15 cm of the calorimeter hits as this is
        // indicative of a reflected electron, which therefore has not deposited its full energy in the
        // calorimeter.
        // (27) No NA prompt Gg hits within 15 cm of scintillator hits
        // This does not appear to be included in the Table from Summers thesis
        // TODO
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1)
            {
                int ssum = 0;
                for(int i{0}; i < NAPromptGgHits; ++ i)
                {
                    if(NAPromptGgHitsDist2Calo[i] < 15.)
                    {
                        ++ ssum;
                    }
                }
                if(ssum > 0) continue;
            }
        }
        ++ cut_counter[cc]; // cut 8 - Gg hit scint dist
        ++ cc;

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
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1)
            {
                int ssum = 0;
                for(int i{0}; i < NAPromptGgHits; ++ i)
                {
                    // assuming units of cm TODO
                    if(NAPromptGgHitsDist2Vertex[i] < 15.)
                    {
                        ++ ssum;
                    }
                }
                if(ssum > 1) continue;
            }
        }
        ++ cut_counter[cc];
        ++ cc; // cut 9 (was cut 10) - Gg hit vertex dist
        
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
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1)
            {
                int ssum = 0;
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
            }
        }
        ++ cut_counter[cc]; // cut 10 (was cut 11) - Gg hit opposite side
        ++ cc;


        
        
        // 23: Both tracks have hit in Gg L0
        // TODO: copied from above, but this looks like something which should
        // go here - might need to check against layer 1 and 11 ?
        //if((electronFirstGG[0] > 1) || (electronFirstGG[1] > 1)) cut = true; TODO: changed from this?
        // (not 4, 18) - position paper number/index, this index appears multiple times
        // (30) Both tracks have hit in Gg L0
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1)
            {
                if((firstGgHitLayer[0] > 0) || (firstGgHitLayer[1] > 0)) continue;
            }
        }
        ++ cut_counter[cc]; // cut 11 - Gg L0
        ++ cc;
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
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1 || mode_flag_2 == 0)
            {
                if(eeInternalProb < 0.01) continue;
                if(eeExternalProb > 0.01) continue;
            }
        }
        ++ cut_counter[cc]; // cut 12 (was cut 13) - prob
        ++ cc;
        

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
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1)
            {
                if(fabs(vertexR[0] - vertexR[1]) > 4.) continue;
                if(fabs(vertexZ[0] - vertexZ[1]) > 8.) continue;
            }
        }
        ++ cut_counter[cc]; // cut 13 (was cut 14) - delta
        ++ cc;
        
        // save data to file if electronEnergy >= 0.2 MeV
        // and phase is P1 or P2
        if(mode_flag == 0)
        {
            if((electronEnergy[0] < 0.2) || (electronEnergy[1] < 0.2))
            {
                // do not save to file
            }
            else
            {
                if(thePhase == 0)
                {
                    //if((1869 <= run) && (run <= 3395)) cut = false;
                    
                    #if TRUTH_ENABLE
                    //outputTree->Fill();
                    outputTree_small->Fill();
                    #endif
                }
                else if(thePhase == 1)
                {
                    //if((3396 <= run) && (run <= 9186)) cut = false;

                    #if TRUTH_ENABLE
                    //outputTree->Fill();
                    outputTree_small->Fill();
                    #endif
                }
            }
        }



        // 26: Each E_e >= 300 keV
        // (33) Each E_e >= 300 keV
        // 7. The energy of each electron is required to be at least 0.3 MeV, as single electron energies
        // below this threshold are not well described by the simulation.
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1)
            {
                if((electronEnergy[0] < 0.3) || (electronEnergy[1] < 0.3)) continue;
            }
        }
        ++ cut_counter[cc]; // cut 14 (was cut 15) - energy
        ++ cc;

        // 27: Each electron track must be associated to a unique scinitllator block
        // (3) - position paper number/index, this index appears multiple times
        //if(electronPMT[0] == electronPMT[1]) continue;
        //++ cut_counter[cc];
        //++ cc;
        // moved above
        

        // Phase cut moved to be final cut
        // 1869 - 3395 is P1
        // 3396 - 9186 is P2
        if(mode_flag == 0)
        {
            cut = true;
            if(thePhase == 0)
            {
                if((1869 <= run) && (run <= 3395)) cut = false;
            }
            else if(thePhase == 1)
            {
                if((3396 <= run) && (run <= 9186)) cut = false;
            }
            else if(nocutonphase == 1)
            {
                cut = false;
            }
            else
            {
                cut = true;
            }

            // accept P1 and P2 for testing purposes
            //cut = false;
            if(cut == true) continue;
        }
        ++ cut_counter[cc]; // cut 1 - phase
        ++ cc;


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

        ++ event_pass_count;

        // moved back here from above
        hRun->Fill(run, weight);
        hNElectrons->Fill(nElectrons, weight);

        //#if TRUTH_ENABLE
        //if((sampleName.CompareTo("nd150_rot_2b2n_m4") == 0)   ||
        //   (sampleName.CompareTo("nd150_rot_2n2b_m4") == 0))
        //{
        //    hMCElectronEnergy->Fill(electronEnergy[0], electronEnergy[1]);
        //    hMCTrueElectronEnergy->Fill(trueElectronEnergy[0], trueElectronEnergy[1]);
        //}
        //#endif

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

        hTrackSignMax->Fill(trackSign[highE_index], weight);
        hTrackSignMin->Fill(trackSign[lowE_index], weight);


        hnGammaClusters->Fill(nGammaClusters, weight);

        hFoilSide->Fill(foilSide, weight);
        hTrueVertexR->Fill(trueVertexR, weight);
        hFoilSideTrueVertexR->Fill(foilSide, trueVertexR, weight);
        hTrueVertexLayer->Fill(trueVertexLayer, weight);
        
        hSingleEnergy->Fill(electronEnergy[lowE_index], weight);
        hSingleEnergy->Fill(electronEnergy[highE_index], weight);
        hHighEnergy->Fill(electronEnergy[highE_index], weight);
        hLowEnergy->Fill(electronEnergy[lowE_index], weight);
        hHighLowEnergy->Fill(electronEnergy[lowE_index], electronEnergy[highE_index], weight);
        hEnergySum->Fill(electronEnergy[highE_index] + electronEnergy[lowE_index], weight);
        hEnergyDiff->Fill(electronEnergy[highE_index] - electronEnergy[lowE_index], weight);


        int hit_counter = 0;
        int nLowEnergyHits = 0;
        double clusterHitEnergyMin;
        double clusterHitEnergyMax;
        double total_cluster_energy;
        for(int i = 0; i < nGammaClusters; ++ i)
        {
            hnInCluster->Fill(nInCluster[i], weight);
            hclusterEnergy->Fill(clusterEnergy[i], weight);

            clusterHitEnergyMin = 0.;
            clusterHitEnergyMax = 0.;
            total_cluster_energy = 0.0;
            for(int j = 0; j < nInCluster[i]; ++ j)
            {
                total_cluster_energy += clusterHitEnergy[hit_counter];

                if(j == 0)
                {
                    clusterHitEnergyMin = clusterHitEnergy[hit_counter];
                    clusterHitEnergyMax = clusterHitEnergy[hit_counter];
                }

                hclusterHitEnergy->Fill(clusterHitEnergy[hit_counter], weight);
                ++ hit_counter;

                if(clusterHitEnergy[hit_counter] < .2)
                {
                    ++ nLowEnergyHits;
                }
            }
            hTotalClusterEnergy->Fill(total_cluster_energy, weight);
            if(nInCluster[i] >= 2)
            {
                // don't store case where no data was set (0)
                // or case where Min = Max (1)
                hclusterHitEnergyMin->Fill(clusterHitEnergyMin, weight);
                hclusterHitEnergyMax->Fill(clusterHitEnergyMax, weight);
            }
        }
        hnLowEnergyHits->Fill(nLowEnergyHits, weight);
        if(hit_counter != nTotalClusterHits)
        {
            std::cout << "hit_counter=" << hit_counter << ", nTotalClusterHits=" << nTotalClusterHits << std::endl;
            std::cin.get();
        }



    } //~event_i
     

    if(mode_flag == 0)
    {
        //std::cout << hSingleEnergy->GetName() << " -> " << hSingleEnergy->GetBinContent(10) << std::endl;
    }


    if(mode_flag == 0)
    {
        ofile_cutcount  << sampleName << "    " << "" << total_event_count << "    " << "" << event_pass_count << std::endl;
        std::cout       << sampleName << "    " << "total_event_count=" << total_event_count << "    " << "event_pass_count=" << event_pass_count << std::endl;
    }
    //std::cout << "electron_cut_A=" << electron_cut_A << std::endl;
    std::cout << std::endl;

    #if TRUTH_ENABLE
    if(mode_flag == 0)
    {
        //if((sampleName.CompareTo("nd150_rot_2b2n_m4") == 0)   ||
        //   (sampleName.CompareTo("nd150_rot_2n2b_m4") == 0))
        //{
        //outputTree->Write();
        //outputFile->Close();

        outputTree_small->Write();
        outputFile_small->Close();
        //}
    }
    #endif

    aFile->Close();
    
    std::string cut_description[nCuts];
    int cc = 0;
    cut_description[cc] = "Input count";
    ++ cc;
    cut_description[cc] = "IN/OUT/Layer Sample Split";
    ++ cc;
    cut_description[cc] = "Vertex sector";
    ++ cc;
    cut_description[cc] = "Block type";
    ++ cc;
    cut_description[cc] = "No gamma with E > 0.2 MeV";
    ++ cc;
    cut_description[cc] = "Tracks negative curvature";
    ++ cc;
    cut_description[cc] = "Track length > 30 cm";
    ++ cc;
    cut_description[cc] = "Tracks not from hot spots";
    ++ cc;
    cut_description[cc] = "0 Prompt Gg hit within 15 cm scint";
    ++ cc;
    cut_description[cc] = "1 Prompt Gg hit within 15 cm vertex";
    ++ cc;
    cut_description[cc] = "Unassoc prompt Gg on opposite side of foil";
    ++ cc;
    cut_description[cc] = "Tracks hits in Gg L0";
    ++ cc;
    cut_description[cc] = "Pint / Pext";
    ++ cc;
    cut_description[cc] = "Delta R / Delta Z";
    ++ cc;
    cut_description[cc] = "Energy 300 keV";
    ++ cc;
    cut_description[cc] = "Phase";
    ++ cc;
    cut_description[cc] = "none";

    //std::cout << "Here are the cut counts" << std::endl;
    std::string of_cutcount_filename;
    if(mode_flag == 0)
    {
        of_cutcount_filename = "of_cutcount_processeddata_" + std::string(sampleName) + ".txt";
    }
    else if(mode_flag == 1)
    {
        of_cutcount_filename = "of_cutcount_rawdata_" + std::string(sampleName) + ".txt";
    }
    std::ofstream of_cutcount(of_cutcount_filename, std::ios::out);
    of_cutcount << sampleName << std::endl;
    for(int i = 0; i < nCuts; i++)
    {
        //std::cout << "Passed cut " << i << ": " << cut_counter[i] << " events. Description=\"" << cut_description[i] << "\"" << std::endl;
        of_cutcount << "Passed cut " << i << ": " << cut_counter[i] << " events. Description=\"" << cut_description[i] << "\"" << std::endl; 
    }
    of_cutcount << std::endl;
    

    for(int i = 0; i < nCuts; i++ )
    {
        // NOTE: bug, should be i + 1, bin 0 is underflow, bin -1 does not exist
        //hNAfterCuts->SetBinContent(i - 1, cut_counter[i]);
        hNAfterCuts->SetBinContent(i + 1, cut_counter[i]);
    }


    TFile *myFile = new TFile("Nd150_2e_P" + Phase + ".root", "UPDATE");
    TDirectory *dir = gDirectory;
    //dir->cd("singleHistos/unscaled");
    //dir->cd("rawdata");
    //dir->mkdir(sampleName);
    //dir->cd(sampleName);
    // TODO: don't want dirs organized by sample name, want them organized
    // by histogram name

    //dir->cd("rawdata/" + sampleName);

    //dir->pwd();
    //dir->print();

    //std::cout << "writing histograms into \"singleHistos/unscaled\'" << std::endl;
    //std::cout << "writing: singleHistos/unscaled/" << hTotalE->GetName() << std::endl;
    //std::cout << "writing histograms into \"rawdata" << "\" or \"processeddata\"" << std::endl;
    //std::cout << "writing: rawdata or processeddata ... " << hTotalE->GetName() << std::endl;

    //for(int i = 0; i < numHistograms; ++ i)
    for(std::map<TString, TH1*>::iterator it{hpmap.begin()}; it != hpmap.end(); ++ it)
    {
        //std::cout << "i=" << i << std::endl;

        //std::string histogram_name = std::string(histogramNames[i].Data());
        //TString histogram_name = histogramNames[i];
        TString histogram_name = it->first;
        //std::cout << "histogram_name=" << histogram_name << std::endl;
        
        // only create the directories once
        // this functions called with DataFile first
        // bit of a dodgy hack
        if(sampleName.CompareTo(DataFile) == 0)
        {
            if(mode_flag == 0)
            {
                TDirectory *dir_histogram_name = dir->mkdir("processeddata/" + histogram_name);
                TDirectory *dir_histogram_name_2 = dir->mkdir("scaled/" + histogram_name);
                TDirectory *dir_histogram_name_3 = dir->mkdir("scaled2/" + histogram_name);
            }
            else if(mode_flag == 1)
            {
                TDirectory *dir_histogram_name = dir->mkdir("rawdata/" + histogram_name);
            }
        }

        if(mode_flag == 0)
        {
            myFile->cd("processeddata/" + histogram_name);
            it->second->Write();
        }
        else if(mode_flag == 1)
        {
            myFile->cd("rawdata/" + histogram_name);
            it->second->Write();
        }
        //myFile->pwd();
        
        //hpmap.at(histogram_name)->Write();
        // same as
        
    }


    myFile->Close();

}



void scale(TFile* myFile,                       // INPUT: unscaled histograms are loaded from this root file
           const std::string& activityfile,     // INPUT: activity values read from this file
           const std::string& typedir,          // INPUT: forms part of path string,
           // eg: "nd150/", "internals/", "externals/", "neighbours/",
           // all radon files are in "externals/"
           const int nSamples,                     // INPUT: number of objects in dataFiles/dataNames
           const TString *const sampleFiles,      // INPUT: list of samples, eg: bi214_sfoil_rot
           const TString *const sampleNames,      // INPUT: corresponding human readable names, Latex format, eg: ^{214}Bi Foil Surface
           const Color_t *const sampleColors,     // INPUT: sample color set using values in this array
           const double *const AcceptedTime,    // INPUT: accepted time for each phase is passed in here 
           double *const sampleNGenMC,            // OUTPUT: number of generated MC events are written to this location
           double *const sampleActivity,          // OUTPUT: activity values written to this location
           double *const sampleEfficiency,        // OUTPUT: efficiency values written to this location
           TObjArray** allSamples,              // OUTPUT: each histogram has colors set and is scaled (using scaling with activity), before being stored as an output here
           const Int_t mode_flag = 0)
{
    std::cout << "scale" << std::endl;

    std::ofstream oftemp("oftemp1.txt", std::ofstream::out | std::ofstream::app);




    TH1D *tmpHist; //we'll be using this later.

    //Let's get the internals/externals/data and add them to an array
    //TObjArray *allSamples[nSamples];
    for(int i = 0; i < nSamples; i++)
    {

        // create new object array
        //std::cout << sampleFiles[i] << std::endl;
        allSamples[i] = new TObjArray();


        // read efficiency from file
        std::ifstream ifile_cutcount("cutcount.txt");
        if(!ifile_cutcount.is_open())
        {
            std::cout << "failed to open file " << "cutcount.txt" << std::endl;
        }

        //double sample_total_event_count;
        double sample_pass_count;
        double sample_efficiency;
        while(!ifile_cutcount.eof())
        {
            std::string sample_name;
            Long64_t total_events;
            Long64_t pass_count;

            ifile_cutcount >> sample_name >> total_events >> pass_count;

            // TODO: removed std::string conversion?
            if(sample_name == std::string(sampleFiles[i]))
            {
                //sample_total_event_count = (double)total_events;
                sample_pass_count = (double)pass_count;
                //sample_efficiency = (double)pass_count / (double)total_events;
            }
        }


        // read in the information needed for calculating efficiencies
        // read in number of generated MC events
        std::ifstream inFile;
        std::string sample_file_name = sampleFiles[i].Data();
        static const std::vector<std::string> search_objects = {"_Lg0", "_INL0", "_OUTL0", "_IN", "_OUT"};
        for(int i = 0; i < search_objects.size(); ++ i)
        {
            std::string search_object = search_objects.at(i);
            if(sample_file_name.find(search_object) != std::string::npos)
            {
                //std::cout << "sample_file_name=" << sample_file_name << std::endl;
                sample_file_name = sample_file_name.substr(0, sample_file_name.size() - search_object.size());
                //std::cout << "sample_file_name=" << sample_file_name << std::endl;
                break;
            }
        }
        inFile.open(filePath + typedir + sample_file_name + "/JobSummary.txt");
        if(!inFile.is_open())
        {
            std::cout << "failed to open " << filePath + typedir + sample_file_name + "/JobSummary.txt" << std::endl;
        }
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
        if(!inFile2.is_open())
        {
            std::cout << "failed to open " << activityfile << std::endl;
        }

        bool found = false;
        while((!found) && (!inFile2.eof()))
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
                
                if(mcname.CompareTo(sampleFiles[i]) == 0)
                {
                    sampleActivity[0 * nSamples + i] = activityPhase1;
                    sampleActivity[1 * nSamples + i] = activityPhase2;
                    found = true;
                }

            }
        }
        if(!found)
        {
            std::cout << "Fault: sampleFiles[" << i << "]=" << sampleFiles[i] << std::endl;
            std::cin.get();
        }

        inFile2.close();

        for(int j = 0; j < numHistograms; j++)
        {

            /*
            TString name_append;
            if(mode_flag == 0)
            {
                //name_append = "";
            }
            else if(mode_flag == 1)
            {
                name_append = "_raw";
            }
            */
            // don't need

            TString histogram_name = histogramNames[j];
            TString new_histogram_name = histogram_name + sampleFiles[i] + "_fit";
            //TString name = "rawdata/" + histogram_name + "/" + histogram_name + sampleFiles[i];

            TString name;
            if(mode_flag == 0)
            {
                /*TString*/ name = "processeddata/" + histogram_name + "/" + histogram_name + sampleFiles[i]; // + name_append
                if(histogram_name.CompareTo("hTotalE_") == 0)
                {
                    std::cout << "Get() : " << name << " from file, Clone() : " << new_histogram_name << std::endl;
                }
                //tmpHist = (TH1D*)myFile->Get(name)->Clone(new_histogram_name);
            }
            else if(mode_flag == 1)
            {
                /*TString*/ name = "rawdata/" + histogram_name + "/" + histogram_name + sampleFiles[i] + "_raw"; // + name_append;
                if(histogram_name.CompareTo("hTotalE_") == 0)
                {
                    std::cout << "Get() : " << name << " from file, Clone() : " << new_histogram_name << std::endl;
                }
                //tmpHist = (TH1D*)myFile->Get(name)->Clone(new_histogram_name);
            }
            try
            {
                tmpHist = (TH1D*)myFile->Get(name)->Clone(new_histogram_name);
            }
            catch(std::exception &e)
            {
                std::cout << "name=" << name << std::endl;
                std::cout << e.what() << std::endl;
                throw e;
            }
            // NOTE: will have 2 new histograms here
            // 1: scaled, processeddata histogram (save this one to file)
            // 2: scaled, rawdata histogram (do not save this one to file)

            // this code was replaced with above: 2020-05-03
            //std::string name("singleHistos/unscaled/" + histogramNames[j] + sampleFiles[i]);
            //std::cout << "Get() : " << name << " from file, Clone() : " << histogramNames[j] + sampleFiles[i] + "_fit" << std::endl;
            //tmpHist = (TH1D*)myFile->Get(name.c_str())->Clone(histogramNames[j] + sampleFiles[i] + "_fit");
            
            //calculate efficiency
            // TODO: this is calculated using the 0th histogram, should perhaps
            // change to the hTotalE ?
            /*
            if(j == 0)
            {
                double NPass = tmpHist->Integral();
                // NPass is phase specific
                // sampleNGenMC[i] is for both phases
                std::cout << "NPass: " << NPass << std::endl;
                std::cout << "sampleNGenMC[i=" << i << "]=" << sampleNGenMC[i] << std::endl;

                sample_efficiency = sample_pass_count / sampleNGenMC[i];
                
                // TODO: don't fully understand this
                // efficiency = (number of events which pass cuts * total time) / 
                // (number of generated events * accepted time)
                // accepted time is like total time - dead time?
                //double eff = NPass / (double)(sampleNGenMC[i] * (AcceptedTime[thePhase]/(double)TotalTime) );
                // this efficiency is phase specific
                double eff = (double)NPass / (double)sampleNGenMC[i];
                //sampleEfficiency[thePhase][i] = eff;
            }
            */
            
            //tmpHist->Sumw2();
            tmpHist->SetFillColor(sampleColors[i]);
            tmpHist->SetLineColor(sampleColors[i]);
            tmpHist->SetLineWidth(0);
            tmpHist->SetTitle(sampleNames[i]);

            //tmpHist->Scale(sampleActivity[thePhase][i] * AcceptedTime[thePhase] * 0.001 / (sampleNGenMC[i] * (AcceptedTime[thePhase]/(double)TotalTime)) );
            //tmpHist->Scale(sampleActivity[thePhase * nSamples + i] * 0.001 * (double)TotalTime / sampleNGenMC[i] );
            //std::cout << "scaling, index=" << thePhase + nSamples + i << std::endl;
            
            // In Position Paper
            // Section 3.4.1, Equation 3.7
            // A = N / (t * eps)
            // A: activity, N: number of events, t: time, eps: efficiency
            // Can calculate number of events backwards from activity, efficiency and time

            
            

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

            
            // moved from stackFunction()
            // 2020-03-30: 
            //(TH1D*)allSamples[j]->At(i)->Write();
            // saving the histogram without scaling it by the activity will
            // cause the fit output (and input) to be in units of activty
            // is this actually the case, what about the factor of the detector
            // livetime which I have left in above?
            /*
            Double_t scale_factor_without_activity = (double)TotalTime / sampleNGenMC[i]; //scale_factor / the_activity;
            tmpHist->Scale(scale_factor_without_activity);
            std::cout << "writing: " << tmpHist->GetName() << std::endl;
            tmpHist->Write();

            tmpHist->Scale(scale_factor / scale_factor_without_activity);
            */



            if(TString(tmpHist->GetName()).Contains("hTotalE"))
            {
                sampleEfficiency[thePhase * nSamples + i] = tmpHist->Integral() / sampleNGenMC[i];
                std::ofstream of_sampleEfficiency("of_sampleEfficiency.txt", std::ios::app);
                of_sampleEfficiency << tmpHist->GetName() << " " << tmpHist->Integral() << " " << sampleNGenMC[i] << " " << sampleEfficiency[0 * nSamples + i] << " " << sampleEfficiency[1 * nSamples + i] << std::endl;
            }

            if(TString(tmpHist->GetName()).Contains("hTotalE"))
            {
                std::cout << "TotalTime=" << TotalTime
                          << ", sampleNGenMC[" << i << "]=" << sampleNGenMC[i]
                          << ", sampleActivity[thePhase="
                              << thePhase
                              << " * nSamples="
                              << nSamples
                              << " + i="
                              << i << "]="
                          << sampleActivity[thePhase * nSamples + i] << std::endl;
                //std::cin.get();
                std::cout << std::endl;
            }



            tmpHist->Scale(TotalTime / sampleNGenMC[i]);
            if(mode_flag == 0)
            {
                if(TString(tmpHist->GetName()).Contains("hSingleEnergy_"))
                {
                    std::cout << tmpHist->GetName() << " -> " << tmpHist->GetBinContent(10) << std::endl;
                    //std::cin.get();
                }
                //std::cout << tmpHist->GetName() << " -> scale=" << TotalTime / sampleNGenMC[i] << std::endl;
            }
            
            if(mode_flag == 0)
            {
                // write into root directory to be read back by newLogLikFitter.C
                myFile->cd("/");
                tmpHist->Write();

                // write into scaled / histogram_name subdirectory
                myFile->cd("scaled/" + histogram_name);
                TString clone_name = TString(tmpHist->GetName()) + TString("_scaled");
                TH1D *tmpHistClone = (TH1D*)tmpHist->Clone(clone_name);
                tmpHistClone->Write();
            }
            //if(sampleActivity[thePhase * nSamples + i] <= 0.0)
            //{
            //    std::cout << "fault" << std::endl;
            //    std::cout << "the name is: " << tmpHist->GetName() << std::endl;
            //    std::cout << sampleFiles[i] << std::endl;
            //    std::cin.get();
            //    //throw "fault";
            //}
            tmpHist->Scale(sampleActivity[thePhase * nSamples + i]);
            if(mode_flag == 0)
            {
                if(TString(tmpHist->GetName()).Contains("hSingleEnergy_"))
                {
                    std::cout << "after scaling by activity A=" << sampleActivity[thePhase * nSamples + i] << " the 10th bin contains:" << std::endl;
                    std::cout << tmpHist->GetName() << " -> " << tmpHist->GetBinContent(10) << std::endl;
                    //std::cin.get();
                }
                //std::cout << tmpHist->GetName() << " -> scale=" << TotalTime / sampleNGenMC[i] << std::endl;
            }
            if(mode_flag == 0)
            {
                // write into scaled2 / histogram_name subdirectory
                myFile->cd("scaled2/" + histogram_name);
                TString clone_name = TString(tmpHist->GetName()) + TString("_scaled2");
                TH1D *tmpHistClone = (TH1D*)tmpHist->Clone(clone_name);
                tmpHistClone->Write();
            }

            std::ofstream of_numberofevents("of_numberofevents.txt", std::ofstream::out | std::ofstream::app);
            if(mode_flag == 0 && TString(tmpHist->GetName()).Contains("hTotalE"))
            {
                std::cout << tmpHist->GetName() << " number of events " << tmpHist->Integral() << std::endl;
                of_numberofevents << tmpHist->GetName() << " number of events " << tmpHist->Integral() << std::endl;
            }
            of_numberofevents.close();

            // NOTE: done AFTER call to Write()
            // so sample is scaled with activity value if used again from allSamples
            allSamples[i]->Add(tmpHist);

        } // all histos
    }   //all of sample type (internal/external/...)

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
                   TH1D **hAllSamples,              // OUTPUT: single histogram stack composed of samples components added together by function call to TH1::Add()
                   const TString *const sampleNames,
                   const TString *const sampleFiles,
                   TH1D **hAllMC,
                   bool allMC_first,
                   bool hAllSamples_addtoexisting,  // INPUT: flag, if true objects in allSamples are added to hAllSamples without call to Clone to change name
                   Color_t majorStackColor,
                   bool hMajorStacks_add)
{

    // for printing statements
    TString h_name;
    try
    {
        h_name = ((TH1D*)allSamples[0]->At(i))->GetName();
    }
    catch(std::exception& e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "h_name is not set" << std::endl;
    }

    // i = histogram index (type, eg; electron energy min/max, total energy, track length...)
    // j = background/data/mc index

    TH1D *tmpHist; //we'll be using this later.

    // this block relating to hAllSamples
    // which is one of hAllExternals, hAllInternals, etc...
    // add all samples of a type into a histogram for the group hAllSamples
    for(int j = 0; j < nSamples; j++)
    {
        if(h_name.Contains("hTotalE"))
        {
            std::cout << "stackfunction: j=" << j << std::endl;
        }

        if(hAllSamples_addtoexisting == false)
        {
            if(j == 0)
            {
                // NOTE: allSamples is either allExternals, allInternals, etc...
                hAllSamples[i] = (TH1D*)allSamples[j]->At(i)->Clone(sampleName);
                if(h_name.Contains("hTotalE"))
                {
                    std::cout << "adding to hAllSamples[" << i << "]: " << ((TH1D*)allSamples[j]->At(i))->GetName() << " (clone)" << std::endl;
                }
            }
            else
            {
                // NOTE: allSamples is either allExternals, allInternals, etc...
                hAllSamples[i]->Add((TH1D*)allSamples[j]->At(i));
                if(h_name.Contains("hTotalE"))
                {
                    std::cout << "adding to hAllSamples[" << i << "]: " << ((TH1D*)allSamples[j]->At(i))->GetName() << std::endl;
                }
            }
        }
        else
        {
            // NOTE: allSamples is either allExternals, allInternals, etc...
            hAllSamples[i]->Add((TH1D*)allSamples[j]->At(i));
            if(h_name.Contains("hTotalE"))
            {
                std::cout << "adding to hAllSamples[" << i << "]: " << ((TH1D*)allSamples[j]->At(i))->GetName() << std::endl;
            }
        }

        // 2020-03-31: moved from below block to this location
        tmpHist = (TH1D*)allSamples[j]->At(i);
        if(tmpHist->Integral() > 0)
        {
            //std::cout << ExternalBkgFiles[j] << " " << tmpHist->Integral() << std::endl;
            tmpHist->SetTitle(sampleNames[j]); // TODO: want this above?
            hMinorStacks[i]->Add((TH1D*)allSamples[j]->At(i)->Clone(sampleFiles[j]));
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
                hAllMC[i] = (TH1D*)allSamples[j]->At(i)->Clone("Total MC");
            }
            else
            {
                // NOTE: allSamples is either allExternals, allInternals, etc...
                hAllMC[i]->Add((TH1D*)allSamples[j]->At(i));
            }
        }
        else
        {
            // NOTE: allSamples is either allExternals, allInternals, etc...
            hAllMC[i]->Add((TH1D*)allSamples[j]->At(i));
        }

   
        //(TH1D*)allSamples[j]->At(i)->Write();
        // 2020-03-30: removed
        

    } // generic


    TString events;
    //std::cout << "150Nd: Add()" << std::endl;
    //events.Form("%i",(int)hAllNd150[i]->Integral());
    //hAllNd150[i]->SetTitle("^{150}Nd 2#nu2#beta ("+events+")");
    //hMajorStacks[i]->Add((TH1D*)hAllNd150[i]);
    
    if(hMajorStacks_add == true)
    {
        if(h_name.Contains("hTotalE"))
        {
            std::cout << sampleName << ": Add()" << std::endl;
        }
        events.Form("%i", (int)hAllSamples[i]->Integral());
        TString title_string = sampleHumanName + " (" + events + ")";
        //tmpHist->SetFillColor(dataColors[i]);
        //tmpHist->SetLineColor(dataColors[i]);
        hAllSamples[i]->SetFillColor(majorStackColor);
        hAllSamples[i]->SetLineColor(majorStackColor);
        hAllSamples[i]->SetLineWidth(0);
        hAllSamples[i]->SetTitle(title_string);
        hMajorStacks[i]->Add((TH1D*)hAllSamples[i]);
        
        //std::cout << "Externals: Add()" << std::endl;
        //events.Form("%i",(int)hAllExternals[i]->Integral());
        //hAllExternals[i]->SetTitle("Externals ("+events+")");
        //hMajorStacks[i]->Add((TH1D*)hAllExternals[i]);
    }
    else
    {
        // do nothing
    }
    

}


// TODO re-enable RAWENABLE ?
// need to save into file in different location, or not save at all
void fitHistograms()
{
    // TODO: draw the 150nd 2D histograms for hMCElectronEnergy and
    // hMCTrueElectronEnergy


    std::cout << "fitHistograms()" << std::endl;

    TFile *myFile = TFile::Open("Nd150_2e_P"+Phase+".root", "UPDATE");
    TDirectory* dir = gDirectory; 
    //std::cout << "dir->ls():" << std::endl;
    //dir->ls();
    TObjArray *allDataHistograms = new TObjArray();    
    TObjArray *allDataHistograms_rawdata = new TObjArray();    
    
    //Get the data histograms
    for(int i = 0; i < numHistograms; i++)
    {

        //TString name_append;
        //if(mode_flag == 0)
        //{
        //    //name_append = "";
        //}
        //else if(mode_flag == 1)
        //{
        //    name_append = "_raw";
        //}

        TString histogram_name = histogramNames[i];

        // Note: DataFile = "data_2e"
        // processed
        {
            //std::string name("singleHistos/unscaled/" + histogramNames[i] + DataFile);
            //std::string name("rawdata/" + histogram_name  + "/" + histogramNames[i] + DataFile);
            //std::string name("processeddata/" + histogram_name  + "/" + histogramNames[i] + DataFile + name_append);
            std::string name("processeddata/" + histogram_name  + "/" + histogramNames[i] + DataFile);
            //std::cout << "name=" << name << std::endl;
            if(histogram_name.CompareTo("hTotalE_") == 0)
            {
                std::cout << "Get() : " << name << " from file, Clone() : " << histogramNames[i] + "data" << std::endl;
            }
            allDataHistograms->Add((TH1D*)myFile->Get(name.c_str())->Clone(histogramNames[i] + "data"));
            // 2020-04-02: changed "_data" to "data"
        }

#if RAWENABLE
        // raw
        {
            //std::string name_rawdata("rawdata/" + histogram_name  + "/" + histogramNames[i] + DataFile + name_append);
            std::string name_rawdata("rawdata/" + histogram_name  + "/" + histogramNames[i] + DataFile + "_raw");
            if(histogram_name.CompareTo("hTotalE_") == 0)
            {
                std::cout << "Get() : " << name_rawdata << " from file, Clone() : " << histogramNames[i] + "data_raw" << std::endl;
            }
            allDataHistograms_rawdata->Add((TH1D*)myFile->Get(name_rawdata.c_str())->Clone(histogramNames[i] + "data_raw"));
        }
#endif
    }





#if 0

            std::cout << "DONT DO THIS" << std::endl;
            std::cin.get();

            TString histogram_name = histogramNames[j];
            TString new_histogram_name = histogram_name + sampleFiles[i] + "_fit";
            //TString name = "rawdata/" + histogram_name + "/" + histogram_name + sampleFiles[i];

            if(mode_flag == 0)
            {
                TString name = "processeddata/" + histogram_name + "/" + histogram_name + sampleFiles[i]; // + name_append
                std::cout << "Get() : " << name << " from file, Clone() : " << new_histogram_name << std::endl;
                tmpHist = (TH1D*)myFile->Get(name)->Clone(new_histogram_name);
            }
            else if(mode_flag == 1)
            {
                TString name = "rawdata/" + histogram_name + "/" + histogram_name + sampleFiles[i] + "_raw"; // + name_append;
                std::cout << "Get() : " << name << " from file, Clone() : " << new_histogram_name << std::endl;
                tmpHist = (TH1D*)myFile->Get(name)->Clone(new_histogram_name);
            }
            // NOTE: will have 2 new histograms here
            // 1: scaled, processeddata histogram (save this one to file)
            // 2: scaled, rawdata histogram (do not save this one to file)

            //tmpHist->Sumw2();
            tmpHist->SetFillColor(sampleColors[i]);
            tmpHist->SetLineColor(sampleColors[i]);
            tmpHist->SetLineWidth(0);
            tmpHist->SetTitle(sampleNames[i]);

            tmpHist->Scale(TotalTime / sampleNGenMC[i]);
            if(mode_flag == 0)
            {
                tmpHist->Write();
            }
            tmpHist->Scale(sampleActivity[thePhase * nSamples + i]);

            // NOTE: done AFTER call to Write()
            // so sample is scaled with activity value if used again from allSamples
            allSamples[i]->Add(tmpHist);


#endif















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
    // TH1D objects
    // the named object below is an array of TObjArray pointers
    // the TObjArray is allocated with new in scale
    // there is one for each background type
    // a single TH1D* object is added to the TObjArray
    // so it is an array of length 1
    TObjArray *allExternals[nExternalBkgs];
    TObjArray *allRn220Bkgs[nRn220Bkgs];
    TObjArray *allRn222Bkgs[nRn222Bkgs];
    TObjArray *allInternals[nInternalBkgs];
    TObjArray *allNeighbours[nNeighbours];
    TObjArray *allNd150Samples[nNd150Samples];
    // rawdata versions
    TObjArray *allExternals_rawdata[nExternalBkgs];
    TObjArray *allRn220Bkgs_rawdata[nRn220Bkgs];
    TObjArray *allRn222Bkgs_rawdata[nRn222Bkgs];
    TObjArray *allInternals_rawdata[nInternalBkgs];
    TObjArray *allNeighbours_rawdata[nNeighbours];
    TObjArray *allNd150Samples_rawdata[nNd150Samples];

    // clear this output file
    std::ofstream of_numberofevents("of_numberofevents.txt", std::ofstream::out | std::ofstream::app);

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
    
    #if RAWENABLE
    scale(myFile,
          "../include/activities.txt",
          "internals/",
          nInternalBkgs,
          InternalBkgFiles,
          InternalBkgNames,
          InternalBkgColors,
          AcceptedTime,
          InternalBkgNGenMC_rawdata,
          &InternalBkgActivity_rawdata[0][0],
          &InternalBkgEfficiency_rawdata[0][0],
          allInternals_rawdata, 1);
    #endif
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
    
    #if RAWENABLE
    scale(myFile,
          "../include/activities.txt",
          "neighbours/",
          nNeighbours,
          NeighbourFiles,
          NeighbourNames,
          NeighbourColors,
          AcceptedTime,
          NeighbourNGenMC_rawdata,
          &NeighbourActivity_rawdata[0][0],
          &NeighbourEfficiency_rawdata[0][0],
          allNeighbours_rawdata, 1);
    #endif
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

    #if RAWENABLE
    scale(myFile,
          "../include/activities.txt",
          "nd150/",
          nNd150Samples,
          Nd150Files,
          Nd150Names,
          Nd150Colors,
          AcceptedTime,
          Nd150NGenMC_rawdata,
          &Nd150Activity_rawdata[0][0],
          &Nd150Efficiency_rawdata[0][0],
          allNd150Samples_rawdata, 1);
    #endif
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

    #if RAWENABLE
    scale(myFile,
          "../include/activities.txt",
          "externals/",
          nExternalBkgs,
          ExternalBkgFiles,
          ExternalBkgNames,
          ExternalBkgColors,
          AcceptedTime,
          ExternalBkgNGenMC_rawdata,
          &ExternalBkgActivity_rawdata[0][0],
          &ExternalBkgEfficiency_rawdata[0][0],
          allExternals_rawdata, 1);
    #endif
#endif

#if RADON_ON
    std::cout << "scale: radon" << std::endl;

    #if 0
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

        #if RAWENABLE
        scale(myFile,
              "../include/activities.txt",
              "externals/",
              nRn222Bkgs,
              Rn222BkgFiles,
              Rn222BkgNames,
              Rn222BkgColors,
              AcceptedTime,
              Rn222BkgNGenMC_rawdata,
              &Rn222BkgActivity_rawdata[0][0],
              &Rn222BkgEfficiency_rawdata[0][0],
              allRn222Bkgs_rawdata, 1);
        #endif
    #else
        scale(myFile,
              "../include/activities.txt",
              "externals/",
              nRn222Bkgs,
              Rn222BkgFilesNew,
              //Rn222BkgNamesNew,
              Rn222BkgNames,
              //Rn222BkgColorsNew,
              Rn222BkgColors,
              AcceptedTime,
              Rn222BkgNGenMC,
              &Rn222BkgActivity[0][0],
              &Rn222BkgEfficiency[0][0],
              allRn222Bkgs);

        #if RAWENABLE
        scale(myFile,
              "../include/activities.txt",
              "externals/",
              nRn222Bkgs,
              Rn222BkgFilesNew,
              //Rn222BkgNamesNew,
              Rn222BkgNames,
              //Rn222BkgColorsNew,
              Rn222BkgColors,
              AcceptedTime,
              Rn222BkgNGenMC_rawdata,
              &Rn222BkgActivity_rawdata[0][0],
              &Rn222BkgEfficiency_rawdata[0][0],
              allRn222Bkgs_rawdata, 1);
        #endif
    #endif

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

    #if RAWENABLE
    scale(myFile,
          "../include/activities.txt",
          "externals/",
          nRn220Bkgs,
          Rn220BkgFiles,
          Rn220BkgNames,
          Rn220BkgColors,
          AcceptedTime,
          Rn220BkgNGenMC_rawdata,
          &Rn220BkgActivity_rawdata[0][0],
          &Rn220BkgEfficiency_rawdata[0][0],
          allRn220Bkgs_rawdata, 1);
    #endif
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
    TH1D *hAllExternals[numHistograms];
    TH1D *hAllRadon[numHistograms];
    TH1D *hAllInternals[numHistograms];
    TH1D *hAllNeighbours[numHistograms];
    TH1D *hAllNd150[numHistograms];
    // rawdata versions
    TH1D *hAllExternals_rawdata[numHistograms];
    TH1D *hAllRadon_rawdata[numHistograms];
    TH1D *hAllInternals_rawdata[numHistograms];
    TH1D *hAllNeighbours_rawdata[numHistograms];
    TH1D *hAllNd150_rawdata[numHistograms];

    TH1D *hAllMC[numHistograms];
    TH1D *hAllMC_rawdata[numHistograms];
    // for all MC samples to calculate chi2
    // also used to draw black line around stack of all MC samples

    // stacks for drawing: hMajorStacks: samples are grouped by type,
    // eg: radon, external, internal, neighbour, 150Nd
    // hMinorStacks: samples are not grouped, each individual sample
    // is displayed seperatly
    THStack *hMajorStacks[numHistograms];
    THStack *hMinorStacks[numHistograms];
    THStack *hMajorStacks_rawdata[numHistograms];
    THStack *hMinorStacks_rawdata[numHistograms];

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
    #if RAWENABLE
        hMajorStacks_rawdata[i] = new THStack("majorstack_rawdata" + i_str, histogramNames[i]);
        hMinorStacks_rawdata[i] = new THStack("minorstack_rawdata" + i_str, histogramNames[i]);
    #endif



#if EXTERNALS_ON
        if(i == 0)
        {
            std::cout << "stack: externals" << std::endl;
        }

        stackfunction(i, nExternalBkgs,
                      hMinorStacks, allExternals, // TObjArray *allExternals[nExternalBkgs]
                      hMajorStacks, "External Backgrounds",
                      "Externals", hAllExternals, // TH1D *hAllExt[numHistograms]
                      ExternalBkgNames, ExternalBkgFiles, hAllMC,
                      true, false,
                      ExternalBkgColor,
                      true);

    #if RAWENABLE
        stackfunction(i, nExternalBkgs,
                      hMinorStacks_rawdata, allExternals_rawdata, // TObjArray *allExternals[nExternalBkgs]
                      hMajorStacks_rawdata, "External Backgrounds",
                      "Externals", hAllExternals_rawdata, // TH1D *hAllExt[numHistograms]
                      ExternalBkgNames, ExternalBkgFiles, hAllMC_rawdata,
                      true, false,
                      ExternalBkgColor,
                      true);
    #endif
#endif

#if RADON_ON
        if(i == 0)
        {
            std::cout << "stack: radon" << std::endl;
        }

        // first
        // add to existing
        // add to majors
        #if 0
        stackfunction(i, nRn222Bkgs,
                      hMinorStacks, allRn222Bkgs,
                      hMajorStacks, "Radon Backgrounds",
                      "Radon", hAllRadon,
                      Rn222BkgNames, Rn222BkgFiles, hAllMC,
                      false, false,
                      RadonBkgColor,
                      false);

            #if RAWENABLE
        stackfunction(i, nRn222Bkgs,
                      hMinorStacks_rawdata, allRn222Bkgs_rawdata,
                      hMajorStacks_rawdata, "Radon Backgrounds",
                      "Radon", hAllRadon_rawdata,
                      Rn222BkgNames, Rn222BkgFiles, hAllMC_rawdata,
                      false, false,
                      RadonBkgColor,
                      false);
            #endif
        #else
        if(i == 0)
        {
            std::cout << "stackfunction Rn222" << std::endl;
        }
        stackfunction(i, nRn222Bkgs,
                      hMinorStacks, allRn222Bkgs,
                      hMajorStacks, "Radon Backgrounds",
                      "Radon", hAllRadon,
                      Rn222BkgNames, Rn222BkgFilesNew, hAllMC,
                      #if EXTERNALS_ON
                      false, false,
                      #else
                      true, false,
                      #endif
                      RadonBkgColor,
                      false);

            #if RAWENABLE
        if(i == 0)
        {   
            std::cout << "stackfunction Rn222 (raw)" << std::endl;
        }
        stackfunction(i, nRn222Bkgs,
                      hMinorStacks_rawdata, allRn222Bkgs_rawdata,
                      hMajorStacks_rawdata, "Radon Backgrounds",
                      "Radon", hAllRadon_rawdata,
                      Rn222BkgNames, Rn222BkgFilesNew, hAllMC_rawdata,
                      #if EXTERNALS_ON
                      false, false,
                      #else
                      true, false,
                      #endif
                      RadonBkgColor,
                      false);
            #endif
        #endif


        if(i == 0)
        {
            std::cout << "stackfunction Rn220" << std::endl;
        }
        stackfunction(i, nRn220Bkgs,
                      hMinorStacks, allRn220Bkgs,
                      hMajorStacks, "Radon Backgrounds",
                      "Radon", hAllRadon,
                      Rn220BkgNames, Rn220BkgFiles, hAllMC,
                      false, true,
                      RadonBkgColor,
                      true);

            #if RAWENABLE
        if(i == 0)
        {
            std::cout << "stackfunction Rn220 (raw)" << std::endl;
        }
        stackfunction(i, nRn220Bkgs,
                      hMinorStacks_rawdata, allRn220Bkgs_rawdata,
                      hMajorStacks_rawdata, "Radon Backgrounds",
                      "Radon", hAllRadon_rawdata,
                      Rn220BkgNames, Rn220BkgFiles, hAllMC_rawdata,
                      false, true,
                      RadonBkgColor,
                      true);
            #endif
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

    #if RAWENABLE
        stackfunction(i, nInternalBkgs,
                      hMinorStacks_rawdata, allInternals_rawdata,
                      hMajorStacks_rawdata, "Internal Backgrounds",
                      "Internals", hAllInternals_rawdata,
                      InternalBkgNames, InternalBkgFiles, hAllMC_rawdata,
                      false, false,
                      InternalBkgColor,
                      true);
    #endif
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

    #if RAWENABLE
        stackfunction(i, nNeighbours,
                      hMinorStacks_rawdata, allNeighbours_rawdata,
                      hMajorStacks_rawdata, "Neighbouring Foils",
                      "Neighbours", hAllNeighbours_rawdata,
                      NeighbourNames, NeighbourFiles, hAllMC_rawdata,
                      false, false,
                      NeighbourColor,
                      true);
    #endif
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

    #if RAWENABLE
        stackfunction(i, nNd150Samples,
                      hMinorStacks_rawdata, allNd150Samples_rawdata,
                      hMajorStacks_rawdata, "^{150}Nd 2#nu2#beta",
                      "Nd150", hAllNd150_rawdata,
                      Nd150Names, Nd150Files, hAllMC_rawdata,
                      false, false,
                      Nd150Color,
                      true);
    #endif
#endif

        
    }

    std::cout << "Completed the numHistograms loop" << std::endl;

    // stack functions for 2D histograms
/*
    TH2D *hAllMC2D[numHistograms]; //for all MC samples to calculate chi2
    TH2D *tmp2D;
    for(int i = numHistograms-1; i < numHistograms; i++)
    {
        //stack doesnt support 3D hist
        TString i_str;
        i_str.Form("%i",i);

        #if EXTERNALS_ON
        for( int j = 0; j < nExternalBkgs; j++)
        {
            (TH2D*)allExternals[j]->At(i)->Write();
            if(j == 0)
            {
                // TODO: this one different because of clone option
                hAllMC2D[i] = (TH2D*)allExternals[j]->At(i)->Clone("allMC"+histogramNames[i]);
            }
            else
            {
                hAllMC2D[i]->Add((TH2D*)allExternals[j]->At(i)->Clone());
            }
        }//externals
        #endif
    
        #if RADON_ON
        for(int j = 0; j < nRn222Bkgs; j++)
        {
            (TH2D*)allRn222Bkgs[j]->At(i)->Write();
            hAllMC2D[i]->Add((TH2D*)allRn222Bkgs[j]->At(i)->Clone());
        }//radon222

        for(int j = 0; j < nRn220Bkgs; j++)
        {
            (TH2D*)allRn220Bkgs[j]->At(i)->Write();
            hAllMC2D[i]->Add((TH2D*)allRn220Bkgs[j]->At(i)->Clone());
        }//radon220
        #endif

        #if INTERNALS_ON
        for(int j = 0; j < nInternalBkgs; j++)
        {
            (TH2D*)allInternals[j]->At(i)->Write();
            hAllMC2D[i]->Add((TH2D*)allInternals[j]->At(i)->Clone());
        }//internals
        #endif

        #if NEIGHBOURS_ON
        for(int j = 0; j < nNeighbours; j++)
        {
            (TH2D*)allNeighbours[j]->At(i)->Write();
            hAllMC2D[i]->Add((TH2D*)allNeighbours[j]->At(i)->Clone());
        }//neighbours
        #endif

        #if ND_ON
        for(int j = 0; j < nNd150Samples; j++)
        {
            (TH2D*)allNd150Samples[j]->At(i)->Write();
            hAllMC2D[i]->Add((TH2D*)allNd150Samples[j]->At(i)->Clone());
        }//150Nd
        #endif

        hAllMC2D[i]->Write();

        tmp2D = (TH2D*)allDataHistograms->At(i);
        tmp2D->SetName(histogramNames[i]+DataName);
        tmp2D->Write();
    }
    std::cout << "Finished the cloning loop" << std::endl;
*/



    // loop to write out histograms before running canvas loop
    // also sets hist structure for use in following loop
    TH1D* hist[numHistograms];
    for(int i = 0; i < numHistograms; i++)
    {

        hist[i] = (TH1D*)allDataHistograms->At(i);
        //hist[i]->SetName(histogramNames[i] + DataName);

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
    #if RAWENABLE
    // same for rawdata versions
    TH1D* hist_rawdata[numHistograms];
    for(int i = 0; i < numHistograms; i++)
    {

        hist_rawdata[i] = (TH1D*)allDataHistograms_rawdata->At(i);
        hist_rawdata[i]->SetName(histogramNames[i] + DataName);
        // TODO: bug here: same name

        // inserted from below loop (2)
        TString Ndata_str;
        Ndata_str.Form("%i", (int)hist_rawdata[i]->Integral());
        hist_rawdata[i]->SetTitle("Data (" + Ndata_str + ")");
        hist_rawdata[i]->SetLineWidth(2);
        hist_rawdata[i]->SetMarkerStyle(20);

        // 2020-04-02: removed
        // 2020-04-02: re-enabled, this is required to write the data
        // histogram to file in the root directory
        // 2020-04-02: disabled again, will read data histogram from
        // /singleHistos/unscaled directory in fitting code
        // 2020-04-02: requires fiddling around with GetListOfKeys() and
        // cd() to different TDirectory, can't be bothered with that
        // so will re-enable
        std::cout << "writing histogram with name " << hist_rawdata[i]->GetName() << " at end of code which was removed" << std::endl;
        // TODO: fix bug by not writing this histogram to file
        // TODO: should it go in the rawdata folder?
        // TODO: it already is
        //hist_rawdata[i]->Write();
        // TODO: don't like that other functions write out the equivalent
        // MC histograms to file in a different place from this Write() call
        // expecially since this call to Write() is inside a block of code
        // which creates a TCanvas - move this and SetName() elsewhere

    }
    #endif


/*
    for(std::map<TString, TH1*>::iterator it{histogramPointers.begin();
        it != histogramPointers.end(); ++ it)
    {

        TH1D *histtmp = (TH1D*)it->second;

        std::ofstream of_hTotalE_numberofevents("hTotalE_numberofevents.txt")
        if(TString(histtmp->GetName()).Contains("hTotalE"))
        {
            std::cout << histtmp->GetName() << " number of events " << histtmp->Integral() << std::endl;
            of_hTotalE_numberofevents << histtmp->GetName() << ", " << histtmp->Integral() << std::endl;
        }
    }
*/


    // TODO: draw raw data histograms 


    TCanvas *c_processeddata[numHistograms];
    TLegend *leg_processeddata[numHistograms];

    #if RAWENABLE
    TCanvas *c_rawdata[numHistograms];
    TLegend *leg_rawdata[numHistograms];
    #endif
    
    std::cout << "About to start the canvas loop" << std::endl;
    for(int j = 0; j < numHistograms; j++)
    {

        ////////////////////
        // processed data //
        ////////////////////
        
        {
            TCanvas *&c = c_processeddata[j];
            TLegend *&leg = leg_processeddata[j];

            if(histogramDrawFlag[j] == 0) continue;

            std::cout << "At start of canvas loop" << endl;

            TString i_str;
            TString Nmc_str;
            
            i_str.Form("%i", j);
            c = new TCanvas("c_processeddata" + i_str + "_");
            c->SetFillColor(kWhite);
            // c[i]->SetGrid(1,1);
            //c[i]->SetTicks(1,1);
            //c[i]->Divide(2,2);

            //c[i]->cd(1);
            //c0 = (TPad*)c[i]->GetPad(1);
            // c1_1->SetBottomMargin(0.1);
            //c0->SetGrid(1,1);
            //c0->SetTicks(1,1);
            c->SetGrid(1, 1);
            c->SetTicks(1, 1);

            // moved to above (1)

            Nmc_str.Form("%i", (int)hAllMC[j]->Integral());
            // moved to above (2)
            // draw twice to set max
            if(j == 58)
            {
            //    hist[j]->Draw("colz");
            }
            else
            {
                hist[j]->Draw("PE");
            }

            //std::cout << " data: " << hist[i]->Integral() << std::endl;
            //hSimpleStacks[j]->Draw("hist same");
            if(j == 58)
            {
            }
            else
            {
                //hMinorStacks[j]->Draw("hist");
                hMinorStacks[j]->Draw("hist same");
                //hMajorStacks[j]->Draw("hist same");
            }


            //hAllMC[j]->SetMaximum(1.0e+04);
            //hAllMC[j]->GetYaxis()->SetMaximum(1.0e+04);
            //hAllMC[j]->GetYaxis()->SetRangeUser(0.0, 1.0e+04);

            hAllMC[j]->SetLineWidth(2);
            hAllMC[j]->SetLineColor(kBlack);
            hAllMC[j]->SetFillColor(kBlack);
            hAllMC[j]->SetFillStyle(0);
            //hAllMC[j]->Sumw2();
            if(j == 58)
            {
            }
            else
            {
                hAllMC[j]->Draw("hist same");
            }

            // moved from above
            if(j == 58)
            {
                hist[j]->Draw("colz");
            }
            else
            {
                hist[j]->Draw("PE same");
            }

            double chi2;
            int ndf, igood;
            TString chi2_str, ndf_str;
            double prob = hist[j]->Chi2TestX(hAllMC[j], chi2, ndf, igood, "UW");
            chi2_str.Form("%4.3f", chi2);
            ndf_str.Form("%i", ndf);
            hAllMC[j]->SetTitle("total MC (" + Nmc_str + ")");
            // TODO: moved Ndata_str out to above loop, consider moving out Nmc_str

            //hAllMC[i]->Write();
            // 2020-03-30: removed

            // TODO: forget difference between hAllMC and hMinorStacks
            // why is only 1 written to file? does this write all the
            // histograms to file? cannot be the case, as some
            // histos would be skipped by the continue statement above

            //leg[i] = c0->BuildLegend();
            leg = c->BuildLegend();
            leg->SetName("leg" + i_str + "_");
            // TODO: build legend using parameter_names.list rather than
            // individual MC samples
            leg->SetFillColor(kWhite);
            leg->AddEntry((TObject*)0, "#chi^{2}/ndf = " + chi2_str + "/" + ndf_str, "");
            //c->cd(2);
            TLegend *tmpLeg = (TLegend*)leg->Clone("leg_copy");
            tmpLeg->Draw();
            //c0->cd();
            leg->Delete();
            hist[j]->SetTitle(histogramNames[j]); //switch to titles once you've written them
            /*
            if(j == 58)
            {
            }
            else
            {
                hist[j]->Draw("PEsame");
            }
            */

            //c[i]->cd(3);
            //c2 = (TPad*)c[i]->GetPad(3);
            // c1_1->SetBottomMargin(0.1);
            //c2->SetGrid(1,1);
            //c2->SetTicks(1,1);

            //make ratio plot
            //TH1D *ratio = (TH1D*)hist[i]->Clone("ratio_" + histogramNames[i]);
            //ratio->Sumw2();
            //ratio->Divide((TH1D*)allMC[i]);
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
        }


        //////////////
        // raw data //
        //////////////
        
    #if RAWENABLE
        if(histogramDrawFlag_rawdata == 1)
        {
            TCanvas *&c = c_rawdata[j];
            TLegend *&leg = leg_rawdata[j];

            if(histogramDrawFlag[j] == 0) continue;

            std::cout << "At start of canvas loop" << endl;

            TString i_str;
            TString Nmc_str;
            
            i_str.Form("%i", j);
            c = new TCanvas("c_rawdata" + i_str + "_");
            c->SetFillColor(kWhite);
            // c[i]->SetGrid(1,1);
            //c[i]->SetTicks(1,1);
            //c[i]->Divide(2,2);

            //c[i]->cd(1);
            //c0 = (TPad*)c[i]->GetPad(1);
            // c1_1->SetBottomMargin(0.1);
            //c0->SetGrid(1,1);
            //c0->SetTicks(1,1);
            c->SetGrid(1, 1);
            c->SetTicks(1, 1);

            // moved to above (1)

            Nmc_str.Form("%i", (int)hAllMC_rawdata[j]->Integral());
            // moved to above (2) 
            if(j == 58)
            {
                //hist_rawdata[j]->Draw("colz");
            }
            else
            {
                hist_rawdata[j]->Draw("PE");
            }

            //std::cout << " data: " << hist[i]->Integral() << std::endl;
            //hSimpleStacks[i]->Draw("hist same");
            if(j == 58)
            {
            }
            else
            {
                //hMinorStacks_rawdata[j]->Draw("hist");
                //hMajorStacks_rawdata[j]->Draw("hist same");
                hMinorStacks_rawdata[j]->Draw("hist same");
            }


            //hAllMC_rawdata[j]->SetMaximum(1.0e+04);
            //hAllMC_rawdata[j]->GetYaxis()->SetRangeUser(0.0, 1.0e+04);

            hAllMC_rawdata[j]->SetLineWidth(2);
            hAllMC_rawdata[j]->SetLineColor(kBlack);
            hAllMC_rawdata[j]->SetFillColor(kBlack);
            hAllMC_rawdata[j]->SetFillStyle(0);
            //hAllMC_rawdata[j]->Sumw2();
            if(j == 58)
            {
            }
            else
            {
                hAllMC_rawdata[j]->Draw("hist same");
                //hAllMC_rawdata[j]->Draw("hist");
            }

            // moved to here
            if(j == 58)
            {
                hist_rawdata[j]->Draw("colz");
            }
            else
            {
                hist_rawdata[j]->Draw("PE same");
            }

            double chi2;
            int ndf, igood;
            TString chi2_str, ndf_str;
            double prob = hist_rawdata[j]->Chi2TestX(hAllMC_rawdata[j], chi2, ndf, igood, "UW");
            chi2_str.Form("%4.3f", chi2);
            ndf_str.Form("%i", ndf);
            hAllMC_rawdata[j]->SetTitle("total MC (" + Nmc_str + ")");
            // TODO: moved Ndata_str out to above loop, consider moving out Nmc_str

            //hAllMC[i]->Write();
            // 2020-03-30: removed

            // TODO: forget difference between hAllMC and hMinorStacks
            // why is only 1 written to file? does this write all the
            // histograms to file? cannot be the case, as some
            // histos would be skipped by the continue statement above

            //leg[i] = c0->BuildLegend();
            leg = c->BuildLegend();
            leg->SetName("leg_rawdata" + i_str + "_");
            leg->SetFillColor(kWhite);
            leg->AddEntry((TObject*)0, "#chi^{2}/ndf = " + chi2_str + "/" + ndf_str, "");
            //c[i]->cd(2);
            TLegend *tmpLeg = (TLegend*)leg->Clone("leg_rawdata_copy");
            tmpLeg->Draw();
            //c0->cd();
            leg->Delete();
            hist_rawdata[j]->SetTitle(histogramNames[j]); //switch to titles once you've written them
            /*
            if()
            {
                hist_rawdata[j]->Draw("colz");
            }
            else
            {
                hist_rawdata[j]->Draw("PEsame");
            }
            */

            //make ratio plot
            //TH1D *ratio = (TH1D*)hist[i]->Clone("ratio_" + histogramNames[i]);
            //ratio->Sumw2();
            //ratio->Divide((TH1D*)allMC[i]);
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
        }
    #endif
   

/*
    TCanvas *c[numHistograms];
    TH1D* hist[numHistograms];
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

        hist[i] = (TH1D*)allDataHistograms->At(i);
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
        TH1D *ratio = (TH1D*)hist[i]->Clone("ratio_" + histogramNames[i]);
        ratio->Sumw2();
        ratio->Divide((TH1D*)allMC[i]);
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


void drawPlots(TH1D *data,TH1D *mc, THStack *hs)
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

    TH1D *ratio = (TH1D*)data->Clone("ratio_" + whichHist); 
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
    data->SetTitle("Phase 1OR2"); // TODO ?
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

double getChi2(TH1D *data, TH1D* mc, Int_t &ndof)
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

