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
void makeHistograms(
    TString thePath,
    TString sampleName,
    std::ofstream& ofile_cutcount,
    const uint64_t sample_split_flags = 0);

void loadFiles()
{

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);

    // map of flags, used for selecting events which are split by
    // L0, L>0, IN, OUT
    // 0 = default
    // 1 = select IN events
    // 2 = select OUT events
    // 4 = select L0 IN events
    // 8 = select L0 OUT events
    // 16 = select L>0 (IN+OUT) events
    std::map<std::string, uint64_t> sample_split_flags;
    sample_split_flags["bi214_sfoil"] = (uint64_t)1 | (uint64_t)2;
    sample_split_flags["pb214_sfoil"] = (uint64_t)1 | (uint64_t)2;
    sample_split_flags["bi214_swire"] = (uint64_t)4 | (uint64_t)8 | (uint64_t)16;
    sample_split_flags["pb214_swire"] = (uint64_t)4 | (uint64_t)8 | (uint64_t)16;


    // open and reset file for storing the number of events that pass cuts
    std::ofstream ofile_cutcount("cutcount_splitadditional.txt", std::ofstream::out | std::ofstream::trunc);
    // TODO: need to count the number of events that pass cuts and store
    // into variable ExternalBkgPassCut[thePhase][i]
    // cannot pass that variable here because this would require
    // loadFiles() to be called before fitHistograms()
    // which is slow
    // instead, store results into a file


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
                    makeHistograms("externals/", Rn222BkgFilesInput[i], ofile_cutcount, flag_bit);
                }

                //flags_shift << 1;
                ++ counter;

                //std::cout << std::hex;
                //std::cout << "0xFFFFFFFFFFFFFFFF << counter = " << (0xFFFFFFFFFFFFFFFF << counter) << " flags=" << flags << std::endl;
                //std::cout << std::dec;
            }
        }
        //else
        //{
        //    makeHistograms("externals/", Rn222BkgFilesInput[i], ofile_cutcount);
        //}

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




    /*
    std::cout << std::endl;
    std::cout << ">>>>> Reading in Rn222 backgrounds..." << std::endl;
    for(int i = 0; i < nRn222Bkgs; i++)
    {
        makeHistograms("externals/", Rn222BkgFilesNew[i], ofile_cutcount, 0);
        makeHistograms("externals/", Rn222BkgFilesNew[i], ofile_cutcount, 1);
    }
    std::cout << "... DONE." << std::endl;
    */



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
    const uint64_t sample_split_flags = 0)
{

    // enable / disable additional cuts
    // by default the basic cuts are enabled
    // this includes sector and phase cut

    TString name_sample_split_additional = "";
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
    /*
    if((sampleName.CompareTo("nd150_rot_2b2n_m4") == 0)   ||
       (sampleName.CompareTo("nd150_rot_2n2b_m4") == 0))
    {
        #if TRUTH_ENABLE
        //aFile = TFile::Open(filePath + thePath + sampleName + "/Nd150_2eNg_output_truth.root");
        aFile = TFile::Open(filePath + thePath + sampleName + "/Nd150_2eNg_output_truth_NEW_all.root");
        #else
        aFile = TFile::Open(filePath + thePath + sampleName + "/Nd150_2eNg_output.root");
        #endif
    }
    else
    {
    */
        aFile = TFile::Open(filePath + thePath + sampleName + "/Nd150_2eNg_output.root");
    /*}*/
    TTree *theTree = (TTree*)aFile->Get("Nd150_2eNg/Nd150_2eNg");

    TFile *outputFile = nullptr;
    TTree *outputTree = nullptr;
    TFile *outputFile_small = nullptr;
    TTree *outputTree_small = nullptr;
    //if(mode_flag == 0)
    //{
        /*
        outputFile_small = new TFile(filePath + thePath + sampleName + "/Nd150_2eNg_output_postprocessed_small.root", "recreate");
        outputFile_small->mkdir("Nd150_2eNg");                                                                           
        outputFile_small->cd("Nd150_2eNg");                                                  
        outputTree_small = new TTree("Nd150_2eNg", "Nd150_2eNg");
        */
        outputFile = new TFile(filePath + thePath + sampleName + name_sample_split_additional + "/Nd150_2eNg_output_.root", "recreate");
        outputFile->mkdir("Nd150_2eNg");                                                                           
        outputFile->cd("Nd150_2eNg");                                                  
        outputTree = new TTree("Nd150_2eNg", "Nd150_2eNg");

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
    //}

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
    theTree->SetBranchAddress("nInCluster"                  , nInCluster);
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
   

        ///////////////////////////////////////////////////////////////////////
        // do this for all sample types
        ///////////////////////////////////////////////////////////////////////

        // big tree, all information

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
        outputTree->Branch("nInCluster"                 , nInCluster                    , "nInCluster[nGammaClusters]/I");
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

        // small tree, essential information

        /*
        outputTree_small->Branch("Run"                          , &run                      , "run/I");
        outputTree_small->Branch("eventTime"                    , &eventTime                , "eventTime/D");
        //outputTree_small->Branch("nElectrons"                 , &nElectrons               , "nElectrons/I"); // could remove
        outputTree_small->Branch("radonWeight"                  , &radonWeight              , "radonWeight/D");
        outputTree_small->Branch("bi210Weight"                  , &bi210Weight              , "bi210Weight/D");
        outputTree_small->Branch("electronEnergy"               , electronEnergy            , "electronEnergy[2]/D");
        // truth information, save in all files to maintain homogenaity
        outputTree_small->Branch("trueElectronEnergy"           , trueElectronEnergy        , "trueElectronEnergy[2]/D");
        */

    std::cout << "Processing: " << sampleName << std::endl;

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
        theTree->GetEvent(event_i);
        ++ total_event_count;

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
        //if(mode_flag_2 == 1)
        if(true)
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


        ++ cut_counter[cc]; // cut 0 - input count
        ++ cc;






        // 1: Accepted Run
        // (1) Accepted Run


        // maybe 5.87
        //if(vertexSec[0] < 5.74 || vertexSec[0] > 5.88 || vertexSec[1] < 5.74 || vertexSec[1] > 5.88)
        // The following values taken from Position Paper Section 2,
        // seems to match data from 2e directory
        // No cut in 2eNg, do not know why
        ++ cut_counter[cc]; // cut 2 - vertex sector
        ++ cc;




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

        // Fill histograms with relevant variables which we want to plot pre-cuts

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
        ++ cut_counter[cc]; // cut 5 - track negative curvature
        ++ cc;
        
        // 19: Each L_e >= 30 cm
        // (25) Track lengths > 30 cm
        // 4. - position paper number/index, this index appears multiple times
        // 4. Each of the electron tracks must be at least 30 cm long, [have a hit in the first Gg cell layer,
        // and be associated to unique scintillator blocks.]
        ++ cut_counter[cc]; // cut 6 - track length
        ++ cc;


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
        ++ cut_counter[cc]; // cut 7 (was cut 9) - track hot spot
        ++ cc;
        
        // 10. No unassociated Geiger hits are allowed within 15 cm of the calorimeter hits as this is
        // indicative of a reflected electron, which therefore has not deposited its full energy in the
        // calorimeter.
        // (27) No NA prompt Gg hits within 15 cm of scintillator hits
        // This does not appear to be included in the Table from Summers thesis
        // TODO
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
        ++ cut_counter[cc]; // cut 10 (was cut 11) - Gg hit opposite side
        ++ cc;


        
        
        // 23: Both tracks have hit in Gg L0
        // TODO: copied from above, but this looks like something which should
        // go here - might need to check against layer 1 and 11 ?
        //if((electronFirstGG[0] > 1) || (electronFirstGG[1] > 1)) cut = true; TODO: changed from this?
        // (not 4, 18) - position paper number/index, this index appears multiple times
        // (30) Both tracks have hit in Gg L0
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



        // 24: P_{int} >= 0.01 and P_{ext} <= 0.01
        // (31) P_int > 0.01 and P_ext < 0.01
        // 12. The internal probability between the two electrons must be at least 1% and the external
        // probability must be no more than 1%
        ++ cut_counter[cc]; // cut 12 (was cut 13) - prob
        ++ cc;
        
        // 25: Delta R <= 4cm and Delta Z <= 8 cm
        // (32) Delta R <= 4 cm and Delta Z <= 8 cm
        // 2. The separation between each individually reconstructed vertex must be ≤ 4 cm radially
        // and ≤ 8 cm longitudinally to ensure that the electrons originate from a common vertex.
        ++ cut_counter[cc]; // cut 13 (was cut 14) - delta
        ++ cc;
        
        // save data to file if electronEnergy >= 0.2 MeV
        // and phase is P1 or P2
        outputTree->Fill();
        //outputTree_small->Fill();


        ++ cut_counter[cc]; // cut 14 (was cut 15) - energy
        ++ cc;

        ++ cut_counter[cc]; // cut 1 - phase
        ++ cc;

        ++ event_pass_count;
        //std::cout << "Total number of cuts recorded in cut_counter: " << cc << std::endl;

    } //~event_i
     

    ofile_cutcount  << sampleName << "    " << "" << total_event_count << "    " << "" << event_pass_count << std::endl;
    std::cout       << sampleName << "    " << "total_event_count=" << total_event_count << "    " << "event_pass_count=" << event_pass_count << std::endl;
    //std::cout << "electron_cut_A=" << electron_cut_A << std::endl;
    std::cout << std::endl;

    //if((sampleName.CompareTo("nd150_rot_2b2n_m4") == 0)   ||
    //   (sampleName.CompareTo("nd150_rot_2n2b_m4") == 0))
    //{
    outputTree->Write();
    outputFile->Close();

    //outputTree_small->Write();
    //outputFile_small->Close();
    //}

    aFile->Close();
    
    std::string cut_description[nCuts];
    int cc = 0;
    cut_description[cc] = "Input count";
    ++ cc;
    cut_description[cc] = "IN/OUT/Layer Sample Split";
    ++ cc;
    cut_description[cc] = "none";
    ++ cc;
    cut_description[cc] = "none";
    ++ cc;
    cut_description[cc] = "none";
    ++ cc;
    cut_description[cc] = "none";
    ++ cc;
    cut_description[cc] = "none";
    ++ cc;
    cut_description[cc] = "none";
    ++ cc;
    cut_description[cc] = "none";
    ++ cc;
    cut_description[cc] = "none";
    ++ cc;
    cut_description[cc] = "none";
    ++ cc;
    cut_description[cc] = "none";
    ++ cc;
    cut_description[cc] = "none";
    ++ cc;
    cut_description[cc] = "none";
    ++ cc;
    cut_description[cc] = "none";
    ++ cc;
    cut_description[cc] = "none";
    ++ cc;
    cut_description[cc] = "none";

    //std::cout << "Here are the cut counts" << std::endl;
    std::string of_cutcount_filename;
    of_cutcount_filename = "of_cutcount_processeddata_splitadditional_" + std::string(sampleName) + ".txt";
    std::ofstream of_cutcount(of_cutcount_filename, std::ios::out);
    of_cutcount << sampleName << std::endl;
    for(int i = 0; i < nCuts; i++)
    {
        //std::cout << "Passed cut " << i << ": " << cut_counter[i] << " events. Description=\"" << cut_description[i] << "\"" << std::endl;
        of_cutcount << "Passed cut " << i << ": " << cut_counter[i] << " events. Description=\"" << cut_description[i] << "\"" << std::endl; 
    }
    of_cutcount << std::endl;
    
}


