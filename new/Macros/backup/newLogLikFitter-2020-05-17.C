#include "TH1F.h"
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
#include <TMinuit.h>


// 
#include <time.h>


#include "newLogLikFitter.h"
#include "newLogLikFitter_read.h"


// TODO:
//
// - call loglikelihood function to calculate chisquare value for different
// values of the 100Mo bb background, to produce chisquare curve plot
// to check penalty terms are working correctly
//
// - do this again with 100Mo as a hard parameter and complete minimization
// each step
//
// - sector change cut or remove a significant background and see if fit
// values come out the same
//
// - remove backgrounds 1 by 1 to see difference in chisquare value
//
// - adjust the cut on sector to see difference in results
//
// - add in spectral shape adjustment for gA measurement
//
// - TODO: minuit conversion problem
//  > not scaling histogram MC samples by activity leads to very large
//    range of parameter values minuit has to work with (large spread
//    in orders of magnitude)
//  > may work better if samples were initially scaled by the activities
//    such that minuit parameters are all of order 1
//  > it is BETTER to do this in this code rather than the previous
//    preprocessing code, because this allows dynamic changing of these
//    activity inputs without resulting in some dodgy/arbitary scaling
//    factor which is present from the previous step
//
// - TODO: if parameter is DISABLED it should NOT appear in the fixed_params,
//   free_params arrays
//
// - TODO: are combined parameters (bi/pb mylar etc) being drawn/fit correctly?
//
// - TODO: job for tomorrow (2020-04-22)
//   create a header file with the definitions of background files, names and
//   colors, and include this as a new header to both the fitting code and
//   pre-processing code
//   
//   TODO: there is another bug, I disabled 210Bi SFoil and 228Ac PMT
//   and yet these still appear in the fit, however the blank col/row
//   still appears in the correlation matrix
//   is this still a bug?
//   TODO: now 234Pam is showing up as a zero events MC
//   fixed: also has zero events in it
//
//   make the order the same as summers correlation matrix
//
//   TODO:
//      > names (for covariance matrix, and legend on Total Energy plot)
//      > set default name using "_0_0_" parameter naming convention
//      > or if this cannot be done (due to disabled parameters) then choose
//      > an alternative logical naming convention
//      > eg; "0" (just the parameter number)
//      > this is the "simple form"
//      > THEN: change parameter name to "nd150_rot_2b2n" etc
//      > since these are available and constructed as unique parameter names
//      > shortly after the parameter number is read
//      > these names are slightly more non-trivial to build, as they are done
//      > inside a for loop, so only change these once they are built
//      > this is the "complex form"
//      > FINALLY: read "NAME" specifiers from file and set the name
//      > only IF the name has not been set before
//      > this will require storing a flag, or storing the origional name(s)
//      > (possibly both the simple and complex forms)
//      > TODO
//   TODO: TLegend names for TotalE
//
//
//  TODO: group all Radon backgrounds together because these might be
//        reducible with improved Radon filtering, this includes the
//        Mylar backgrounds


//
// Dear Future Programmer,
//
// You are probably wondering why is does this code resemble a contraption
// of random arrays, maps and indices? The answer is that the previous PhD
// student who worked on this analysis cera 2016 at some point towards 
// the end of their PhD copied their code to their own personal laptop,
// before returning this laptop to the department at the end of their PhD.
// They claim that the laptop, or at least the data on it was destroyed, and
// that they made, however subsequently lost, another backup which was on
// an external hard drive.
//
// Therefore I started with a half finished code, which was riddled with
// obscure bugs, contained no error-sanity checking, and was completely void
// of the use of functions, as if the student did not posess the technical
// ability to write a function for themselves.
//
// I have cleaned up the code substantially, mostly though the use of
// functions, which allowed me to track down many obscure problems. However,
// the insane system of using fixed size arrays to hold all relevant data
// still exists, because I did not change it early on during the development,
// as I saw this as a lower priority task. Although I hate it because it is
// clearly an inflexible way to program such a thing, in some ways I like the
// simplicity of it as a solution, although it is extremely confusing
// syntactically to see things like arrays of strings and arrays of vectors.
//
// I improved things substantially by including new std::map objects, since
// these arrays were essentially used like a map, to map integer parameter
// index objects to other objects, such as strings, or vectors of strings.
//
// However I never fundamentally managed to remove all the fixed length arrays.
// Part of the reason why this is difficult to do is that there is another file
// called fit_2e.C, as well as many other similar files, which are written
// in a similar way to how this code was written. (The same student wrote them
// so this is expected.) Due to this there are some dependencies between how
// the different files are expected to work, and so in order not to risk
// breaking things in a way that I could not spot or quickly repair, I left
// this "feature" in the code...
//
// I would recommend leaving it alone. I have many times considered creating
// a set of classes to manage parameters, which would make a lot more sense.
//
// Additionally, I was the one who introduced the new "parameter_names.lst"
// file format. Unfortunatly, this file format, while a lot more feature rich
// and easy to use than the previous one I inhereted, has not been made to work
// with the pre-processing fit_2e.C code.
//
// I am writing this note because I am at a point in the development where
// making everything work "properly" would require me to rip out the entire
// foundations of the code, and re-write most of the lines from scratch, which
// I have not the time to do.
//
// Sorry.
//
// Ed Bird
// 2020-04-22
//


//-----------------------------------------
//    Functions
//----------------------------------------
void loadFiles();

void book1DHistograms(Int_t channel_counter, TString theChannel,TString thePhase, TString theHistogram);
void book2DHistograms(Int_t channel_counter, TString theChannel,TString thePhase, TString theHistogram);
void logLikelihood(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */); 
Double_t getNumberMC1D(Int_t channel, Int_t binx, Double_t *p);
Double_t getNumberMC2D(Int_t channel, Int_t binx, Int_t biny, Double_t *p);

//void fitBackgrounds(double *AdjustActs, double *AdjustActs_Err, double*& CovMatrix, int& number_free_params, Int_t thePhase);
TMinuit * fitBackgrounds(double *AdjustActs, double *AdjustActs_Err, double*& CovMatrix, int& number_free_params, Int_t thePhase);














void loadFiles()
{

    // load data
    std::cout << "attempting to load spectral data from file" << std::endl;
    //TTree *data_G0;
    //TTree *data_G2;
    //Long64_t count_G0;
    //Long64_t count_G2;
    //LD_REWEIGHT_DATA(data_G0, data_G2, count_G0, count_G2);

    // global
    //h_nEqNull = nullptr;
    //h_nEqTwo = nullptr;
    //h_nEqNull = new TH2D("h_nEqNull", "nEqNull", dimension_xy, 0.0, bb_Q, dimension_xy, 0.0, bb_Q);
    //h_nEqTwo = new TH2D("h_nEqTwo", "nEqTwo", dimension_xy, 0.0, bb_Q, dimension_xy, 0.0, bb_Q);
    //h_nEqNull->SetStats(0);
    //h_nEqTwo->SetStats(0);
    bb_Q = 3.368;
    //std::size_t count_G0 = LD_REWEIGHT_DATA_2(h_nEqNull, "../../data/150Nd-data/dG150Nd/G0/dG0.dat", "h_nEqNull", "nEqNull", 0.0, bb_Q);
    //std::size_t count_G2 = LD_REWEIGHT_DATA_2(h_nEqTwo,  "../../data/150Nd-data/dG150Nd/G2/dG2.dat", "h_nEqNull", "nEqTwo",  0.0, bb_Q);
    // ramdisk, should be faster?
    //std::size_t count_G0 = LD_REWEIGHT_DATA_2(h_nEqNull, "/mnt/ramdisknd150/data/150Nd-data/dG150Nd/G0/dG0.dat", "h_nEqNull", "nEqNull", 0.0, bb_Q);
    //std::size_t count_G2 = LD_REWEIGHT_DATA_2(h_nEqTwo,  "/mnt/ramdisknd150/data/150Nd-data/dG150Nd/G2/dG2.dat", "h_nEqNull", "nEqTwo",  0.0, bb_Q);
    std::size_t count_G0 = LD_REWEIGHT_DATA_2(h_nEqNull, "/home/ecb/150Nd/150Nd-data/dG150Nd/G0/dG0.dat", "h_nEqNull", "nEqNull", 0.0, bb_Q);
    std::size_t count_G2 = LD_REWEIGHT_DATA_2(h_nEqTwo,  "/home/ecb/150Nd/150Nd-data/dG150Nd/G2/dG2.dat", "h_nEqTwo", "nEqTwo",  0.0, bb_Q);

    // phase space integrals
    const Double_t G0_ps_integral_MeV = 0.420438E-45;
    const Double_t G0_ps_integral_yrinv = 0.201577E-16;
    const Double_t G2_ps_integral_MeV = 0.744684E-45;
    const Double_t G2_ps_integral_yrinv = 0.357034E-16;

    // TODO: what is the value for the basline for Nd150? is this for 100Mo
    const Double_t xi_31_baseline{0.368}; // TODO: this is WRONG change it
    // TODO: also change in parameter list file
    //Double_t xi_31_init = 0.8;
    Double_t xi_31_init = 0.368; // change to baseline value for testing purposes

    ///*const Double_t*/ bb_Q = 3.368;
    double count = 0;
    if(count_G0 == count_G2)
    {
        count = count_G0;
    }
    else
    {
        std::cout << "error: count_G0=" << count_G0 << ", count_G2=" << count_G2 << std::endl;
    }
    /*
    //const double A = 0.5;
    //const double B = 0.5;
    //double C = -(double)count;
    //double dimension_xy = (-B + std::sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
    double dimension_xy = (int)(-0.5 + std::sqrt(0.25 - 2.0 * (-(double)count)));
    if(dimension_xy != 3371)
    {
        std::cout << "error! dimension_xy != 3371, dimension_xy=" << dimension_xy << std::endl;
    }
    //const double dimension_xy{3371};

    std::cout << "dimension_xy=" << dimension_xy << std::endl;
    */

    //std::vector<std::vector<double>> data_nEqNull;
    //std::vector<std::vector<double>> data_nEqTwo;
    //read_data(data_G0, data_nEqNull);
    //read_data(data_G2, data_nEqTwo);
    //std::cout << "vector format constructed" << std::endl;
    
    // convert data to histogram format
    /*TH2D **///h_nEqNull = nullptr;
    /*TH2D **///h_nEqTwo = nullptr;
    //convert_data_to_histogram_format(data_nEqNull, data_nEqTwo, h_nEqNull, h_nEqTwo, dimension_xy, bb_Q);
    psiN0 = G0_ps_integral_MeV;
    psiN2 = G2_ps_integral_MeV; // TODO: check this is the correct option
    std::cout << "histogram format constructed" << std::endl;


    // debug
    //TCanvas *c_nEqNull;
    //c_nEqNull = new TCanvas("c_nEqNull", "c_nEqNull"); //, 4000, 3000);
    //h_nEqNull->Draw("colz");
    //c_nEqNull->SaveAs("c_nEqNull.png");
    //c_nEqNull->SaveAs("c_nEqNull.pdf");
    //c_nEqNull->SaveAs("c_nEqNull.C");
    //delete c_nEqNull;

    //TCanvas *c_nEqTwo;
    //c_nEqTwo = new TCanvas("c_nEqTwo", "c_nEqTwo"); //, 4000, 3000);
    //h_nEqTwo->Draw("colz");
    //c_nEqTwo->SaveAs("c_nEqTwo.png");
    //c_nEqTwo->SaveAs("c_nEqTwo.pdf");
    //c_nEqTwo->SaveAs("c_nEqTwo.C");
    //delete c_nEqTwo;
    // debug

    //std::cin.get();


    std::cout << std::fixed << std::endl;

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);

    allDataSamples1D = new TObjArray();
    allDataSamples2D = new TObjArray();



    // read parameter_name.lst file
    //for(int i = 0; i < numberParams; i++)
    //{
        // paramNameMap[i] = std::vector<TString>();
        //paramActMap[i] = std::vector<double>();
        //paramActErrMap[i] = std::vector<double>();

        /*paramNameMap[i].clear();*/
        //paramActMap[i].clear();
        //paramActErrMap[i].clear();
        //paramActP2Map[i].clear();
        //paramActErrP2Map[i].clear();
        // old names above, new below
        /*paramInitValueP1Map[i].clear();
        paramInitErrorP1Map[i].clear();
        paramInitValueP2Map[i].clear();
        paramInitErrorP2Map[i].clear();
        paramConstraintValueP1Map[i].clear();
        paramConstraintErrorP1Map[i].clear();
        paramConstraintValueP2Map[i].clear();
        paramConstraintErrorP2Map[i].clear();*/

    MCNameToParamNameMap.clear();
    MCNameToParamNumberMap.clear();

    paramNameToNumberMap.clear();

    paramNumberToMinuitParamNumberMap.clear();
    minuitParamNumberToParamNumberMap.clear();
    minuitParamNumberCounter = 0;

    numberEnabledParams = 0;

    fixed_params.clear();
    free_params.clear();
    //fixed_param_names.clear();
    //free_param_names.clear();
    //index_free_params.clear();

    //paramNameToHumanReadableParamNameMap.clear();
    MCSampleNameToHumanReadableMCSampleNameMap.clear();



    enabled_params.clear();
    disabled_params.clear();
    //enabled_param_names.clear();
    //disabled_param_names.clear();
    //index_enabled_params.clear();

        // TODO: this is done by calling "new"
        /*
        for(int i = 0; i < number1DHists; ++ i)
        {
            allMCSamples1D[i]->Clear();
        }
        for(int i = 0; i < number2DHists; ++ i)
        {
            allMCSamples2D[i]->Clear();
        }
        allDataSamples1D->Clear();
        allDataSamples2D->Clear();
        */
    //}

    

    // read parameter list file
    read_parameter_list_file();





    // First, create a root file to hold all of the histograms
    TFile *myFile = TFile::Open("Nd150_loglikResults.root", "RECREATE");
    //myFile->Close();
    // 2020-04-02: removed call to close, no idea why this was here
    // or if it is supposed to be here for some reason?

    // Read in all the histograms we want
    //book2DHistograms(0,"1e1g_","P1","hEeVEg_");
    //book2DHistograms(0,"1e1g_","P2","hEeVEg_");
    //book1DHistograms(0,"1e1g_","P1","hClusterEnergy_");
    //book1DHistograms(1,"1e1g_","P2","hClusterEnergy_");
    
    //book1DHistograms(0,"1e_","P1","hEe_");
    //book1DHistograms(0,"1e_","P2","hEe_");
    //book1DHistograms(2,"1e1g_","P1","hTotalE_");
    //book1DHistograms(3,"1e1g_","P2","hTotalE_");
    //book1DHistograms(4,"1e2g_","P1","hClusterEnergyMax_");
    //book1DHistograms(5,"1e2g_","P2","hClusterEnergyMax_");
    //book1DHistograms(2,"1e1gExt_","P1","hClusterEnergy_");
    //book1DHistograms(3,"1e1gExt_","P2","hClusterEnergy_");
    //book1DHistograms(6,"2e_","P1","hTotalE_");
    //book1DHistograms(7,"2e_","P2","hTotalE_");
    //book1DHistograms(4,"OCE_","P1","hTotalE_");
    //book1DHistograms(5,"OCE_","P2","hTotalE_");
    //book1DHistograms(0, "2e_", "P1", "hTotalE_");
    //book1DHistograms(0, "2e_", "P2", "hTotalE_");
    book1DHistograms(0, "2e_", "P" + Phase, "hTotalE_");
  
    // Array to hold activity adjustment parameters
    //Double_t AdjustActs[numberParams];
    //Double_t AdjustActs_Err[numberParams];
    //for(Int_t ix{0}; ix < numberParams; ++ ix)
    Double_t AdjustActs[numberEnabledParams];
    Double_t AdjustActs_Err[numberEnabledParams];
    for(Int_t ix{0}; ix < numberEnabledParams; ++ ix)
    {
        AdjustActs[ix] = 1.0;
        AdjustActs_Err[ix] = 1.0;
    }
    // TODO: fix this
    AdjustActs[1] = xi_31_init;

    // TODO: MARKER1

    /*
    Int_t number_free_params = minuit->GetNumFreePars();
    Double_t *CovMatrix = new Double_t[number_free_params * number_free_params];
    for(int ix{0}; ix < number_free_params * number_free_params; ++ ix)
    {
        CovMatrix[ix] = 0.;
    }
    */
    
    int number_free_params = -1;
    double *CovMatrix = nullptr;

    //fitBackgrounds(AdjustActs, AdjustActs_Err, CovMatrix, number_free_params, thePhase);
    // fit of backgrounds disabled.
    // need to fit Nd150 parameter only
    // disabled in fitBackgrounds function
    
    TMinuit *minuit = fitBackgrounds(AdjustActs, AdjustActs_Err, CovMatrix, number_free_params, thePhase);


    #if 0
    ///////////////////////////////////////////////////////////////////////////
    // testing
    
    // run chisquare tests

    std::cout << "running chi-square tests (gA): " << "variable: g_A parameter (1)" << std::endl;

    int n_tests = 10;
    // 100 Mo
    int axial_vector_parameter_0_index = paramNumberToMinuitParamNumberMap.at(1);
    std::cout << "the internal index for parameter 1 is " << axial_vector_parameter_0_index << std::endl;
    // These are in units of minuit internal parameter units
    // To convert to external parameter units, multiply by the value of the
    // external input parameter initial activity
    // Caution: For cases where the fitted parameter minimum is not at 1.0
    // the errors must be treated as upper and lower bound separatly by adding
    // them (internal param & error) to the central value fit parameter
    // external_param_error_lowerbound = (internal_param_CV - internal_param_error) * external_param_init_value
    // similar for upperbound, then subtract and / 2.0
    double test_central_value = AdjustActs[axial_vector_parameter_0_index];
    double test_range = 0.0; //10.0 * AdjustActs_Err[axial_vector_parameter_0_index];
    // this range should hit delta sigma = 1.0 at 66 % of the width, but it
    // doesn't.
    double test_start = test_central_value - 0.5 * test_range;
    double test_end =   test_central_value + 0.5 * test_range;
    double *test_values = new double[n_tests];
    double test_step = test_range / (double)n_tests;
    std::cout << "test_central_value=" << test_central_value << "\n"
              << "test_range=" << test_range << "\n"
              << "test_start=" << test_start << "\n"
              << "test_end=" << test_end
              << std::endl;
    int n_params = minuit->GetNumPars();
    double *params = new double[n_params];
    double *param_errs = new double[n_params];
    for(int jx = 0; jx < n_params; ++ jx)
    {
        minuit->GetParameter(jx, params[jx], param_errs[jx]);
    }
    std::ofstream ofstream_testvalue("testvalue_gA.txt");
    for(int ix = 0; ix < n_tests; ++ ix)
    {
        test_values[ix] = test_start + test_step * ix;

        // get chisquare value for test
        double fval = 0.;

        // set parameter for 100Mo
        double test_value = test_values[ix];
        params[axial_vector_parameter_0_index] = test_value;
        
        std::cout << "test: ix=" << ix << ", " << "test_value=" << test_value << ", " << "fval=" << fval << std::endl;

        // TODO: reenable
        logLikelihood(n_params, nullptr, fval, params, 0);

        ofstream_testvalue << "value," << test_value << ",chisquare," << fval << std::endl;

        //void logLikelihood(Int_t & /*nPar*/, Double_t* /*grad*/, Double_t &fval, Double_t *p, Int_t /*iflag */)

    }
    ofstream_testvalue.close();
    delete [] test_values;
    //delete [] params;
    //delete [] param_errs;
    #endif


    ///////////////////////////////////////////////////////////////////////////

    std::cout << "The following adjustments (in minuit parameter units) should be made:" << std::endl;
    std::cout << "Note that gA (1) is a special parameter" << std::endl;
    //for(int i = 0; i < numberParams; i++)
    for(int i = 0; i < numberEnabledParams; i++)
    {
        std::cout << i << " :\t" << AdjustActs[i] << " +- " << AdjustActs_Err[i] << std::endl;  
        std::ofstream myFileFitResults("fit_results.txt", std::ios::out | std::ios::app);
        myFileFitResults << AdjustActs[i] << " +- " << AdjustActs_Err[i] << std::endl;

        myFileFitResults.close();
    }

    std::cout << "The following adjustments (in units of Bq) should be made:" << std::endl;
    //for(int i = 0; i < numberParams; i++)
    //std::cout << "numberEnabledParams=" << numberEnabledParams << std::endl;
    for(int i = 0; i < numberEnabledParams; i++)
    {
        //if(i == 1)
        //{
        //    std::cout << "skipping gA parameter TODO FIX" << std::endl;
        //}

        int j = minuitParamNumberToParamNumberMap.at(i);
        //std::cout << "i=" << i << " j=" << j << std::endl;
        Double_t param_init_value = 0.;
        Double_t param_init_error = 0.; 
        if(thePhase == 0)
        {
            param_init_value = paramInitValueP1Map[j];
            param_init_error = paramInitErrorP1Map[j];
        }
        else if(thePhase == 1)
        {
            param_init_value = paramInitValueP2Map[j];
            param_init_error = paramInitErrorP2Map[j];
        }
        else
        {
            std::cout << "ERROR: Invalid value for thePhase: thePhase=" << thePhase << std::endl;
        }
        //std::cout << "value=" << param_init_value << " err=" << param_init_error << " AdjustActs[i]=" << AdjustActs[i] << std::endl;
        if(i != 1)
        {
            std::cout << i << " :\t" << AdjustActs[i] * param_init_value
                           << " +- " << AdjustActs_Err[i] * param_init_value;
            // TODO: put the mutiplication by xi_31_init INSIDE the reweight/fit functions,
            // to restore uniformity in minuit parameters
        }
        else
        {
            std::cout << i << " :\t" << AdjustActs[i]
                           << " +- " << AdjustActs_Err[i];
        }
        Double_t change = 0.0;
        if(i != 1)
        {
            change = 100.0 * (AdjustActs[i] - 1.0);
        }
        else
        {
            change = 100.0 * (AdjustActs[i] - xi_31_init);
        }
        if(change >= 0)
        {
            std::cout << " -> +";
        }
        else
        {
            std::cout << " -> -";
        }
        std::cout << std::abs(change) << " %" << std::endl;  
    }

    {
        std::string summary_fname("fitsummary.txt");
        std::cout << "writing to " << summary_fname << std::endl;
        std::ofstream summary_ofstream(summary_fname.c_str(), std::ios::out);
        
        //for(int i = 0; i < numberParams; ++ i)
        for(int i = 0; i < numberEnabledParams; ++ i)
        {
            int j = minuitParamNumberToParamNumberMap.at(i);

            if(thePhase == 0)
            {
                summary_ofstream << i << ","
                                 << j << ","
                                 << paramInitValueP1Map[j] << ","
                                 << paramInitErrorP1Map[j] << ","
                                 << AdjustActs[i] << ","
                                 << AdjustActs_Err[i] << ","
                                 << (AdjustActs[i] - paramInitValueP1Map[j]) / paramInitValueP1Map[j] * 100.0 << " %" << std::endl;
            }
            else if(thePhase == 1)
            {
                summary_ofstream << i << ","
                                 << j << ","
                                 << paramInitValueP2Map[j] << ","
                                 << paramInitErrorP2Map[j] << ","
                                 << AdjustActs[i] << ","
                                 << AdjustActs_Err[i] << ","
                                 << (AdjustActs[i] - paramInitValueP2Map[j]) / paramInitValueP2Map[j] * 100.0 << " %" << std::endl;
            }
            else
            {
                std::cout << "ERROR: Invalid value for thePhase: thePhase=" << thePhase << std::endl;
            }
        }
    }

    
    {
        std::string chisquarelog_fname("chisquarelog.txt");
        std::cout << "writing to " << chisquarelog_fname << std::endl;
        std::ofstream chisquarelog_ofstream(chisquarelog_fname.c_str(), std::ios::out | std::ios::app);
        chisquarelog_ofstream << std::endl;
        time_t rawtime;
        time(&rawtime);
        struct tm *info;
        info = localtime(&rawtime);
        char timestringbuf[100];
        size_t ncount = strftime(timestringbuf, sizeof(timestringbuf), "%F %T", info);
        std::string datetimestamp(timestringbuf);
        chisquarelog_ofstream << datetimestamp << std::endl;
        int mo100bbnumber = paramNumberToMinuitParamNumberMap.at(10);
        //std::cout << "CHECK CHECK CHECK: mo100bbnumber=" << mo100bbnumber << std::endl;
        chisquarelog_ofstream << paramNameMap[mo100bbnumber] << ", " << AdjustActs[mo100bbnumber] << ", " << global_chisquare << ", (should be values for 100Mo bb)" << std::endl;
    }

  


    ///////////////////////////////////////////////////////////////////////////
    // draw result
    ///////////////////////////////////////////////////////////////////////////


    THStack *stacks1D[number1DHists];
    TCanvas *c;
    TH1F *hAllMC1D[number1DHists];
    TH1F *data1D[number1DHists];

    // TODO: copy over legend code from fit_2e.C

    std::cout << "debug: number of data samples: " << allDataSamples1D->GetEntries() << std::endl;
    std::cout << "debug: number of MC samples: " << allMCSamples1D[0]->GetEntries() << std::endl;


    // each channel 1D hists
    // this is for(i = 0; i < 1; ++ i)
    for(int i = 0; i < allDataSamples1D->GetEntries(); i++)
    {
        if(i != 0)
        {
            std::cout << "debug: ERROR: i = " << i << " should ONLY be 0" << std::endl;
        }

        data1D[i] = (TH1F*)allDataSamples1D->At(i)->Clone();

        TString i_str;
        i_str.Form("%i", i);
    
        c = new TCanvas("c" + i_str);
        c->SetFillColor(kWhite);

        
        stacks1D[i] = new THStack("stacks1D" + i_str, i_str);


        TH1F *tmpHist_draw1D;

        // j list MC samples for this channel i
        //std::cout << "debug: number of MC samples (i=" << i << "): " << allMCSamples1D[i]->GetEntries() << std::endl;

        /*
        std::cout << "list of all objects in allMCSamples1D[i=" << i << "]" << std::endl;
        for(int j = 0; j < allMCSamples1D[i]->GetEntries(); j++)
        {
            //std::cout << (((TH1F*)allMCSamples1D[i])->At(j))->GetName() << std::endl;
            tmpHist_draw1D = (TH1F*)allMCSamples1D[i]->At(j)->Clone();
            TString tmpName = tmpHist_draw1D->GetName();
            std::cout << tmpName << std::endl;
        }
        */

        // what does this code actually do?
        // loops through all MC samples using index j
        // gets name of MC sample into tmpName
        // eg: WHAT
        // loops over 0..numberParams-1 with k
        // loops over 0..paramNameMap[k].size()-1 with n
        // checks if tmpName contains paramNameMap[k].at(n)
        // what is paramNameMap[k].at(n)?
        // it is: std::string paramNameMap[numberParams]
        // this code appears to be broken - was not updated from a previous
        // version perhaps?
        // allMCSamples1D[0] contains objects such as: "zr96_rot_k40_2e_P2"
        
        

        /*
        std::cout << "contents of map: MCNameToParamNumberMap" << std::endl;
        for(auto it = MCNameToParamNumberMap.cbegin(); it != MCNameToParamNumberMap.cend(); ++ it)
        {
            std::cout << it->first << " -> " << it->second << std::endl;
        }
        */

        // MARKER3
        for(int j = 0; j < allMCSamples1D[i]->GetEntries(); j++)
        {

            //std::cout << "j=" << j << std::endl;

            TString j_str;
            j_str.Form("%i", j);

            tmpHist_draw1D = (TH1F*)allMCSamples1D[i]->At(j)->Clone();
            TString tmpName = tmpHist_draw1D->GetName();

            //std::cout << "(1) tmpName=" << tmpName << std::endl;

            //std::cout << "looking for " << tmpName << std::endl;
            int which_param = -1;
            bool found_param = false;

            // search through parameters to find right one
            // this code is broken: it is SUPPOSED to select k such that
            // which_param = k
            // using tmpName.Contains(paramNameMap[k].at(i))
            // which_param corresponds to an MC sample (it is an index, or
            // in other words it is paramNumber)
            // paramNameMap was expected to be a list of parameter names/MC
            // samples corresponding to parameter number/index k
            // in other words, I have a data structure called
            // MCNameToParamNumberMap which maps a parameterMCSample to the
            // parameter index. let's try it and see if it works
            // what is it trying to do? use the names of the histograms to 
            // map what to what? name of histogram to parameter number?
            // so that parameter number can be used to get fitted activity
            // for scaling ? - yes this is what it is trying to do, but
            // includes an error check for if the parameter is not found
            // the histogram names are formatted like:
            // hTotalE_bi214_mylar_fit
            // histogram_name + "_" + mc_sample_name + "_fit"
            
            //std::cout << "NEW CODE" << std::endl;
            //try
            //{
            // TODO: remove TString

            // used later
            double activity_scale_branching_ratio = 1.0;

            {
                std::string tmp_hist_name(tmpName);
                auto i_start = tmp_hist_name.find('_') + 1;
                auto i_end = tmp_hist_name.rfind('_');
                if(i_end - i_start > 0)
                {
                    std::string tmp_sample_name = tmp_hist_name.substr(i_start, i_end - i_start);

                    // set branching ratio fraction
                    if(tmp_sample_name == std::string("tl208_int_rot") ||
                       tmp_sample_name == std::string("tl208_feShield") ||
                       tmp_sample_name == std::string("tl208_pmt"))
                    {
                        //std::cout << "tmp_sample_name=" << tmp_sample_name << std::endl;
                        //std::cin.get();
                        activity_scale_branching_ratio = 0.36;
                    }

                    //std::cout << "tmp_sample_name=" << tmp_sample_name << std::endl;
                    if(MCNameToParamNumberMap.count(tmp_sample_name) > 0)
                    {
                        int paramNumber = MCNameToParamNumberMap.at(tmp_sample_name);
                        // TODO: removed std::string, change tmpName type to be std::string from TString
                    
                        //std::cout << "paramNumber=" << paramNumber << " -> tmp_sample_name=" << tmp_sample_name << " ~> tmpName=" << tmpName << std::endl;                    
                        //which_param = paramNumber;
                        which_param = paramNumberToMinuitParamNumberMap.at(paramNumber);
                        found_param = true;
                    }
                    else
                    {
                       std::cout << "ERROR: could not find " << tmp_sample_name << " in MCNameToParamNumberMap" << std::endl;
                    }
                }
            }
            /*
            }
            catch(std::exception &e)
            {
                std::cout << "e.what(): " << e.what() << std::endl;
                std::cout << "tmpName=" << tmpName << std::endl;
                std::cout << "contents of map" << std::endl;
                for(auto it = MCNameToParamNumberMap.cbegin(); it != MCNameToParamNumberMap.cend(); ++ it)
                {
                    std::cout << it->first << " -> " << it->second << std::endl;
                }
            }
            std::cin.get();
            */



            if(found_param == true)
            {
                //std::cout << "found histogram: tmpName=" << tmpName << " which_param=" << which_param << std::endl;

                // scale histogram to correct size using output parameter
                // from fit
                if(which_param >= numberEnabledParams)
                {
                    std::cout << "throwing exception, which_param=" << which_param << std::endl;
                    throw std::runtime_error("which_param invalid value");
                }


            /* if I rebuild the MC in the fitHistograms loop
             * then I don't need to reweight it here
                // check if MC sample is 150Nd, if so need to apply xi
                // reweighting
                TString tmpHist_name = tmpHist_draw1D->GetName();
                if(tmpHist_name.CompareTo("nd150_rot_2n2b_m4") == 0)
                {
                    std::cout << "found the 150Nd MC" << std::endl;
                    std::cin.get();

                    //TH1F *tmpHist_draw1D_clone = nullptr;
                    TH1F *tmpHist_reweight = nullptr;
                    //reweight_apply(tmpHist_draw1D_clone, tmpHist_draw1D, ... );
                    reweight_apply(tmpHist_reweight, "Nd150_2eNg_output_truth_postprocessed.root", xi_31, xi_31_baseline, h_nEqNull, h_nEqTwo, psiN0, psiN2, bb_Q);
                    // TODO: after reweight function called, replace 150nd MC
                    // in containers, or add a _reweight version to containers

                    tmpHist_drawpointer = tmpHist_draw1D_clone;
                }
                else
                {
                    std::cout << "it is not " << tmpHist_draw1D->GetName() << std::endl;
                    std::cin.get();
                }
              */

                // no error thrown, which_param is presumably the correct index
                Double_t activity_scale = AdjustActs[which_param] * activity_scale_branching_ratio;
                //if()
                //{
                //    activity_scale *= 0.36;
                //}
                tmpHist_draw1D->Scale(activity_scale);
                // TODO: fix this

                //TH1F *tmpHist_drawpointer = tmpHist_draw1D;

                
                if(tmpHist_draw1D->Integral() > 0)
                //if(tmpHist_drawpointer->Integral() > 0)
                {
                    stacks1D[i]->Add(tmpHist_draw1D);
                    //stacks1D[i]->Add(tmpHist_drawpointer);

                    if(j == 0)
                    {
                        hAllMC1D[i] = (TH1F*)tmpHist_draw1D->Clone("Total MC");
                        //hAllMC1D[i] = (TH1F*)tmpHist_drawpointer->Clone("Total MC");
                    }
                    else
                    {
                        hAllMC1D[i]->Add((TH1F*)tmpHist_draw1D);
                        //hAllMC1D[i]->Add((TH1F*)tmpHist_drawpointer);
                    }
	            }
                else
                {
                    std::cout << "not adding to stack, Integral() <= 0: " << tmpHist_draw1D->GetName() << std::endl;
                }
            }
            else
            {
                std::cout << "error could not find histogram: tmpName=" << tmpName << std::endl;
            } 

        }


        stacks1D[i]->SetMaximum(350.);
        stacks1D[i]->Draw("hist");
        hAllMC1D[i]->SetLineWidth(2);
        hAllMC1D[i]->SetLineColor(kBlack);
        hAllMC1D[i]->SetFillColor(kWhite);
        hAllMC1D[i]->SetFillStyle(0);
        hAllMC1D[i]->Sumw2();
        //hAllMC1D[i]->Draw("hist sames");
        TString Nmc_str;
        Nmc_str.Form("%i", (int)hAllMC1D[i]->Integral()); // TODO: float?
        hAllMC1D[i]->SetTitle("Total MC (" + Nmc_str + ")");
        hAllMC1D[i]->Draw("hist same");
        data1D[i]->SetLineWidth(2);
        data1D[i]->SetMarkerStyle(20);
        data1D[i]->SetMarkerSize(0.5);
        TString Ndata_str;
        Ndata_str.Form("%i", (int)data1D[i]->Integral()); // TODO: float?
        data1D[i]->SetTitle("Data (" + Ndata_str + ")");
        //data1D[i]->Draw("PEsames");
        data1D[i]->Draw("PEsame");

        double chi2;
        int ndf;
        int igood;
        TString chi2_str;
        TString ndf_str;

        // TODO: should chisquare value include the constraints? because at
        // the moment it does not

        // TODO: chi2 value is different from fit_2e code
        double prob = data1D[i]->Chi2TestX(hAllMC1D[i], chi2, ndf, igood, "UW");
        // TODO: check if I can get fcn value from the minuit fit object
        chi2_str.Form("%4.3f", chi2);
        ndf_str.Form("%i", ndf);
        TLegend *leg;
        leg = c->BuildLegend();
        leg->SetName("leg" + i_str + "_");
        leg->SetFillColor(kWhite);
        leg->AddEntry((TObject*)0, "#chi^{2}/ndf = " + chi2_str + "/" + ndf_str, "");
        leg->Draw();

        //std::cout << "saving to file (canvas)" << std::endl;
        //c->SaveAs("finalHisto1D_" + i_str + ".C");
        //c->SaveAs("finalHisto1D_" + i_str + ".eps");
        //c->SaveAs("finalHisto1D_" + i_str + ".png");
    
    }


    ///////////////////////////////////////////////////////////////////////////
    // draw covarience matrix
    ///////////////////////////////////////////////////////////////////////////

    std::cout << "draw: covariance matrix" << std::endl;

    /*
    std::cout << "contents of free_params:" << std::endl;
    for(int i = 0; i < free_params.size(); ++ i)
    {
        std::cout << free_params.at(i) << std::endl;
    }
    */

    // do not store names elsewhere, create them here
    TH2D* hCorrMatrix = new TH2D("hCorrMatrix", "Correlation Matrix",
                                 number_free_params, 0, (double)number_free_params,
                                 number_free_params, 0, (double)number_free_params);

    std::cout << "number_free_params=" << number_free_params << std::endl;
    for(int i = 0; i < number_free_params; i++)
    {
        for(int j = 0; j < number_free_params; j++)
        {
            //std::cout << index_free_params.at(i) << "  " << index_free_params.at(j)<< std::endl;
            //CovMatrix[i][j] = minuit->GetCovarianceMatrixElement(i,j);

            // std::cout << CovMatrix[i][j] << std::endl;
            //Double_t value = CovMatrix[i][j] / (std::sqrt(CovMatrix[i][i]) * std::sqrt(CovMatrix[j][j]));
            //
            Double_t CovMatrix_i_j = CovMatrix[i * number_free_params + j];
            Double_t CovMatrix_i_i = CovMatrix[i * number_free_params + i];
            Double_t CovMatrix_j_j = CovMatrix[j * number_free_params + j];

            Double_t value = CovMatrix_i_j / (std::sqrt(CovMatrix_i_i) * std::sqrt(CovMatrix_j_j));
            //hCorrMatrix->Fill(free_params_names.at(i), free_params_names.at(j), value);
            //hCorrMatrix->Fill(i, j, value);
            
            // need to ignore all HARD constrained (minuit->FixParameter)
            // parameters
            // I don't know if I like this solution, might be better to
            // iterate over the vectors themselves
            // TODO: what order is this vector in and how does that affect
            // output?
            // TODO: change to another incremental index, and have an
            // if(find(fixed params)) statement
            //std::cout << "i=" << i << " j=" << j << std::endl;

            int i_free_param_index = free_params.at(i);
            int j_free_param_index = free_params.at(j);
            //std::cout << "i_free_param_index=" << i_free_param_index << " j_free_param_index=" << j_free_param_index << std::endl;

            // convert minuit/internal parameter
            // ("free parameter" meaning not SOFT or FREE constrain mode)
            // to external parameter number
            //int i_external = minuitParamNumberToParamNumberMap.at(i_free_param_index);
            //int j_external = minuitParamNumberToParamNumberMap.at(j_free_param_index);
            //std::cout << "i_external=" << i_external << " j_external=" << j_external << std::endl;
            // NOTE: using method of accessing vector elements in free_params
            // gives index which is already in external index format
            int i_external = i_free_param_index;
            int j_external = j_free_param_index;

            // convert external parameter number to parameter name
            TString free_param_name_i = paramNameMap[i_external];
            TString free_param_name_j = paramNameMap[j_external];
            //std::cout << "free_param_name_i=" << free_param_name_i << " free_param_name_j=" << free_param_name_j << std::endl;

            // convert parameter name to human readable name
            //TString free_param_human_name_i = paramNameToHumanParamNameMap.at(free_param_name_i);
            //TString free_param_human_name_j = paramNameToHumanParamNameMap.at(free_param_name_j);
            //TString free_param_human_name_i = paramNameToHumanReadableParamNameMap.at(free_param_name_i);
            //TString free_param_human_name_j = paramNameToHumanReadableParamNameMap.at(free_param_name_j);
            //std::cout << "free_param_human_name_i=" << free_param_human_name_i << " free_param_human_name_j=" << free_param_human_name_j << std::endl;

            TString free_param_human_name_i = "NONAME";
            TString free_param_human_name_j = "NONAME";
            // if available, use name from manually assigned names in parameter_names.lst file
            if(paramNumberToHumanReadableParamNameMap.count(i_external) == 1)
            {
                free_param_human_name_i = paramNumberToHumanReadableParamNameMap.at(i_external);
            }
            else
            {
                free_param_human_name_i = paramNameToHumanReadableParamNameMap.at(free_param_name_i);
            }
            if(paramNumberToHumanReadableParamNameMap.count(j_external) == 1)
            {
                free_param_human_name_j = paramNumberToHumanReadableParamNameMap.at(j_external);
            }
            else
            {
                free_param_human_name_j = paramNameToHumanReadableParamNameMap.at(free_param_name_j);
            }

            // fill
            hCorrMatrix->Fill(free_param_human_name_i, free_param_human_name_j, value);

        }
    }
  
    // TODO: re-enable this check chi2 value is the same
    // remove global variable for chi2?
    // these were commented out for some reason, I did not remove them
    /*
    double chi2, edm, errdef; 
    int nvpar, nparx;
    minuit->GetStats(chi2, edm, errdef, nvpar, nparx);
    int ndf = npfits-nvpar;
    */


    // TODO: margin
    //c->SetLeftMargin(0.15)
    //c->SetTopMargin(0.15);
    
    // create a function to draw the final stacked histograms
    gStyle->SetPaintTextFormat("4.2f");
    gStyle->SetGridStyle(0);

    hCorrMatrix->SetTitle("Parameter Correlation Matrix");
    hCorrMatrix->LabelsDeflate();
    hCorrMatrix->GetXaxis()->LabelsOption("v");
    hCorrMatrix->GetXaxis()->SetLabelSize(0.035);
    hCorrMatrix->GetYaxis()->SetLabelSize(0.035);
    hCorrMatrix->GetZaxis()->SetRangeUser(-1.0, 1.0);

    //TCanvas *c2 = new TCanvas("c2");
    TCanvas *c2 = new TCanvas("c2", "Parameter Correlation Matrix", 200, 10, 1013, 885);
    c2->SetLeftMargin(0.18);
    c2->SetRightMargin(0.14);
    c2->SetTopMargin(0.14);
    c2->SetBottomMargin(0.21);
    c2->SetFillColor(kWhite);
    c2->SetGrid();
    hCorrMatrix->Draw("colz text");
   
    // clean
    delete [] CovMatrix;

    // TODO: save plots

    ///////////////////////////////////////////////////////////////////////////
    // draw results - 2d
    ///////////////////////////////////////////////////////////////////////////

#if 0
    // each channel 1D hists
    THStack *stacks2D[number2DHists];
    TH2F *data2D[number2DHists];
    for(int i = 0; i < allDataSamples2D->GetEntries(); i++)
    {
        data2D[i] = (TH2F*)allDataSamples2D->At(i)->Clone();
        TString i_str;
        i_str.Form("%i",i);
        c = new TCanvas("c_"+i_str);
        c->SetFillColor(kWhite);
        stacks2D[i] = new THStack("stacks2D"+i_str,i_str);

        // j list MC samples for this channel i
        TH2F *tmpHist_draw2D;
        for(int j = 0; j < allMCSamples2D[i]->GetEntries(); j++)
        {
            TString j_str;
            j_str.Form("%i",j);
            tmpHist_draw2D = (TH2F*)allMCSamples2D[i]->At(j)->Clone();
            TString tmpName = tmpHist_draw2D->GetName();

            //std::cout << "looking for " << tmpName << std::endl;
            
            // k searching through array of params for the right one
            int which_param;
            bool foundParam = false;
            for(int k = 0; (k < numberParams) && !foundParam; k++)
            {
            
                // match up the isotope with parm
                for(int n = 0; (n < paramNameMap[k].size()) && !foundParam; n++)
                {
	                //std::cout <<"is it...  " <<paramNameMap[k].at(j) << std::endl;

	                if(tmpName.Contains(paramNameMap[k].at(n)))
                    {
	                    foundParam = true;
                        which_param = k;
	                }
	            }
            }

            if(foundParam)
            {
                tmpHist_draw2D->Scale(AdjustActs[which_param]);
     
                if(tmpHist_draw2D->Integral() > 0)
                {
                    stacks2D[i]->Add(tmpHist_draw2D);
                }
            }
            else
            {
                std::cout << "error could not find histograms" << std::endl;
            }

        }

        stacks2D[i]->Draw("");
        data2D[i]->SetLineWidth(2);
        data2D[i]->SetMarkerStyle(20);
        data2D[i]->SetMarkerSize(0.5);
        data2D[i]->Draw("PEsames");

        c->SaveAs("finalHisto2D_" + i_str + ".C");
        c->SaveAs("finalHisto2D_" + i_str + ".eps");
        c->SaveAs("finalHisto2D_" + i_str + ".png");
        std::cout << "saved 2D histogram to file" << std::endl;

    }
#endif


}


///////////////////////////////////////////////////////////////////////////////
// book1DHistograms
///////////////////////////////////////////////////////////////////////////////

void book1DHistograms_helper(TFile *myFile, Int_t channel_counter, TString theChannel, TString thePhase_arg, TString theHistogram,
    const int nBkgs, TString *BkgFiles)//, TH1F *tmpHist)
{
        
    TH1F *tmpHist = nullptr;


    for(int i = 0; i < nBkgs; i++)
    {
        // check if parameter is enabled
        // convert parameter string name to index
        //std::cout << "searching for string " << ExternalBkgOneParamFiles[i] << " in paramNameMap" << std::endl;
        //TString search_object = ExternalBkgOneParamFiles[i];
        std::string mc_name = std::string(BkgFiles[i]);
        std::string search_object = MCNameToParamNameMap.at(mc_name);
        if(paramNameToNumberMap.count(search_object) > 0)
        {
            // convert from mc sample name to param number

            int param_number = paramNameToNumberMap.at(search_object);
            //std::cout << "parameber number " << param_number << " is in the paramNameToNumberMap" << std::endl;
            if(std::find(enabled_params.begin(), enabled_params.end(), param_number) != enabled_params.end())
            {
                // check if param number is enabled

                std::string directory("scaled/hTotalE_/");
                std::string name(theHistogram + BkgFiles[i] + "_fit_scaled");
                std::string fullname = directory + name;
                std::string new_name(theHistogram + BkgFiles[i] + "_fit");
                std::cout << "fullname=" << fullname << std::endl;

                //gDirectory->GetListOfKeys();

                //tmpHist = (TH1F*)gDirectory->Get(fullname.c_str())->Clone();
                tmpHist = (TH1F*)myFile->Get(fullname.c_str())->Clone(new_name.c_str());

                if(tmpHist != nullptr)
                //if(gDirectory->GetListOfKeys()->Contains(fullname.c_str()))
                //std::cout << "parameter number " << param_number << " is enabled" << std::endl;
                //std::string name(theHistogram + BkgFiles[i] + "_fit");
                //if(gDirectory->GetListOfKeys()->Contains(name.c_str()))
                {
                    // load sample

                    // 2020-04-03: removed changing of histogram name
                    //check if the histograms exists 
                    //std::string hist_name(BkgFiles[i] + "_" + theChannel + thePhase_arg);
                    //std::cout << "Get() : " << name << " from file, Clone() : " << hist_name << std::endl;
                    //tmpHist = (TH1F*)gDirectory->Get(name.c_str())->Clone(hist_name.c_str());
                    //tmpHist = (TH1F*)gDirectory->Get(fullname.c_str())->Clone();

                    // scale by activity

                    // convert parameter number to minuit parameter number
                    //minuit_param_number = paramNumberToMinuitParamNumberMap.at(param_number);

                    // TODO: change such that samples are pre-scaled by activity input value
                    // get initial parameter values and error
                    Double_t param_init_value = 0.;
                    Double_t param_init_error = 0.; 
                    if(thePhase == 0)
                    {
                        int i = param_number;
                        param_init_value = paramInitValueP1Map[i];
                        param_init_error = paramInitErrorP1Map[i];
                    }
                    else if(thePhase == 1)
                    {
                        int i = param_number;
                        param_init_value = paramInitValueP2Map[i];
                        param_init_error = paramInitErrorP2Map[i];
                    }
                    else
                    {
                        std::cout << "ERROR: Invalid value for thePhase: thePhase=" << thePhase << std::endl;
                    }
                    Double_t scale_factor = param_init_value;

                    // account for 208 Tl branching ratio of 36 %
                    if(mc_name == std::string("tl208_int_rot") ||
                       mc_name == std::string("tl208_feShield") || // TODO: this doesn't seem to work
                       mc_name == std::string("tl208_pmt"))
                       // TODO: do not apply to tl208_air ?
                    {
                        //std::cout << "mc_name=" << mc_name << " applying additional scaling factor of 0.36" << std::endl;
                        //std::cin.get();
                        scale_factor *= 0.36;
                    }

                    // NOTE: TODO
                    // possible flaw with this method: error is no longer
                    // pre-set using values from input file
                    // TODO: note this in input file documentation
                    // however, this may be an improvement because it
                    // guarantees minuit is responsible for error estimation
                    tmpHist->Scale(scale_factor);
                    // samples are now scaled by activity
                    // changed input, and pre-scaling, now need to change output


// NOTE: do NOT apply xi reweighting here
// this section just LOADS histograms from file and we want to LOAD
// the default (not reweighted) nd150 spectra
// TODO: this may no longer be true
/*
        // MARKER 10
                // check if MC sample is 150Nd, if so need to apply xi
                // reweighting
                if(tmpHist_draw1D->GetName().CompareTo("nd150_rot_2n2b_m4") == 0)
                {
                    std::cout << "found the 150Nd MC" << std::endl;
                    std::cin.get();

                    TH1F *tmpHist_draw_1d_clone = nullptr;
                    reweight_apply(tmpHist_draw1D_clone, tmpHist_draw1D, ... );
                }
                else
                {
                    std::cout << "it is not " << tmpHist_draw1D->GetName() << std::endl;
                    std::cin.get();
                }
*/


                    allMCSamples1D[channel_counter]->Add(tmpHist);
                    // TODO: does this work as expected for secular equlibrium samples?

                    //std::cout << tmpHist->GetName() << std::endl;

                }
                else
                {
                    std::cout << "gDirectory->GetListOfKeys() does not contain " << fullname << " - disabling parameter number " << param_number << std::endl;
                    // cannot find histogram input data, so disable parameter
                    std::remove(enabled_params.begin(), enabled_params.end(), param_number);
                }
            }
            else
            {
                // paramter not enabled, do not load histogram/sample
                std::cout << "parameter number " << param_number << " is not enabled (not found in vector)" << std::endl;
            }
        }
        else
        {
            std::cout << "!!!!! ERROR: search_object=" << search_object << " not found in paramNameToNumberMap" << std::endl;
            std::cout << "mc_name=" << mc_name << std::endl;

            std::cout << "contents of map paramNameToNumberMap:" << std::endl;
            for(auto it = paramNameToNumberMap.cbegin(); it != paramNameToNumberMap.cend(); ++ it)
            {
                std::cout << it->first << " -> " << it->second << std::endl;
            }
            std::cout << "contents of map MCNameToParamNameMap:" << std::endl;
            for(auto it = MCNameToParamNameMap.cbegin(); it != MCNameToParamNameMap.cend(); ++ it)
            {
                std::cout << it->first << " -> " << it->second << std::endl;
            }
        }
    }
}

// channel_counter = 0
// theChannel = "2e_"
// thePhase = "P1"
// theHistogram = "hTotalE_"
void book1DHistograms(Int_t channel_counter, TString theChannel, TString thePhase_arg, TString theHistogram) {

    std::cout << "booking 1D hists for " << theChannel << " " << thePhase_arg << std::endl;
    allMCSamples1D[channel_counter] = new TObjArray();

    //TFile *aFile = TFile::Open("/home/ebirdsall/NEMO3/Nd150_analysis/MeasureStuff/new/Macros/Nd150_" + theChannel + thePhase_arg + ".root");
    TFile *aFile = TFile::Open("Nd150_" + theChannel + thePhase_arg + ".root");
    //gDirectory->cd("singleHistos");
    //gDirectory->ls();


    //TH1F *tmpHist = nullptr; //new TH1F("tmpHist_" + theChannel + thePhase_arg, "" , 1, 0, 1);
    
    std::cout << "External" << std::endl;
    book1DHistograms_helper(aFile, channel_counter, theChannel,
                            thePhase_arg, theHistogram,
                            nExternalBkgs,
                            ExternalBkgFiles);//,
                            //tmpHist);

    // TODO: does this work as expected for secular equlibrium samples?

    std::cout << "Internal" << std::endl;
    book1DHistograms_helper(aFile, channel_counter, theChannel,
                            thePhase_arg, theHistogram,
                            nInternalBkgs,
                            InternalBkgFiles);//,
                            //tmpHist);

    std::cout << "Rn 222" << std::endl;
    book1DHistograms_helper(aFile, channel_counter, theChannel,
                            thePhase_arg, theHistogram,
                            nRn222Bkgs,
                            Rn222BkgFiles);//,
                            //tmpHist);

    std::cout << "Rn 220" << std::endl;
    book1DHistograms_helper(aFile, channel_counter, theChannel,
                            thePhase_arg, theHistogram,
                            nRn220Bkgs,
                            Rn220BkgFiles);//,
                            //tmpHist);

    std::cout << "Neighbour" << std::endl;
    book1DHistograms_helper(aFile, channel_counter, theChannel,
                            thePhase_arg, theHistogram,
                            nNeighbours,
                            NeighbourFiles);//,
                            //tmpHist);

    std::cout << "Nd150" << std::endl;
    book1DHistograms_helper(aFile, channel_counter, theChannel,
                            thePhase_arg, theHistogram,
                            nNd150Samples,
                            Nd150Files);//,
                            //tmpHist);

    // TODO here
    // what is name in other section of code
    //std::string name(theHistogram + "data_2e");
    std::string directory("processeddata/hTotalE_/");
    std::string name(theHistogram + "data_2e");
    std::string fullname = directory + name;
    std::cout << "fullname=" << fullname << std::endl;
    //if(gDirectory->GetListOfKeys()->Contains(fullname.c_str()))
    //TH1F *tmpHist = (TH1F*)gDirectory->Get(fullname.c_str())->Clone();
    TH1F *tmpHist = (TH1F*)aFile->Get(fullname.c_str())->Clone();
    if(tmpHist != nullptr)
    {
        //TH1F *tmpHist = nullptr;
        // 2020-04-03: removed changing of histogram name
        //std::string hist_name("data_" + theChannel + thePhase_arg);
        //std::cout << "Get() : " << name << " from file, Clone() : " << hist_name << std::endl;
        //tmpHist = (TH1F*)gDirectory->Get(name.c_str())->Clone(hist_name.c_str());
        //tmpHist = (TH1F*)gDirectory->Get(fullname.c_str())->Clone();
        allDataSamples1D->Add((TH1F*)tmpHist);
    }
    else
    {
        std::cout << "gDirectory->GetListOfKeys() does not contain " << fullname << std::endl;
    }
    /*
    if(gDirectory->GetListOfKeys()->Contains(theHistogram + "Data"))
    {
        std::string name(theHistogram + "Data");
        std::cout << "Get() : " << name << " from file, Clone() : " << "Data_" + theChannel + thePhase_arg << std::endl;
        tmpHist = (TH1F*)gDirectory->Get(name.c_str())->Clone("Data_" + theChannel + thePhase_arg);
        allDataSamples1D->Add((TH1F*)tmpHist);
    }
    else
    {
        std::cout << "gDirectory->GetListOfKeys() does not contain " << theHistogram + "Data" << std::endl;
    }
    */

    // std::cout << tmpHist->GetName() << std::endl;
    // tmpHist->Delete();
    // aFile->Close();
    // aFile->Delete();
}


///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////

void book2DHistograms(Int_t channel_counter, TString theChannel, TString thePhase, TString theHistogram) {

  std::cout << "booking 2D hists for "<< theChannel<<" " <<thePhase<<std::endl;
  allMCSamples2D[channel_counter] = new TObjArray();

  TFile *aFile = TFile::Open("/home/blotsd/NEMO3/Nd150_analysis/MeasureStuff/Macros/Nd150_"+theChannel+thePhase+".root");
  gDirectory->cd("singleHistos");
  // gDirectory->ls();
  TH2F *tmpHist2  = new TH2F("tmpHist2_"+theChannel+thePhase,"",1,0,1,1,0,1);
  for ( int i = 0; i < nExternalBkgs; i++ ) {
    if ( gDirectory->GetListOfKeys()->Contains(theHistogram+ExternalBkgFiles[i]+"_fit") ) {//check if the histograms exists
      tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+ExternalBkgFiles[i]+"_fit")->Clone(ExternalBkgFiles[i]+"_"+theChannel+thePhase);
      allMCSamples2D[channel_counter]->Add(tmpHist2);

      //std::cout << tmpHist->GetName() << std::endl;

    }
  }
  /*
  for ( int i = 0; i < nExternalBkgs_twoParam; i++ ) {
    if ( gDirectory->GetListOfKeys()->Contains(theHistogram+ExternalBkgTwoParamFiles[i]+"_fit") ) {//check if the histograms exists
      tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+ExternalBkgTwoParamFiles[i]+"_fit")->Clone(ExternalBkgTwoParamFiles[i]+"_"+thePhase);
      allMCSamples2D[channel_counter]->Add(tmpHist2);

      //std::cout << tmpHist->GetName() << std::endl;

    }
  }*/

  for ( int i = 0; i < nInternalBkgs; i++ ) {
    if ( gDirectory->GetListOfKeys()->Contains(theHistogram+InternalBkgFiles[i]+"_fit") ) {//check if the histograms exists
      tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+InternalBkgFiles[i]+"_fit")->Clone(InternalBkgFiles[i]+"_"+theChannel+thePhase);
      allMCSamples2D[channel_counter]->Add(tmpHist2);
    }
  }

  for ( int i = 0; i < nRn222Bkgs; i++ ) {
    if ( gDirectory->GetListOfKeys()->Contains(theHistogram+Rn222BkgFiles[i]+"_fit") ) {//check if the histograms exists
      tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+Rn222BkgFiles[i]+"_fit")->Clone(Rn222BkgFiles[i]+"_"+theChannel+thePhase);
      allMCSamples2D[channel_counter]->Add(tmpHist2);
    }
  }

  for ( int i = 0; i < nRn220Bkgs; i++ ) {
    if ( gDirectory->GetListOfKeys()->Contains(theHistogram+Rn220BkgFiles[i]+"_fit") ) {//check if the histograms exists
      tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+Rn220BkgFiles[i]+"_fit")->Clone(Rn220BkgFiles[i]+"_"+theChannel+thePhase);
      allMCSamples2D[channel_counter]->Add(tmpHist2);
    }
  }

  for ( int i = 0; i < nNd150Samples; i++ ) {
    if ( gDirectory->GetListOfKeys()->Contains(theHistogram+Nd150Files[i]+"_fit") ) {//check if the histograms exists
      tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+Nd150Files[i]+"_fit")->Clone(Nd150Files[i]+"_"+theChannel+thePhase);
      allMCSamples2D[channel_counter]->Add(tmpHist2);
    }
  }

  for ( int i = 0; i < nNeighbours; i++ ) {
    if ( gDirectory->GetListOfKeys()->Contains(theHistogram+NeighbourFiles[i]+"_fit") ) {//check if the histograms exists
      tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+NeighbourFiles[i]+"_fit")->Clone(NeighbourFiles[i]+"_"+theChannel+thePhase);
      allMCSamples2D[channel_counter]->Add(tmpHist2);
    }
  }

  if ( gDirectory->GetListOfKeys()->Contains(theHistogram+"Data") ) {
    tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+"Data")->Clone("Data_"+theChannel+thePhase);
    allDataSamples2D->Add((TH2F*)tmpHist2);
  }
  // NOTE: Data changed to data but not updated here

  std::cout << tmpHist2->GetName() << std::endl;
  // tmpHist2->Delete();
  //aFile->Close();
  // aFile->Delete();
}


///////////////////////////////////////////////////////////////////////////////
// fitBackgrounds
///////////////////////////////////////////////////////////////////////////////

//void fitBackgrounds(double *AdjustActs, double *AdjustActs_Err, double *&CovMatrix, int& number_free_params, Int_t thePhase)
TMinuit * fitBackgrounds(double *AdjustActs, double *AdjustActs_Err, double *&CovMatrix, int& number_free_params, Int_t thePhase)
{

    
    

    std::cout << ">>>>> fitBackgrounds()" << std::endl;

    //TVirtualFitter::SetDefaultFitter("Minuit2");
    
    // TODO:
    // want to change Minuit so that it does not store disabled parameters
    // do this using a map to translate from parameter number to minuit
    // parameter numbers
    // and change numberParams to numberEnabled params
    //TMinuit *minuit = new TMinuit(numberParams);
    std::cout << "numberEnabledParams=" << numberEnabledParams << std::endl;
    TMinuit *minuit = new TMinuit(numberEnabledParams);
    // TODO working here need to check all instances of numberParams


    std::cout << "Fit created" << std::endl;

    //std::cout << "minuit tests" << std::endl;
    //minuit->DefineParameter(1, "test1", 100.0, 10.0, -500.0, 500.0);
    //std::cin.get();

    // set debug level
    minuit->SetPrintLevel(1);

    // moved reading of parameter list file to before book histogram function
    // calls
  
    TString sample_names[numberParams];
    // read samples names out again just to make sure its correct (useful when changing the number of params)
    for(int i = 0; i < numberParams; i++)
    {
        // tmpStr doesn't appear to do anything ?
        // TODO
        TString tmpStr = "";

        TString i_str;
        i_str.Form("%i", i);

        sample_names[i] = i_str;
        
        std::cout << "param " << i << ": sample_names[" << i << "]=" << sample_names[i] << " -> ";
        std::cout << paramNameMap[i] << " : ";
        bool first = true;
        for(auto it = MCNameToParamNameMap.cbegin(); it != MCNameToParamNameMap.cend(); ++ it)
        {
            if(it->second == paramNameMap[i])
            {
                // match for MC name
                if(!first) std::cout << ", ";
                std::cout << it->first;
                first = false;
            }
        }
        std::cout << std::endl;

        /*
        for(int j = 0; j < paramNameMap[i].size(); j++)
        { 
            std::cout << paramNameMap[i].at(j);
            tmpStr += paramNameMap[i].at(j);
            
            if(j + 1 < paramNameMap[i].size())
            {
                std::cout << ", ";
                tmpStr += ",";
            }
        }
        */
        
        std::cout << std::endl;
    }

    std::cout << "set errors" << std::endl;

    // Set the parameters for the fit, and give them an arbitrary 10% error to start.
    //for(int i = 0; i < numberParams; i++)
    //int i_minuit = 0; // internal minuit counter
    for(int i = 0; i < numberParams; i++)
    {
        // internal (minuit) parameter number
        int minuit_param_number = -1;

        // check if parameter enabled
        if(std::find(enabled_params.begin(), enabled_params.end(), i) == enabled_params.end())
        {
            // NOT enabled
            // ignore
            continue;
        }
        else
        {
            // is enabled
            // do nothing (exec code in following block)

            // set internal parameter number
            minuit_param_number = paramNumberToMinuitParamNumberMap.at(i);

            /*
            std::cout << "contents of paramNumberToMinuitParamNumberMap" << std::endl;
            for(auto it = paramNumberToMinuitParamNumberMap.cbegin(); it != paramNumberToMinuitParamNumberMap.cend(); ++ it)
            {
                std::cout << it->first << " -> " << it->second << std::endl;
            }
            std::cin.get();
            */
        }
            

        // AdjustActs[i] = 1.;
        //AdjustActs_Err[i] = 1.;
        TString i_str;
        i_str.Form("%i", i);
        TString minuit_param_number_str;
        minuit_param_number_str.Form("%i", minuit_param_number);
        std::cout << "DefineParameter: i=" << i << " -> minuit_param_number=" << minuit_param_number << std::endl;
        //minuit->DefineParameter(i, "_" + i_str + "_", 1.0, 0.1, 0.0, 1000.0);

        // TODO: change such that samples are pre-scaled by activity input value
        #if 0
        // get initial parameter values and error
        Double_t param_init_value = 0.;
        Double_t param_init_error = 0.; 
        if(thePhase == 0)
        {
            /*
            for(int j = 0; j < paramInitValueP1Map[i].size(); ++ j)
            {
                param_init_value += paramInitValueP1Map[i].at(j);
            }
            for(int j = 0; j < paramInitErrorP1Map[i].size(); ++ j)
            {
                param_init_error += paramInitErrorP1Map[i].at(j);
            }
            */

            param_init_value = paramInitValueP1Map[i];
            param_init_error = paramInitErrorP1Map[i];
        }
        else if(thePhase == 1)
        {
            /*
            for(int j = 0; j < paramInitValueP2Map[i].size(); ++ j)
            {
                param_init_value += paramInitValueP2Map[i].at(j);
            }
            for(int j = 0; j < paramInitErrorP2Map[i].size(); ++ j)
            {
                param_init_error += paramInitErrorP2Map[i].at(j);
            }
            */

            param_init_value = paramInitValueP2Map[i];
            param_init_error = paramInitErrorP2Map[i];
        }
        else
        {
            std::cout << "ERROR: Invalid value for thePhase: thePhase=" << thePhase << std::endl;
        }
        #endif


        if(std::find(fixed_params.begin(), fixed_params.end(), i) != fixed_params.end())
        {
            // define parameter using constrained value if hard constrained

            std::cout << "minuit: fixed parameter i=" << i << std::endl;
            
            //minuit->DefineParameter(i, "_" + i_str + "_", 1.0, 0.1, 0.0, 2.0);
            //minuit->DefineParameter(i, "_" + i_str + "_FIXED", param_init_value, param_init_error, -2.0, 1.0e+06);
            //minuit->FixParameter(i);
            
            //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_FIXED", param_init_value, param_init_error);
            // TODO: change such that samples are pre-scaled by activity input value
            if(i == 1)
            {
                minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_FIXED", AdjustActs[i], 0.5, 0.0, 50.0);
            }
            else
            {
                minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_FIXED", 1.0, 0.5, 0.0, 50.0);
            }
            minuit->FixParameter(minuit_param_number);
        }
        else
        {
            // define parameter using initial value if free/soft constrained
            
            std::cout << "minuit: parameter i=" << i << " is enabled and not fixed, leaving free" << std::endl;

            //minuit->DefineParameter(i, "_" + i_str + "_", param_init_value, param_init_error, 0.0, 1.0e+06);
            
            //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", param_init_value, param_init_error, 0.0, 1.0e+05);
            // TODO: change such that samples are pre-scaled by activity input value
            //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", 1.0, 0.5, 0.0, 50.0);
            //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", 1.0, 0.5, 0.0, 1000.0);
            if(i == 1)
            {
                // TODO: fix this
                // does not work if xi_31 paramter is not number 1
                //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", xi_31_init, 0.5, 0.0, 1000.0);
                minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", AdjustActs[i], 0.5, 0.0, 1000.0);
            }
            else
            {
                minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", 1.0, 0.5, 0.0, 1000.0);
            }
            // TODO: set initial error using initError/initValue
            // TODO: limits were set to 50. minuit trying to exceed 50 for some
            // backgrounds when only using nd, mo bb, zr bb, ca bb an Ca Y90
            // allow larger range of values, check fit result
            // then put back and include more backgrounds, does fit result
            // still try and exceed 50x
            // TODO: there is still the issue of the discontinuity in
            // the chisquare plots. MPS. is this due to Poisson function
            // failing with negative events?
            
        }
        // TODO: I have changed the input, now I have to change the output
        // and also the initial scaling given to histograms

        // note, no longer need this section, because only enabled parameters
        // get through to this execution block
        // if parameter is enabled, then it is either fixed or free
        // it has been copied and edited above
    /*
        // TODO: not sure if this is the correct thing to do
        // I think that disabled params do not necessarily appear in the
        // fixed params list, and therefore I have to fix them (not sure
        // about that) here if they are disabled
        if(std::find(enabled_params.begin(), enabled_params.end(), i) == enabled_params.end())
        {
            // NOT enabled, therefore fix
            
            std::cout << "minuit: disabled parameter i=" << i << std::endl;
         
            // leave parameter (amplitude) at default value of 1.0, which
            // is more likely to show up any errors in code
            minuit->DefineParameter(i, "_" + i_str + "_DISABLED", 1.0, 0.1, 0.0, 2.0);
            //minuit->DefineParameter(i, "_" + i_str + "_", param_init_value, param_init_error, 0.0, 2.0);
            minuit->FixParameter(i);

            // TODO: this code never executed, because disabled parameters are
            // always fixed? or always set as "hard" in parameter list file?
        }
        else if(std::find(fixed_params.begin(), fixed_params.end(), i) != fixed_params.end())
        {
            // define parameter using constrained value if hard constrained

            std::cout << "minuit: fixed parameter i=" << i << std::endl;
            
            //minuit->DefineParameter(i, "_" + i_str + "_", 1.0, 0.1, 0.0, 2.0);
            minuit->DefineParameter(i, "_" + i_str + "_FIXED", param_init_value, param_init_error, -2.0, 50000.0);
            minuit->FixParameter(i);
        }
        else
        {
            // define parameter using initial value if free/soft constrained
            
            std::cout << "minuit: parameter i=" << i << " is enabled and not fixed, leaving free" << std::endl;

            minuit->DefineParameter(i, "_" + i_str + "_", param_init_value, param_init_error, 0.0, 10000.0);
            
        }
        // TODO: what about disabled parameters
    */

        //if ((i != 17 ) && (i < 32) || (i > 38)  ) {
        //if ( (i > 32) ) {
        //  minuit->FixParameter(i);
        //  fixed_params.push_back(i);
        //}

        //if ( i == 24 ) {
        //  minuit->FixParameter(i);
        //  fixed_params.push_back(i);
        //}
    }

    std::cout << "all parameters fixed" << std::endl;

    // TODO: check what parameters Summer had fixed and check the MC cuts
    // Fix the tracker radon activities
    //minuit->FixParameter(34);
    //minuit->FixParameter(35);
    //minuit->FixParameter(23);
    // minuit->FixParameter(24);
    /*
    minuit->FixParameter(1);
    minuit->FixParameter(3);
    minuit->FixParameter(7);
    minuit->FixParameter(11);
    minuit->FixParameter(21);
    minuit->FixParameter(22);
    minuit->FixParameter(29);
    minuit->FixParameter(30);
    //minuit->FixParameter(31);
    minuit->FixParameter(34);
    minuit->FixParameter(35);
    //minuit->FixParameter(48);
    minuit->FixParameter(49);
    minuit->FixParameter(50);
    minuit->FixParameter(51);
    minuit->FixParameter(52);
    minuit->FixParameter(53);
    minuit->FixParameter(54);

    fixed_params.push_back(1);
    fixed_params.push_back(3);
    fixed_params.push_back(7);
    fixed_params.push_back(11);
    fixed_params.push_back(21);
    fixed_params.push_back(22);
    fixed_params.push_back(29);
    fixed_params.push_back(30);
    //fixed_params.push_back(31);
    fixed_params.push_back(34);
    fixed_params.push_back(35);
    fixed_params.push_back(49);
    fixed_params.push_back(50);
    fixed_params.push_back(51);
    fixed_params.push_back(52);
    fixed_params.push_back(53);
    fixed_params.push_back(54);
    */
    //fixed_params.push_back(48);




    // Array to pass arguments to fitter
    //double arglist[100];
    //arglist[0] = 0;

    // Set print level
    //  minuit->ExecuteCommand("SET PRINT",arglist,2);
    std::cout << "lets do the fit" << std::endl;

    // Do minimisation
    // arglist[0] = 50000;  // number of function calls
    // arglist[1] = 0.1;  // tolerance

    minuit->SetErrorDef(0.5);
    // 1.0 = chisquare
    // 0.5 = negative log likelihood
    

    // mnsimp()?

    
    // MARKER
    // disable the 150 Nd gA parameter
    //minuit->FixParameter(1);




    //minuit->SetMaxIterations(50000);
    minuit->SetMaxIterations(1000);
    //give it the function to use
    minuit->SetFCN(logLikelihood);
    //std::cout << "calling: minuit->mnsimp()" << std::endl;
    //minuit->mnsimp();
    std::cout << "calling: minuit->Migrad()" << std::endl;
    minuit->Migrad();
    // don't bother calling Migrad() for now
    // TODO
    //minuit->ExecuteCommand("SIMPLEX",arglist,2);

    // Then get results
    //for(int i = 0; i < numberParams; i++)
    for(int i = 0; i < numberEnabledParams; i++)
    {
        minuit->GetParameter(i, AdjustActs[i], AdjustActs_Err[i]);
        //AdjustActs_Err[i] = minuit->GetParError(i);
    }

   
#if 0
    ///////////////////////////////////////////////////////////////////////////
    // testing
    
    // run chisquare tests

    std::cout << "running chi-square tests: " << "variable: 100Mo (10)" << std::endl;

    int n_tests = 250;
    // 100 Mo
    int mo100_99_rot_2n2b_m14_index = paramNumberToMinuitParamNumberMap.at(10);
    std::cout << "the internal index for parameter 10 is " << mo100_99_rot_2n2b_m14_index << std::endl;
    // These are in units of minuit internal parameter units
    // To convert to external parameter units, multiply by the value of the
    // external input parameter initial activity
    // Caution: For cases where the fitted parameter minimum is not at 1.0
    // the errors must be treated as upper and lower bound separatly by adding
    // them (internal param & error) to the central value fit parameter
    // external_param_error_lowerbound = (internal_param_CV - internal_param_error) * external_param_init_value
    // similar for upperbound, then subtract and / 2.0
    double test_central_value = AdjustActs[mo100_99_rot_2n2b_m14_index];
    double test_range = 10.0 * AdjustActs_Err[mo100_99_rot_2n2b_m14_index];
    // this range should hit delta sigma = 1.0 at 66 % of the width, but it
    // doesn't.
    double test_start = test_central_value - 0.5 * test_range;
    double test_end =   test_central_value + 0.5 * test_range;
    double *test_values = new double[n_tests];
    double test_step = test_range / (double)n_tests;
    std::cout << "test_central_value=" << test_central_value << "\n"
              << "test_range=" << test_range << "\n"
              << "test_start=" << test_start << "\n"
              << "test_end=" << test_end
              << std::endl;
    int n_params = minuit->GetNumPars();
    double *params = new double[n_params];
    double *param_errs = new double[n_params];
    for(int jx = 0; jx < n_params; ++ jx)
    {
        minuit->GetParameter(jx, params[jx], param_errs[jx]);
    }
    std::ofstream ofstream_testvalue("testvalue.txt");
    for(int ix = 0; ix < n_tests; ++ ix)
    {

        test_values[ix] = test_start + test_step * ix;

        // get chisquare value for test
        double fval = 0.;

        // set parameter for 100Mo
        double test_value = test_values[ix];
        params[mo100_99_rot_2n2b_m14_index] = test_value;

        std::cout << "test: ix=" << ix << ", " << "test_value=" << test_value << std::endl;

        logLikelihood(n_params, nullptr, fval, params, 0);

        ofstream_testvalue << "value=," << test_value << ",chisquare=," << fval << std::endl;

        //void logLikelihood(Int_t & /*nPar*/, Double_t* /*grad*/, Double_t &fval, Double_t *p, Int_t /*iflag */)

    }
    ofstream_testvalue.close();
    delete [] test_values;
    //delete [] params;
    //delete [] param_errs;
    

    ///////////////////////////////////////////////////////////////////////////
#endif


    ///////////////////////////////////////////////////////////////////////////

    // testing - my phase space
    // plot phase space for Nd150 and Mo100 parameters
    


    // Some people have seen further by standing on the shoulders of giants.
    // In my case, my vision has been obscured by floundering hopeless idiots
    // doing about as much as it was possible to do to inhibit my ability to
    // see anything, either simply by being incompetent, or by actively
    // doing everything possible to throw as many spanners into the works as
    // time would allow for.

#if 0
    std::vector<TCanvas*> c_mps_v;
    std::vector<TH2D*> h_mps_v;

    // loop over all combinations of parameters
    for(int param_1_ix = 0; param_1_ix < numberEnabledParams; ++ param_1_ix)
    {
        for(int param_2_ix = 0; param_2_ix < param_1_ix; ++ param_2_ix)
        {

            int param_1_ix_external = minuitParamNumberToParamNumberMap.at(param_1_ix);
            int param_2_ix_external = minuitParamNumberToParamNumberMap.at(param_2_ix);
            
            TString param_1_ix_str_external;
            param_1_ix_str_external.Form("%i", param_1_ix_external);
            TString param_2_ix_str_external;
            param_2_ix_str_external.Form("%i", param_2_ix_external);

            TString c_mps_name_base = "c_mps";
            TString c_mps_name = c_mps_name_base + "_" + param_1_ix_str_external + "_" + param_2_ix_str_external;
            
            std::cout << "rendering: " << c_mps_name << std::endl;

            TCanvas *c_mps = new TCanvas(c_mps_name, c_mps_name);
            c_mps_v.push_back(c_mps);
            c_mps = nullptr;

            int n_param_1 = 100;
            int n_param_2 = 100;

            double param_1 = AdjustActs[param_1_ix];
            double sigma_1 = AdjustActs_Err[param_1_ix];
            double width_1 = 5.0;
            double param_1_min = param_1 + width_1 * sigma_1 * (-0.5); //(-n_param_1 / 2);
            double param_1_max = param_1 + width_1 * sigma_1 * 0.5; //(n_param_1 - n_param_1 / 2);

            double param_2 = AdjustActs[param_2_ix];
            double sigma_2 = AdjustActs_Err[param_2_ix];
            double width_2 = 5.0;
            double param_2_min = param_2 + width_2 * sigma_2 * (-0.5); //(-n_param_2 / 2);
            double param_2_max = param_2 + width_2 * sigma_2 * 0.5; //(n_param_2 - n_param_2 / 2);
            
            TString h_mps_name_base = "h_mps";
            TString h_mps_name = h_mps_name_base + "_" + param_1_ix_str_external + "_" + param_2_ix_str_external;

            std::cout << h_mps_name << " param_1=" << param_1 << " sigma_1=" << sigma_1
                                    << " param_1_min=" << param_1_min << " param_1_max=" << param_1_max
                                    << " param_2=" << param_2 << " sigma_2=" << sigma_2
                                    << " param_2_min=" << param_2_min << " param_2_max=" << param_2_max
                                    << std::endl;

            TH2D *h_mps = new TH2D(h_mps_name, h_mps_name,
                                   n_param_1, param_1_min, param_1_max,
                                   n_param_2, param_2_min, param_2_max); 
            h_mps_v.push_back(h_mps);
            //h_mps = nullptr;

            //h_mps->GetZaxis()->SetRangeUser(0.0, 1.0e+04);
            h_mps->SetContour(1000);
            
            TString param_1_name_str = TString(paramNameMap[param_1_ix_external].c_str());
            TString param_2_name_str = TString(paramNameMap[param_2_ix_external]);

            h_mps->GetXaxis()->SetTitle(param_1_name_str);
            h_mps->GetYaxis()->SetTitle(param_2_name_str);

            // reset params array
            for(int jx = 0; jx < n_params; ++ jx)
            {
                minuit->GetParameter(jx, params[jx], param_errs[jx]);
            }

            // get minimum
            double fval_min = 0.0;
            logLikelihood(n_params, nullptr, fval_min, params, 0);

            // modify parameters
            //for(int n_1 = 0; n_1 <= n_param_1; ++ n_1)
            for(int n_1 = 0; n_1 < n_param_1; ++ n_1)
            {
                //for(int n_2 = 0; n_2 <= n_param_2; ++ n_2)
                for(int n_2 = 0; n_2 < n_param_2; ++ n_2)
                {
                    // TODO: try using GetBinCenter() and looping over bins
                    // in combination with Fill method

                    double fval = 0.;

                    double a_1 = (double)n_1 / (double)n_param_1 - 0.5;
                    double a_2 = (double)n_2 / (double)n_param_2 - 0.5;

                    double t_param_1 = param_1 + width_1 * sigma_1 * a_1;
                    double t_param_2 = param_2 + width_2 * sigma_2 * a_2;

                    params[param_1_ix] = t_param_1;
                    params[param_2_ix] = t_param_2;

                    logLikelihood(n_params, nullptr, fval, params, 0);

                    /*
                    if(m == 50)
                    {
                        std::cout << "n=" << n << " a_nd150=" << a_nd150 << " p_nd150=" << p_nd150 << " fval=" << fval << std::endl;
                    }
                    */

                    //std::cout << n << " " << m << " " << fval << std::endl;
                    //std::cin.get();

                    //h_mps->Fill(n, m, fval);
                    //h_mps->SetBinContent(n_1 + 1, n_2 + 1, fval - fval_min);
                    double step_1 = width_1 * sigma_1 * (double)1 / (double)n_param_1;
                    double step_2 = width_2 * sigma_2 * (double)1 / (double)n_param_2;
                    h_mps->Fill(t_param_1 + step_1 / 2.0, t_param_2 + step_2 / 2.0, fval - fval_min);
                    // TODO: fval_min does not appear to always be the minimum

                    if(fval - fval_min <= 0.0)
                    {
                        std::cout << "dbg1: " << n_1 << " " << n_2 << " " << h_mps->GetBinContent(n_1, n_2) << std::endl;
                    }
                    if(n_1 == n_param_1 / 2)
                    {
                        if(n_2 == n_param_2 / 2)
                        {
                            std::cout << "dbg2: " << n_param_1 / 2 << " " << n_param_2 / 2 << " " << h_mps->GetBinContent(n_1, n_2) << std::endl;
                        }
                    }
                }
            }

            h_mps->Draw("colz");

            h_mps = nullptr;
            
        }
    }



    // reset params array
    for(int jx = 0; jx < n_params; ++ jx)
    {
        minuit->GetParameter(jx, params[jx], param_errs[jx]);
    }

    TCanvas *c_mps = new TCanvas("c_mps", "c_mps");

    int nd150_rot_2n2b_m4_index = paramNumberToMinuitParamNumberMap.at(0);

    int n_nd150 = 100;
    int n_mo100 = 100;

    double p_nd150_cv = AdjustActs[nd150_rot_2n2b_m4_index];
    double p_nd150_sigma = AdjustActs_Err[nd150_rot_2n2b_m4_index];
    double w_nd150 = 5.0;
    double p_min_nd150 = p_nd150_cv + w_nd150 * p_nd150_sigma * (-n_nd150 / 2);
    double p_max_nd150 = p_nd150_cv + w_nd150 * p_nd150_sigma * (n_nd150 - n_nd150 / 2);
    
    double p_mo100_cv = AdjustActs[mo100_99_rot_2n2b_m14_index];
    double p_mo100_sigma = AdjustActs_Err[mo100_99_rot_2n2b_m14_index];
    double w_mo100 = 5.0;
    double p_min_mo100 = p_mo100_cv + w_mo100 * p_mo100_sigma * (double)(-n_mo100 / 2);
    double p_max_mo100 = p_mo100_cv + w_mo100 * p_mo100_sigma * (double)(n_mo100 - n_mo100 / 2);
    
    std::cout << p_nd150_cv << ", " << p_nd150_sigma << ", " << p_min_nd150 << ", " << p_max_nd150 << std::endl;
    std::cout << p_mo100_cv << ", " << p_mo100_sigma << ", " << p_min_mo100 << ", " << p_max_mo100 << std::endl;

    TH2F *h_mps = new TH2F("h_mps", "h_mps", n_nd150, p_min_nd150, p_max_nd150, n_mo100, p_min_mo100, p_max_mo100); 
    //h_mps->GetZaxis()->SetRangeUser(0.0, 1.0e+04);
    h_mps->SetContour(256);
    h_mps->GetXaxis()->SetTitle("^{150}Nd");
    h_mps->GetYaxis()->SetTitle("^{100}Mo");

    // get minimum
    double fval_min = 0.0;
    logLikelihood(n_params, nullptr, fval_min, params, 0);

    for(int n = 0; n < n_nd150; ++ n)
    {
        for(int m = 0; m < n_mo100; ++ m)
        {
            double fval = 0.;

            double a_nd150 = (double)n / (double)n_nd150 - 0.5;
            double a_mo100 = (double)m / (double)n_mo100 - 0.5;

            double p_nd150 = p_nd150_cv + w_nd150 * p_nd150_sigma * a_nd150;
            double p_mo100 = p_mo100_cv + w_mo100 * p_mo100_sigma * a_mo100;

            params[nd150_rot_2n2b_m4_index] = p_nd150;
            params[mo100_99_rot_2n2b_m14_index] = p_mo100;

            logLikelihood(n_params, nullptr, fval, params, 0);

            /*
            if(m == 50)
            {
                std::cout << "n=" << n << " a_nd150=" << a_nd150 << " p_nd150=" << p_nd150 << " fval=" << fval << std::endl;
            }
            */

            //std::cout << n << " " << m << " " << fval << std::endl;
            //std::cin.get();

            //h_mps->Fill(n, m, fval);
            h_mps->SetBinContent(n + 1, m + 1, fval - fval_min);
        }
    }

    h_mps->Draw("colz");


    
    delete [] params;
    delete [] param_errs;
#endif
    ///////////////////////////////////////////////////////////////////////////



    // TODO: this no longer works, or does it?
    // needs to take into account the number of ENABLED free params
    // NOTE: 2020-04-16 fixed
    //Double_t CovMatrix[free_params.size()][free_params.size()];
    //minuit->mnemat(&CovMatrix[0][0],free_params.size());
    //Double_t CovMatrix[minuit->GetNumFreePars()][minuit->GetNumFreePars()];
    //Double_t *CovMatrix = new Double_t[minuit->GetNumFreePars() * minuit->GetNumFreePars()];
    number_free_params = minuit->GetNumFreePars();
    CovMatrix = new Double_t[number_free_params * number_free_params];
    for(int ix{0}; ix < number_free_params * number_free_params; ++ ix)
    {
        CovMatrix[ix] = 0.;
    }
    //minuit->mnemat(CovMatrix, minuit->GetNumFreePars());
    minuit->mnemat(CovMatrix, number_free_params);

  
    return minuit;

}



///////////////////////////////////////////////////////////////////////////////
// logLikelihood
///////////////////////////////////////////////////////////////////////////////


// TODO don't appear to work with parameters with more than one MC
void logLikelihood(Int_t & nPar, Double_t* /*grad*/, Double_t &fval, Double_t *p, Int_t /*iflag */)
{
    std::cout << "logLikelihood" << std::endl;
    std::cout << "p[0]=" << p[0] << " p[1]=" << p[1] << std::endl;


    // TODO: rebuild nd150 xi_31 paramter histogram here

    // there are i samples for each channel
    for(int i = 0; i < allDataSamples1D->GetEntries(); ++ i)
    {
        // note: I have really no idea how this is supposed to work
        // with multiple data channels, i = 0 here, for one data channel

        int channel = i;

        TH1F *h_before_reweight = nullptr;

        std::cout << "there are " << allMCSamples1D[channel]->GetEntries() << " objects" << std::endl;
        // new code to reweight 150Nd by xi_{31} parameter
        for(int i = 0; i < allMCSamples1D[channel]->GetEntries(); ++ i)
        {
            TH1F *tmpHist = (TH1F*)allMCSamples1D[channel]->At(i);
            TString tmpHist_name = tmpHist->GetName();
            // TODO: had to add "_fit" here - might not work after 1 iteration
            //if(tmpHist_name.CompareTo("hTotalE_nd150_rot_2n2b_m4_fit") == 0 ||
            //   tmpHist_name.CompareTo("hTotalE_nd150_rot_2b2n_m4_fit") == 0)
            if(tmpHist_name.Contains("nd150_rot_2n2b_m4") ||
               tmpHist_name.Contains("nd150_rot_2b2n_m4"))
            {
                //std::cout << "found the 150Nd MC" << std::endl;
                //std::cin.get();

                // TODO: this is very slow, gets re-built for each bin_ix


                //TH1F *tmpHist_draw1D_clone = nullptr;
                TH1F *tmpHist_reweight = nullptr;
                //reweight_apply(tmpHist_draw1D_clone, tmpHist_draw1D, ... );
                const double xi_31{p[1]};
                // TODO: this will not work if parameter number changes
                std::cout << "xi_31=" << xi_31 << std::endl;
                //const double xi_31{p[nPar-1]};
                //std::cout << "check that nPar=" << nPar << "=29/28 ?" << std::endl;
                // note: it isn't
                const double xi_31_baseline{0.368}; // TODO: change to actual value and pass in as argument somehow
                // fixed in parameter list file?

                // some debug stuff
                //TCanvas *ctmp = new TCanvas("ctmp", "ctmp");
                //h_nEqNull->Draw();
                //std::cin.get();

                //reweight_apply(tmpHist_reweight, "/mnt/ramdisknd150/Nd150_2eNg_output_truth_postprocessed.root", xi_31, xi_31_baseline, h_nEqNull, h_nEqTwo, psiN0, psiN2, bb_Q);
                // line below disabled
                //reweight_apply(tmpHist_reweight, "Nd150_2eNg_output_truth_postprocessed.root", xi_31, xi_31_baseline, h_nEqNull, h_nEqTwo, psiN0, psiN2, bb_Q);
                // TODO: after reweight function called, replace 150nd MC
                // in containers, or add a _reweight version to containers

                //std::cout << "name before: " << tmpHist->GetName() << std::endl;
                //std::cout << "name after: " << tmpHist_reweight->GetName() << std::endl;

                //std::cin.get();
                //std::cout << "removing i=" << i << std::endl;
                //allMCSamples1D[channel]->RemoveAt(i);
                //std::cout << "now there are " << allMCSamples1D[channel]->GetEntries() << " objects" << std::endl;
                //allMCSamples1D[channel]->Add(tmpHist_reweight);
                //std::cout << "and now there are " << allMCSamples1D[channel]->GetEntries() << " objects" << std::endl;
                //allMCSamples1D[i] = tmpHist_reweight;

                h_before_reweight = (TH1F*)tmpHist->Clone("h_before");

                std::cout << "tmpHist->GetName() -> " << tmpHist->GetName() << std::endl;
                reweight_apply(tmpHist_reweight, "Nd150_2eNg_output_truth_postprocessed.root", xi_31, xi_31_baseline, h_nEqNull, h_nEqTwo, psiN0, psiN2, bb_Q);
                allMCSamples1D[channel]->RemoveAt(i);
                allMCSamples1D[channel]->Add(tmpHist_reweight);

                /*
                for(int i{1}; i < h_before_reweight->GetNbinsX(); ++ i)
                {
                    Double_t bin1 = h_before_reweight->GetBinContent(i);
                    Double_t bin2 = tmpHist_reweight->GetBinContent(i);
                    if(std::abs(bin1 - bin2) > 1.0e-2)
                    {
                        std::cout << "bin1=" << bin1 << " bin2=" << bin2 << std::endl;
                    }
                }
                */

            }
            else
            {
                //std::cout << "it is not " << tmpHist->GetName() << std::endl;
                //std::cin.get();
            }

        }
   }




    // TODO: add check here to see if any disabled parameters are accessed


    //std::cout << "ll" << std::endl;

    double loglik = 0.; 
    //double tmp;


//   std::cout << "getting 1D histograms" << std::endl;

    TH1F *tmpData1D;
    // std::cout << allDataSamples1D->GetEntries()  << std::endl;

    // there are i samples for each channel
    for(int i = 0; i < allDataSamples1D->GetEntries(); ++ i)
    {

        // allDataSamples1D only contains one object
        
        TString i_str;
        i_str.Form("%i", i);
        //std::cout << i << std::endl;

        // TODO: can I remove this Clone() call safely to improve speed?
        //tmpData1D = (TH1F*)allDataSamples1D->At(i)->Clone("tmpData1D" + i_str + "_");
        tmpData1D = (TH1F*)allDataSamples1D->At(i);

        // std::cout << tmpData1D->Integral() << std::endl;

        int nBinsX = tmpData1D->GetNbinsX();
        for(int bin_ix = 1; bin_ix <= nBinsX; ++ bin_ix)
        {
            Int_t nData = (Int_t)tmpData1D->GetBinContent(bin_ix);
            // i is the index of the sample ?
            // ix is the bin index
            // p is a pointer to an array of parameter values

            // TODO:
            // think there is a bug here
            // i appears to be the index of the data sample not the MC sample?
            // it becomes channel index
            double nMC = getNumberMC1D(i, bin_ix, p);

            //std::cout << "for bin_ix=" << bin_ix << " nMC=" << nMC << " nData=" << nData << std::endl;

            //Int_t new_i = -1;
            //TString name = names.at(i);
            //new_i = paramNameToNumberMap[name];
            //std::cout << "the new i value is new_i=" << new_i << std::endl;
            //double nMC = getNumberMC1D(new_i, ix, p);

            //std::cin.get();

            
            // 2020-04-17 Notes:
            // output histograms do not look right (2d MPS plots)
            //
            // range of values for parameter include values such as
            // -8 to 10
            // -15 to 15
            // these ranges seem too large? / are indicating very large uncertainty
            // update: caused by a bug in the min/max parameter settings, now fixed
            // (not fixed in h_mps)
            //
            // there is a white square in the center. what value does this have?
            // is it negative, or zero?
            //
            // h_mps looks different to h_mps_10_0, they should be the same!
            // update: they are the same if ranges set the same and fval_min
            // is subtracted from fval before filling
            //
            // my guess was that the change in parameter values and thus n_mc
            // is having a much weaker effect compared to the penalty term
            // why is this?
            // this may not be correct, since all plots appear identical
            // indicating that something is not being computed correctly

            if(nMC > 0.)
            {
                Double_t poisson = TMath::Poisson(nData, nMC);
                if(poisson > 0.)
                {
                    //std::cout << "adding loglik value : " << TMath::Log(poisson) << " bin_ix=" << bin_ix << " poisson=" << poisson << std::endl;
	                loglik += TMath::Log(poisson);
                    // adding positive number makes fval go down
                    // NOTE: Log(poisson) is always negative! so fval goes UP NOT DOWN
                    // log is taken here, should I take log of the penalty term?
                    // TODO: answer above question
                    // TODO: are there any conditions for which this can be negative?
                    // poisson is a probability, so values are between 0 and 1 which means that
                    // log of this value is always negative
                }
                else
                {
                    // MARK have not tested this yet
                    // can this ever happen? is this the correct way to deal
                    // with the problem?
                    //std::cout << "catch: poisson" << std::endl;
                    // this does not appear to happen

                    //std::cout << "adding penalty of -10. to loglik bin_ix=" << bin_ix << " nMC=" << nMC << " poisson=" << poisson << std::endl;
                    // TODO: should this be removed? check for bins where ndata = 0?

                    // TODO: there were a lot of failures here
                    //std::cout << "ERROR: failed to evaluate TMath::Poisson()=" << poisson << " -> nData=" << nData << ", nMC=" << nMC << "; bin_ix=" << bin_ix << std::endl;
	                //loglik -= 10.;
                    // subtracting positive number makes fval go up
                    // TODO: this may not be a large enough penalty

                    // TODO: removed this
                }
            }
            else
            {
                // MARK have not tested this yet
                // not sure we are dealing with zero bins correctly, should
                // ignore?
                //std::cout << "catch2: poisson" << std::endl;
                // this appears to happen a lot

                // if nMC <= 0., just add penalty and cout nothing
                // subtracting positive number makes fval go up
                //loglik -= 10.;
                // 2020-04-17: removed, should I have something here?
            }
        
        } //~bins
    } //channels
    
    
    //  std::cout << "getting 2D histograms" << std::endl;


    /*
    TH2F *tmpData2D;
    // std::cout << allDataSamples->GetEntries()  << std::endl;
    for(int i = 0; i < allDataSamples2D->GetEntries(); ++ i)
    {
        // there are i samples for each channel
        TString i_str;
        i_str.Form("%i", i);

        tmpData2D = (TH2F*)allDataSamples2D->At(i)->Clone("tmpData2D" + i_str + "_");

        int nBinsX = tmpData2D->GetNbinsX();
        int nBinsY = tmpData2D->GetNbinsY();

        for(int ix = 1; ix <= nBinsX; ++ ix)
        {
            for(int iy = 1; iy <= nBinsY; ++ iy)
            {
                Int_t nData = (Int_t)tmpData2D->GetBinContent(ix, iy);
	            double nMC = getNumberMC2D(i, ix, iy, p);

	            if(nMC > 0 && TMath::Poisson(nData, nMC) > 0)
                {
	                loglik += TMath::Log(TMath::Poisson(nData, nMC));
	            }
                else
                {
	                loglik -= 10.;
	            }
            } //~ybins
        } //~xbins
    } //channels
    */
 
    // add constraints to improve likelihood
    //will eventually add gaussian constraint
    // disabled
    /*
    for ( int i = 0; i < 14; i++ ) {
        double constraint = 0.;
        //constraint = pow( (p[i] - 1) / (0.01), 2);
        // std::cout << i << "   " << p[i];
        bool fixed = false;
        for ( int j = 0; j < fixed_params.size(); j++ ) {
            if ( i == fixed_params.at(j) ) fixed = true;
        }

        //if ( (i <2) || (i == 7 ) || (i==10) || (i == 14) || (i==15) || (i==16) || (i == 18) || (i == 19) || (i == 21) || 9i==22) ) {
        // constraint = pow( (p[i] - 1.) / (paramActErrMap[i].front() / paramActMap[i].front()), 2);
        constraint = TMath::Gaus(p[i],1.,(paramActErrMap[i].front() / paramActMap[i].front()),true);
      

        // std::cout << "   " << constraint << std::endl;
        //  loglik += log(constraint);  
	    // }
    }
    */

    int mode = MODE_PARAM_UNDEFINED;

    // penalty terms section
    double penalty_term = 0.0;

    // TODO: I don't like this - should loop over the enabled params
    // however, this should still work as it is
    for(int i = 0; i < numberParams; ++ i)
    {

        if(std::find(enabled_params.begin(), enabled_params.end(), i) == enabled_params.end())
        {
            //std::cout << "parameter number " << param_number << " is disabled" << std::endl;
            //std::cout << "ERROR: i=" << i << " - parameter is DISABLED" << std::endl;
            //std::cin.get();
            continue;
        }


        if(thePhase == 0)
        {
            mode = paramConstrainModeP1Map[i];

            //if(paramConstrainModeP1Map[i] == MODE_PARAM_SOFT)
            if(mode == MODE_PARAM_SOFT)
            {
                // do nothing, soft constraint will be applied below
            }
            //else if(paramConstrainModeP1Map[i] == MODE_PARAM_HARD)
            else if(mode == MODE_PARAM_HARD)
            {
                // parameter fixed by minuit, continue to next param
                //continue;
                // although param is fixed by minuit, we still want to add
                // the penalty term, if available (TODO: will it always
                // be available?)
                // can check if error == 0.0 first
                
                // NOTE: changed to ignore HARD
                continue;
            }
            //else if(paramConstrainModeP1Map[i] == MODE_PARAM_FREE)
            else if(mode == MODE_PARAM_FREE)
            {
                // no constraint to apply, continue to next param
                continue;
            }
            else
            {
                std::cout << "ERROR: Invalid value in paramConstrainModeP1Map: paramConstrainModeP1Map[" << i << "]=" << paramConstrainModeP1Map[i] << std::endl;
            }
        }
        else if(thePhase == 1)
        {
            mode = paramConstrainModeP2Map[i];

            //if(paramConstrainModeP2Map[i] == MODE_PARAM_SOFT)
            if(mode == MODE_PARAM_SOFT)
            {
                // do nothing, soft constraint will be applied below
            }
            //else if(paramConstrainModeP2Map[i] == MODE_PARAM_HARD)
            else if(mode == MODE_PARAM_HARD)
            {
                // parameter fixed by minuit, continue to next param
                //continue;
                // although param is fixed by minuit, we still want to add
                // the penalty term, if available (TODO: will it always
                // be available?)
                // can check if error == 0.0 first
                
                // NOTE: changed to ignore HARD
                continue;
            }
            //else if(paramConstrainModeP2Map[i] == MODE_PARAM_FREE)
            else if(mode == MODE_PARAM_FREE)
            {
                // no constraint to apply, continue to next param
                continue;
            }
            else
            {
                std::cout << "ERROR: Invalid value in paramConstrainModeP2Map: paramConstrainModeP2Map[" << i << "]=" << paramConstrainModeP2Map[i] << std::endl;
            }
        }
        else
        {
            std::cout << "ERROR: Invalid value for thePhase: thePhase=" << thePhase << std::endl;
        }
        
        // soft constraint is applied here
        double constraint = 0.;
        double error = 0.;
        // NOTE: these values read from parameter list file and thus are in
        // units of activity (Bq)

        /*
        if(thePhase == 0)
        {
            for(int j = 0; j < paramInitValueP1Map[i].size(); ++ j)
            {
                param_init_value += paramInitValueP1Map[i].at(j);
            }
            for(int j = 0; j < paramInitErrorP1Map[i].size(); ++ j)
            {
                param_init_error += paramInitErrorP1Map[i].at(j);
            }
        }
        else if(thePhase == 1)
        {
            for(int j = 0; j < paramInitValueP2Map[i].size(); ++ j)
            {
                param_init_value += paramInitValueP2Map[i].at(j);
            }
            for(int j = 0; j < paramInitErrorP2Map[i].size(); ++ j)
            {
                param_init_error += paramInitErrorP2Map[i].at(j);
            }
        }
        */
        if(thePhase == 0)
        {
            //for(int j = 0; j < paramConstraintValueP1Map[i].size(); ++ j)
            //{
                //constraint += paramConstraintValueP1Map[i].at(j);
            //}
            //for(int j = 0; j < paramConstraintErrorP1Map[i].size(); ++ j)
            //{
                //error += paramConstraintErrorP1Map[i].at(j);
            //}
            constraint = paramConstraintValueP1Map[i];
            error = paramConstraintErrorP1Map[i];
            // TODO: this is a bit weird and perhaps not what would
            // be expected looking at the input file
            // perhaps it would be better to use .at(0)
            // and ignore all later values in parameters
            // in particular, errors do not usually add, but usually
            // add in quadrature, which is not what happens here
            // TODO: consider adding errors in quadrature, or
            // changing input file to a list of samples as a single parameter
            // rather than having multiple parameter numbers, each with a
            // single sample
            //
            //constraint = paramConstraintValueP1Map[i];
            //error = paramConstraintErrorP1Map[i];
        }
        else if(thePhase == 1)
        {
            //for(int j = 0; j < paramConstraintValueP2Map[i].size(); ++ j)
            //{
            //    constraint += paramConstraintValueP2Map[i].at(j);
            //}
            //for(int j = 0; j < paramConstraintErrorP2Map[i].size(); ++ j)
            //{
            //    error += paramConstraintErrorP2Map[i].at(j);
            //}
            constraint = paramConstraintValueP2Map[i];
            error = paramConstraintErrorP2Map[i];

            //constraint = paramConstraintValueP2Map[i];
            //error = paramConstraintErrorP2Map[i];
        }
        else
        {
            std::cout << "ERROR: Invalid value for thePhase: thePhase=" << thePhase << std::endl;
        }
        // TODO: can optimize this code

        if(error < 0.0)
        {
            std::cout << "ERROR: Invalid error value: error=" << error << std::endl;
        }

        // check if hard parameter and error == 0.0
        if(mode == MODE_PARAM_HARD)
        {
            if(error == 0.0)
            {
                // this parameter is a "constant" (according to minuit)
                // so ignore
                continue;
            }
            else
            {
                // do nothing, add constraint for hard param
            }
        }

        //double param_value = p[i];
        // i is index of the parameter number (external / non minuit)
        // convert to internal parameter number (minuit param number)
        int j = paramNumberToMinuitParamNumberMap.at(i);
        if(j < numberEnabledParams)
        {
            // ok
        }
        else
        {
            throw std::runtime_error("error: invalid value of j (internal param number)");
        }
        double param_value = p[j];
        // this parameter is from minuit internal and thus is in minuit
        // internal units (not Bq)
        // have to convert to Bq units
    
        // convert to Bq
        // multiply by the initial value
        double activity_value_Bq = 0.0;
        if(thePhase == 0)
        {
            activity_value_Bq = paramInitValueP1Map[i];
        }
        else if(thePhase == 1)
        {
            activity_value_Bq = paramInitValueP2Map[i];
        }
        else
        {
            std::cout << "ERROR: Invalid value for thePhase: thePhase=" << thePhase << std::endl;
        }
        //param_value *= activity_value_Bq;

        double value = param_value * activity_value_Bq;
        //double penalty = std::pow((param_value - constraint) / error, 2.0);
        double penalty = std::pow((value - constraint) / error, 2.0);
        // TODO: penalty term should be a Gaussian constraint?
        // NOTE: gaussian constraint, after taking log, is the same as
        // quadratic constraint - however there is the issue of a remaining
        // constant term which I do not yet fully understand
        // TODO: is this the correct error term?
        // error on constraint rather than error on current fit value?
        // TODO: is the value correct?
        // NOTE: I think it's all correct

        //double value = param_value * activity_value_Bq;
        //double penalty = std::pow((value - constraint) / error, 2.0);

        //std::cout << "adding (but it has to be subtracting!) penalty (i=" << i << ") : " << penalty << std::endl;

        // subtracting positive number makes fval go up
        //loglik -= penalty;
        penalty_term += penalty;
    }
  
    //fval = -2.0 * loglik; 
    // equivalent to
    //fval = -2.0 * (loglik_no_penalty_terms + penalty_terms); 
    //fval = -2.0 * loglik_no_penalty_terms - 2.0 * penalty_terms; 
    //fval = -2.0 * loglik_no_penalty_terms + 2.0 * penalty_terms_positive_sign;
    // then fix factor of 2.0 bug to get
    //fval = -2.0 * loglik_no_penalty_terms + penalty_terms_positive_sign;
    fval = -2.0 * loglik + penalty_term;
    //tmpData->Delete();

    // hook
    global_chisquare = fval;

    return;

}


Double_t getNumberMC1D(Int_t channel, Int_t bin_ix, Double_t *p) {


    //std::cout << "getNumberMC1D() called with channel=" << channel << " binx=" << binx << std::endl;
    //std::cout << "printing contents of parameter array" << std::endl;
    //for(int i = 0; i < numberParams; ++ i)
    //{
    //    std::cout << "i=" << i << " p[i]=" << p[i] << std::endl;
    //}


    // std::cout <<"getting number of1D MC... "  <<channel << std::endl;

    double nMC = 0.;

    // (1) grab a hist from the sample list of this channel
    // (2) figure out which parameter it corresponds to
    TH1F *tmpHist;
    //int which_param;

    //std::cout << "getting number of MC... "  << channel << std::endl;

    //std::cout << allMCSamples[channel]->GetEntries() << std::endl;


    // copied from above
    //for(int k = 0; k < allMCSamples1D[channel]->GetEntries(); k++)
    for(int j = 0; j < allMCSamples1D[channel]->GetEntries(); j++)
    {

        //std::cout << "DEBUG: allMCSamples1D[channel]->GetEntries() -> " << allMCSamples1D[channel]->GetEntries() << std::endl;
        //std::cin.get();
        // disabled params do not appear in this array

        //tmpHist = (TH1F*)allMCSamples1D[channel]->At(k);
        //tmpHist = (TH1F*)allMCSamples1D[channel]->At(j)->Clone();
        tmpHist = (TH1F*)allMCSamples1D[channel]->At(j);
        //tmpHist_draw1D = (TH1F*)allMCSamples1D[channel]->At(j)->Clone();

        //if ( tmpHist->Integral() == 0 ) continue;
       
        TString tmpName = tmpHist->GetName();
        //TString tmpName;// = tmpHist->GetName();
        //try
        //{
        //tmpName = tmpHist->GetName();
        //}
        //catch(std::exception &e)
        //{
        //    std::cout << "j=" << j << std::endl;
        //    std::cout << e.what() << std::endl;
        //    throw e;
        //}
        //TString tmpName = tmpHist_draw1D->GetName();
        
        #if PRINT_MO
        if(tmpName.Contains("mo100_99_rot_2n2b_m14") && bin_ix == 1)
        {
            std::cout << "caught mo100" << std::endl;
            //std::cin.get();
        }
        #endif
        

        //std::cout << "looking for " << tmpName << std::endl;
        //bool foundParam = false;
        int which_param = -1;
        bool found_param = false;

        /*
        for(int i = 0; (i < numberParams) && !foundParam; i++)
        {
            // std::vector<int>::iterator it = find(fixed_params.begin(), fixed_params.end(), i);
            //if ( it != fixed_params.end() ) continue;      

            for(int j = 0; (j < paramNameMap[i].size()) && !foundParam; j++)
            {
                //std::cout << "is it: " << paramNameMap[i].at(j) << std::endl;

	            if(tmpName.Contains(paramNameMap[i].at(j)))
                {
	                foundParam = true;
                    which_param = i;
                    //std::cout << "match found: " << tmpName << " -> " << paramNameMap[i].at(j) << ", which_param=" << which_param << std::endl;
	            }

            } //~j searching through array of params for the right one
        } //~i list of isotopes with same parameter
        */

        //std::cout << "NEW CODE" << std::endl;
        //try
        //{
        // TODO: remove TString
        {
            std::string tmp_hist_name(tmpName);
            auto i_start = tmp_hist_name.find('_') + 1;
            auto i_end = tmp_hist_name.rfind('_');
            if(i_end - i_start > 0)
            {
                std::string tmp_sample_name = tmp_hist_name.substr(i_start, i_end - i_start);
                //std::cout << "tmp_sample_name=" << tmp_sample_name << std::endl;
                if(MCNameToParamNumberMap.count(tmp_sample_name) > 0)
                {
                    int paramNumber = MCNameToParamNumberMap.at(tmp_sample_name);
                    // TODO: removed std::string, change tmpName type to be std::string from TString
                
                    //std::cout << "paramNumber=" << paramNumber << " -> tmp_sample_name=" << tmp_sample_name << " ~> tmpName=" << tmpName << std::endl;                    
                    //which_param = paramNumber;
                    which_param = paramNumberToMinuitParamNumberMap.at(paramNumber);
                    found_param = true;

                    //std::cout << "DEBUG: found parameter with minuit (internal) number: " << which_param << std::endl;
                    //std::cin.get();
                }
                else
                {
                   std::cout << "ERROR: could not find " << tmp_sample_name << " in MCNameToParamNumberMap" << std::endl;
                }
            }
        }
        /*
        }
        catch(std::exception &e)
        {
            std::cout << "e.what(): " << e.what() << std::endl;
            std::cout << "tmpName=" << tmpName << std::endl;
            std::cout << "contents of map" << std::endl;
            for(auto it = MCNameToParamNumberMap.cbegin(); it != MCNameToParamNumberMap.cend(); ++ it)
            {
                std::cout << it->first << " -> " << it->second << std::endl;
            }
        }
        std::cin.get();
        */

        //std::cout << "bin_ix=" << bin_ix << " tmpHist->GetName()=" << tmpHist->GetName() << " which_param=" << which_param << std::endl;

        /*
        if(foundParam)
        {
            //std::cout << "adding to nMC with index of which_param=" << which_param << std::endl;
            // TODO: think this is collecting the wrong parameter? or is it?
            nMC += p[which_param] * tmpHist->GetBinContent(binx);
        }
        else
        {
            std::cout << "error could not find histogram: " << tmpName << std::endl;
        }
        */

        if(found_param == true)
        {
        #if PRINT_MO
            if(tmpName.Contains("mo100_99_rot_2n2b_m14") && bin_ix == 1)
            {
                //std::cout << "caught mo100" << std::endl;
                //std::cin.get();
                std::cout << "which_param=" << which_param << std::endl;
            }
        #endif

            // check here to see if param is disabled
            // TODO
            //
            // since which_param must be an enabled parameter, no longer need this
            // leave for now as a check?
            //if(std::find(enabled_params.begin(), enabled_params.end(), which_param) == enabled_params.end())
            if(std::find(enabled_params.begin(), enabled_params.end(), minuitParamNumberToParamNumberMap.at(which_param)) == enabled_params.end())
            {
                //std::cout << "parameter number " << param_number << " is disabled" << std::endl;
                std::cout << "ERROR: which_param=" << which_param << " - parameter is DISABLED" << std::endl;
                std::cin.get();
            }

            //std::cout << "found histogram: tmpName=" << tmpName << " which_param=" << which_param << std::endl;

            /*
            // scale histogram to correct size using output parameter
            // from fit
            //tmpHist_draw1D->Scale(AdjustActs[which_param]);
            tmpHist->Scale(AdjustActs[which_param]);
            
            //if(tmpHist_draw1D->Integral() > 0)
            if(tmpHist->Integral() > 0)
            {
                //stacks1D[channel]->Add(tmpHist_draw1D);
                stacks1D[channel]->Add(tmpHist);
            }
            else
            {
                std::cout << "not adding to stack, Integral() <= 0" << std::endl;
            }
            */

        #if PRINT_MO
            if(tmpName.Contains("mo100_99_rot_2n2b_m14"))
            {
                //std::cout << "caught mo100" << std::endl;
                //std::cin.get();
                std::cout << "p[" << which_param << "]=" << p[which_param] << ", tmpHist->GetBinContent(" << bin_ix << ")=" << tmpHist->GetBinContent(bin_ix) << std::endl;
            }
        #endif
            nMC += p[which_param] * (double)tmpHist->GetBinContent(bin_ix);
            



            //std::cout << "DEBUG: tmpName=" << tmpName << " which_param=" << which_param << std::endl;
            //std::cin.get();


            //std::cout << "contents of map" << std::endl;
            //for(auto it = MCNameToParamNumberMap.cbegin(); it != MCNameToParamNumberMap.cend(); ++ it)
            //{
            //    std::cout << it->first << " -> " << it->second << std::endl;
            //}

            //std::cout << "adding to nMC: which_param=" << which_param << ", bin_ix=" << bin_ix << ", " << p[which_param] * tmpHist->GetBinContent(bin_ix) << std::endl;
            //std::cin.get();
        }
        else
        {
            std::cout << "error could not find histogram: tmpName=" << tmpName << std::endl;
        } 

        #if PRINT_MO
        if(tmpName.Contains("mo100_99_rot_2n2b_m14") && bin_ix == 50)
        {
            //std::cout << "caught mo100" << std::endl;
            std::cin.get();
        }
        #endif

    }

    return nMC;


    /*
    for(int k = 0; k < allMCSamples1D[channel]->GetEntries(); k++)
    {
        tmpHist = (TH1F*)allMCSamples1D[channel]->At(k);

        //if ( tmpHist->Integral() == 0 ) continue;
       
        TString tmpName = tmpHist->GetName();

        //std::cout << "looking for " << tmpName << std::endl;
        bool foundParam = false;

        for(int i = 0; (i < numberParams) && !foundParam; i++)
        {
            // std::vector<int>::iterator it = find(fixed_params.begin(), fixed_params.end(), i);
            //if ( it != fixed_params.end() ) continue;      

            for(int j = 0; (j < paramNameMap[i].size()) && !foundParam; j++)
            {
                //std::cout << "is it: " << paramNameMap[i].at(j) << std::endl;

	            if(tmpName.Contains(paramNameMap[i].at(j)))
                {
	                foundParam = true;
                    which_param = i;
                    //std::cout << "match found: " << tmpName << " -> " << paramNameMap[i].at(j) << ", which_param=" << which_param << std::endl;
	            }

            } //~j searching through array of params for the right one
        } //~i list of isotopes with same parameter

        if(foundParam)
        {
            //std::cout << "adding to nMC with index of which_param=" << which_param << std::endl;
            // TODO: think this is collecting the wrong parameter? or is it?
            nMC += p[which_param] * tmpHist->GetBinContent(binx);
        }
        else
        {
            std::cout << "error could not find histogram: " << tmpName << std::endl;
        }	 
    } //each histogram

    return nMC;
    */
}


// TODO: has not been updated
Double_t getNumberMC2D(Int_t channel, Int_t binx, Int_t biny, Double_t *p) {

  // std::cout <<"getting number of 2D MC... "  <<channel << std::endl;

  double nMC = 0;

  // (1) grab a hist from the sample list of this channel
  // (2) figure out which parameter it corresponds to
  TH2F *tmpHist2;
  int which_param;
  //  std::cout <<"getting number of MC... "  <<channel << std::endl;

  //std::cout << allMCSamples[channel]->GetEntries() << std::endl;

  for ( int k = 0; k < allMCSamples2D[channel]->GetEntries(); k++ ) {
    tmpHist2 = (TH2F*)allMCSamples2D[channel]->At(k);
    TString tmpName = tmpHist2->GetName();

    //  std::cout << "looking for " << tmpName << std::endl;
    bool foundParam = false;

    for ( int i = 0; (i < numberParams) && !foundParam; i++ ) {

      std::vector<int>::iterator it = find(fixed_params.begin(), fixed_params.end(), i);
      if ( it != fixed_params.end() ) continue;    

      for ( int j = 0; (j < paramNameMap[i].size()) && !foundParam; j++ ) {
	//	   std::cout <<"is it...  " <<paramNameMap[i].at(j) << std::endl;

	if ( tmpName.Contains(paramNameMap[i].at(j)) ) {
	  foundParam = true;
          which_param = i;
	}

      }//~j searching through array of params for the right one
    }//~i list of isotopes with same parameter

    if ( foundParam ) {nMC += p[which_param]*tmpHist2->GetBinContent(binx,biny);}
    // else {std::cout << "error could not find histogram: " << tmpName << std::endl;}	 
  }//each histogram

  return nMC;
}
