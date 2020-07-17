#include "TH1.h"
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
#include <vector>
//#include <TMinuit.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnMigrad.h>


// 
#include <time.h>


// note: these must appear in correct order and after general includes above
#include "newLogLikFitter.h"
#include "newLogLikFitter_aux.h"
#include "newLogLikFitter_print.h"
#include "newLogLikFitter_read_parameternames_lst.h"
#include "newLogLikFitter_book1DHistograms.h"
#include "newLogLikFitter_book2DHistograms.h"
#include "newLogLikFitter_reweight.h"
#include "newLogLikFitter_stacker_helper.h"
#include "newLogLikFitter_drawchannel.h"
#include "newLogLikFitter_logpoisson.h"
#include "newLogLikFitter_buildfakedata.h"
#include "newLogLikFitter_draw.h"
#include "newLogLikFitter_draw_2D.h"
#include "newLogLikFitter_draw_outputdiff.h"
#include "newLogLikFitter_draw_all.h"
#include "MinimizeFCNAxialVector.h"
#include "newLogLikFitter_chisquaretest.h"
#include "newLogLikFitter_test.h"
#include "newLogLikFitter_fitBackgrounds.h"


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
//  TODO: rescale xi_31 parameter to restore homogenaity in AdjustActs


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



// NOTE: fake data
// to switch between real data and fake data mode change switches
// in draw function as well as log likelihood function


//-----------------------------------------
//    Functions
//----------------------------------------
void loadFiles();


//void fitBackgrounds(double *AdjustActs, double *AdjustActs_Err, double*& CovMatrix, int& number_free_params, Int_t thePhase);
TMinuit * fitBackgrounds(double *AdjustActs, double *AdjustActs_Err, double*& CovMatrix, int& number_free_params, Int_t thePhase);















void loadFiles()
{



    std::cout << std::scientific;

    // load data
    std::cout << "attempting to load spectral data from file" << std::endl;

    bb_Q = 3.368;
    // ramdisk, should be faster?
    std::size_t count_G0 = LD_REWEIGHT_DATA_2(h_nEqNull, "/home/ecb/150Nd/150Nd-data/dG150Nd/G0/dG0.dat", "h_nEqNull", "nEqNull", 0.0, bb_Q);
    std::size_t count_G2 = LD_REWEIGHT_DATA_2(h_nEqTwo,  "/home/ecb/150Nd/150Nd-data/dG150Nd/G2/dG2.dat", "h_nEqTwo", "nEqTwo",  0.0, bb_Q);

    // phase space integrals
    const Double_t G0_ps_integral_MeV = 0.420438E-45;
    const Double_t G0_ps_integral_yrinv = 0.201577E-16;
    const Double_t G2_ps_integral_MeV = 0.744684E-45;
    const Double_t G2_ps_integral_yrinv = 0.357034E-16;

    // TODO: what is the value for the basline for Nd150? is this for 100Mo
    //const Double_t xi_31_baseline{0.296}; // TODO: this is WRONG change it
//    const Double_t xi_31_baseline{0.0}; // TODO: this is WRONG change it
    xi_31_baseline = 0.0;
    // TODO: multiple instances, move to header file
    // TODO: also change in parameter list file
    //Double_t xi_31_init = 0.8;
    //Double_t xi_31_init = 0.0; //1.13; //0.296; // change to baseline value for testing purposes
    // TODO: remove xi_31_init variable and replace with paramInitValuePXMap[1]
    // and then change index 1 for a sensible way of looking up the index
//    Double_t xi_31_init = 2.63e-01; //1.13; //0.296; // change to baseline value for testing purposes
    last_xi_31_parameter_value = 0.0;


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

    psiN0 = G0_ps_integral_MeV;
    psiN2 = G2_ps_integral_MeV; // TODO: check this is the correct option
    
    std::cout << "histogram format constructed" << std::endl;


    if(0)
    {
        draw_inputdata();
    }

    std::cout << std::fixed << std::endl;

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPalette(kBird);
    gStyle->SetPalette(kBrownCyan);
    gStyle->SetPalette(kLightTemperature);
    //gStyle->SetNumberContours(1000);

    // TODO: C++ syntax, should be TObjArray?
    // without parens
    allDataSamples1D = new TObjArray();
    allDataSamples2D = new TObjArray();
    // TODO could call build_fake_data() in this function after load?
    allFakeDataSamples1D = nullptr;
    allFakeDataSamples2D = nullptr;



    

    // read parameter_name.lst file
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
    book1DHistograms(1, "2e_", "P" + Phase, "hSingleEnergy_");
    book1DHistograms(2, "2e_", "P" + Phase, "hHighEnergy_");
    book1DHistograms(3, "2e_", "P" + Phase, "hLowEnergy_");
    book1DHistograms(4, "2e_", "P" + Phase, "hEnergySum_");
    book1DHistograms(5, "2e_", "P" + Phase, "hEnergyDiff_");
    book2DHistograms(0, "2e_", "P" + Phase, "hHighLowEnergy_");
  
    // Array to hold activity adjustment parameters
    //Double_t AdjustActs[numberParams];
    //Double_t AdjustActs_Err[numberParams];
    //for(Int_t ix{0}; ix < numberParams; ++ ix)
    Double_t AdjustActs[numberEnabledParams];
    Double_t AdjustActs_Err[numberEnabledParams];
    for(Int_t ix{0}; ix < numberEnabledParams; ++ ix)
    {
        AdjustActs[ix] = 1.0;
        AdjustActs_Err[ix] = 0.5; // TODO: set using parameter_names.list
    }
    // TODO: fix this
    double xi_31_init_value = 0.0;
    double xi_31_init_error = 0.0;
    get_paramInitValueError(thePhase, 1, xi_31_init_value, xi_31_init_error);
    AdjustActs[1] = xi_31_init_value;

    std::cout << "initial value is set to " << xi_31_init_value << std::endl;

    std::vector<double> init_par;
    std::vector<double> init_err;
    for(Int_t ix{0}; ix < numberEnabledParams; ++ ix)
    {
        init_par.push_back(AdjustActs[ix]);
        init_err.push_back(AdjustActs_Err[ix]);
    }

    /*
    Int_t number_free_params = minuit->GetNumFreePars();
    Double_t *CovMatrix = new Double_t[number_free_params * number_free_params];
    for(int ix{0}; ix < number_free_params * number_free_params; ++ ix)
    {
        CovMatrix[ix] = 0.;
    }
    */
    

    if(0)
    {

        // TODO: somewhere I copy using numberParams instead of numberEnabledParams
        // if I remember correctly, need to find and fix this

        Double_t AdjustActs_copy[numberEnabledParams];
        Double_t AdjustActs_Err_copy[numberEnabledParams];

        for(int i = 0; i < numberEnabledParams; ++ i)
        {
            AdjustActs_copy[i] = AdjustActs[i];
            AdjustActs_Err_copy[i] = AdjustActs_Err[i];
        }
        

        //do_test_xi_31_test1(AdjustActs_copy, AdjustActs_Err_copy);
        return 0;
    }



    int number_free_params = -1;
    double *CovMatrix = nullptr;




    //fitBackgrounds(AdjustActs, AdjustActs_Err, CovMatrix, number_free_params, thePhase);
    // fit of backgrounds disabled.
    // need to fit Nd150 parameter only
    // disabled in fitBackgrounds function
    
    // needs to remain enabled to define parameters in minuit
//    TMinuit *minuit = nullptr;
//    TMinuit *minuit = fitBackgrounds(AdjustActs, AdjustActs_Err, CovMatrix, number_free_params, thePhase);
//    fitBackgrounds_exec(minuit, AdjustActs, AdjustActs_Err);

    // create minimizer
    ROOT::Minuit2::MnUserParameterState theParameterState;
    ROOT::Minuit2::VariableMetricMinimizer theMinimizer;
    MinimizeFCNAxialVector theFCN;

    ROOT::Minuit2::FunctionMinimum FCN_min =
        fitBackgrounds(
            theParameterState,
            theMinimizer,
            theFCN,
            AdjustActs,
            AdjustActs_Err,
            CovMatrix,
            number_free_params,
            thePhase);

//    theParameterState.add();

    // minimize
    //ROOT::Minuit2::FunctionMinimum FCN_min = theMinimizer.Minimize(theFCN, init_par, init_err);
    //ROOT::Minuit2::FunctionMinimum FCN_min = theMinimizer.Minimize(theFCN, init_par, init_err);
    std::cout << "Minimization finished" << std::endl;
    std::cout << "minimum: " << FCN_min << std::endl;
    std::cout << "chi2: " << FCN_min.Fval() << std::endl;
    std::cout << "edm: " << FCN_min.Edm() << std::endl;


    // draw my data for a parameter xi = 0.0
    double fval = 0.;
    int n_params = theParameterState.Params().size();
    std::vector<double> params = theParameterState.Params();
    std::vector<double> param_errs = theParameterState.Errors();
    
    std::cout << params[0] << " " << params[1] << std::endl;

    //logLikelihood(n_params, nullptr, fval, params, 0);
    fval = theFCN.operator()(params);

    std::cout << "fval=" << fval << " for params[0]=" << params[0] << " params[1]=" << params[1] << std::endl;
//    draw_channel(1, params, fval, "test_channel_1.png");
    TH1D *j1, *j2, *j3, *j4;
    draw(params, param_errs, fval, j1, j2, j3, j4, "minuit_output.*", ".", false);
    //draw(params, nullptr, "NOSAVE", fval, j1, j2, j3, j4, true);

    std::cin.get();

/*
//    params[0] = 1.0;
    params[1] = 0.0 + 0.001;
    std::cout << params[0] << " " << params[1] << std::endl;
    logLikelihood(n_params, nullptr, fval, params, 0);
    std::cout << "fval=" << fval << " for params[0]=" << params[0] << " params[1]=" << params[1] << std::endl;
//    draw_channel(1, params, fval, "test_channel_2.png");
    draw(params, nullptr, "NOSAVE", fval, j1, j2, j3, j4, false);
    //draw(params, nullptr, "NOSAVE", fval, j1, j2, j3, j4, true);

    std::cin.get();

    //params[0] = 1.0;
    params[1] = 0.296;
    std::cout << params[0] << " " << params[1] << std::endl;
    logLikelihood(n_params, nullptr, fval, params, 0);
    std::cout << "fval=" << fval << " for params[0]=" << params[0] << " params[1]=" << params[1] << std::endl;
//    draw_channel(1, params, fval, "test_channel_3.png");
    draw(params, nullptr, "NOSAVE", fval, j1, j2, j3, j4, false);
    //draw(params, nullptr, "NOSAVE", fval, j1, j2, j3, j4, true);

    std::cin.get();

//    params[0] = 0.8;
    params[1] = 0.0;
    std::cout << params[0] << " " << params[1] << std::endl;
    logLikelihood(n_params, nullptr, fval, params, 0);
    std::cout << "fval=" << fval << " for params[0]=" << params[0] << " params[1]=" << params[1] << std::endl;
//    draw_channel(1, params, fval, "test_channel_4.png");
    draw(params, nullptr, "NOSAVE", fval, j1, j2, j3, j4, false);
    //draw(params, nullptr, "NOSAVE", fval, j1, j2, j3, j4, true);

    std::cin.get();
    //return 0;
*/


#if 0
    std::ofstream os("parab.csv");
    const double gA_param_min = -0.5;
    const double gA_param_max = 3.0;
    const int nsteps = 20;
    for(int i = 0; i <= nsteps; ++ i)
    {
        int number_free_params = -1;
        double *CovMatrix = nullptr;

        double gA_param = gA_param_min + (gA_param_max - gA_param_min) * ((double)i / (double)(nsteps));
        AdjustActs[1] = gA_param;
//        std::cout << "gA=" << gA_param << std::endl;
//        std::cin.get();
        //fitBackgrounds_setparams(minuit, AdjustActs, AdjustActs_Err);
        TMinuit *new_minuit = fitBackgrounds(AdjustActs, AdjustActs_Err, CovMatrix, number_free_params, thePhase);


        const double width = 0.1;
        const double nnsteps = 20;
        int n_params = new_minuit->GetNumPars();
        double *params = new double[n_params];
        double *params_err = new double[n_params];
        for(int jx = 0; jx < n_params; ++ jx)
        {
            new_minuit->GetParameter(jx, params[jx], params_err[jx]);
        }
//        std::cout << "params[1]=" << params[1] << std::endl;
//        std::cin.get();
        const double min = params[0] - width;
        const double max = params[0] + width;
        os << "i=" << i << std::endl;
        for(int i = 0; i <= nnsteps; ++ i)
        {
            double fval = 0.0;
            double a = (double)i / (double)nnsteps;
            params[0] = min + (max - min) * a;
            logLikelihood(n_params, nullptr, fval, params, 0);
            os << params[0] << "," << params[1] << "," << fval << std::endl;
        }
    }
    return 0;
#endif


    if(1)
    {
        std::cout << "ready to test MPS" << std::endl;
        newloglikfitter_testmyphasespace(theParameterState, theFCN, FCN_min);
        std::cout << "MPS done" << std::endl;
    }
   

    if(0)
    {
        //newloglikfitter_gA_chisquaretest(minuit, AdjustActs, AdjustActs_Err);
    }


    ///////////////////////////////////////////////////////////////////////////

    print_adjustacts(std::cout);
    std::ofstream myFileFitResults("fit_results.txt", std::ios::out | std::ios::app);
    print_adjustacts(myFileFitResults);
    myFileFitResults.close();

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
        double param_init_value = 0.;
        double param_init_error = 0.; 
        get_paramInitValueError(thePhase, j, param_init_value, param_init_error);
        //std::cout << "value=" << param_init_value << " err=" << param_init_error << " AdjustActs[i]=" << AdjustActs[i] << std::endl;
        
        // 2020-06-17
        
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
        
        /*
        std::cout << i << " :\t" << AdjustActs[i] * param_init_value
                       << " +- " << AdjustActs_Err[i] * param_init_value;
        */
        // TODO: put the mutiplication by xi_31_init INSIDE the reweight/fit functions,
        // to restore uniformity in minuit parameters
        // xi works as an additive parameter so cannot do this, in case where HSD
        // is zero does not work


        Double_t change = 0.0;
        // 2020-06-17
        if(i != 1)
        {
            change = 100.0 * (AdjustActs[i] - 1.0);
        }
        else
        {
            // 2020-06-17
            //change = 100.0 * (AdjustActs[i] - xi_31_init);
            change = 100.0 * (AdjustActs[i] - xi_31_init_value);
            //change = 100.0 * (AdjustActs[i] - 1.0);
        }
        /*change = 100.0 * (AdjustActs[i] - 1.0);*/
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
        timestamp(summary_ofstream);
        
        //for(int i = 0; i < numberParams; ++ i)
        for(int i = 0; i < numberEnabledParams; ++ i)
        {
            int j = minuitParamNumberToParamNumberMap.at(i);

            double param_init_value;
            double param_init_error;
            get_paramInitValueError(thePhase, j, param_init_value, param_init_error);
            // TODO: summary_ofstream is printing nonsense for percentage change
            summary_ofstream << i << ","
                             << j << ","
                             << param_init_value << ","
                             << param_init_error << ","
                             << AdjustActs[i] << ","
                             << AdjustActs_Err[i] << ","
                             << (AdjustActs[i] - param_init_value) / param_init_error * 100.0 << " %" << std::endl;
            /*
            if(thePhase == 0)
            {
                // TODO: summary_ofstream is printing nonsense for percentage change
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
                // TODO: summary_ofstream is printing nonsense for percentage change
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
            */
        }
        summary_ofstream.close();
    }

    
    /*
    {
        std::string chisquarelog_fname("chisquarelog.txt");
        std::cout << "writing to " << chisquarelog_fname << std::endl;
        std::ofstream chisquarelog_ofstream(chisquarelog_fname.c_str(), std::ios::out | std::ios::app);
        chisquarelog_ofstream << std::endl;
        timestamp(chisquarelog_ofstream);
        int mo100bbnumber = paramNumberToMinuitParamNumberMap.at(10);
        //std::cout << "CHECK CHECK CHECK: mo100bbnumber=" << mo100bbnumber << std::endl;
        chisquarelog_ofstream << paramNameMap[mo100bbnumber] << ", " << AdjustActs[mo100bbnumber] << ", " << global_chisquare << ", (should be values for 100Mo bb)" << std::endl;
    }
    */


    std::ofstream of_numberofeventsafterfit("of_numberofeventsafterfit.txt", std::ofstream::out | std::ofstream::app);
    timestamp(of_numberofeventsafterfit);
    for(int i = 0; i < allMCSamples1D[0]->GetEntries(); ++ i)
    {
        TH1D *tmpHist = (TH1D*)allMCSamples1D[0]->At(i);
        Double_t integral = tmpHist->Integral();
        of_numberofeventsafterfit << tmpHist->GetName() << " number of events " << integral << std::endl;
    }
    of_numberofeventsafterfit.close();
 

    ///////////////////////////////////////////////////////////////////////////
    // write to output file
    ///////////////////////////////////////////////////////////////////////////

    /*
    TFile *fout = new TFile("Nd150_2e_P" + Phase + "_fit_histograms.root", "UPDATE");
    //fout->mkdir("1D");
    //fout->mkdir("2D");
    std::cout << "writing histograms into \"\"" << std::endl;

    std::ofstream ofoutaux("Nd150_2e_P" + Phase + "_fit_histograms.txt", std::ofstream::out);
    
    for(int channel = 0; channel < number1DHists; ++ channel)
    {
        std::cout << "writing: 1D: channel=" << channel << std::endl;
        
        TString channel_str;
        channel_str.Form("%i", channel);
        //fout->cd("/1D");
        fout->cd("/");
        TDirectory *dir = gDirectory;
        //TDirectory *dir_histogram = dir->mkdir("channel_" + channel_str);
        TDirectory *dir_histogram = dir->mkdir("channel_1D_" + channel_str);

        allDataSamples1D->At(channel)->Write();

        std::string auxname = "channel_1D_" + std::string(channel_str) + "/" + std::string(allDataSamples1D->At(channel)->GetName());
        ofoutaux << auxname << std::endl;

        for(int j = 0; j < allMCSamples1D[channel]->GetEntries(); ++ j)
        {
            allMCSamples1D[channel]->Write();

            std::string auxname = "channel_1D_" + std::string(channel_str) + "/" + std::string(allMCSamples1D[channel]->At(j)->GetName());
            ofoutaux << auxname << std::endl;
        }
    }

    for(int channel = 0; channel < number2DHists; ++ channel)
    {
        std::cout << "writing: 2D: channel=" << channel << std::endl;

        TString channel_str;
        channel_str.Form("%i", channel);
        //fout->cd("/2D");
        fout->cd("/");
        TDirectory *dir = gDirectory;
        //TDirectory *dir_histogram = dir->mkdir("channel_" + channel_str);
        TDirectory *dir_histogram = dir->mkdir("channel_2D_" + channel_str);
        
        allDataSamples2D->At(channel)->Write();

        std::string auxname = "channel_2D_" + std::string(channel_str) + "/" + std::string(allDataSamples2D->At(channel)->GetName());
        ofoutaux << auxname << std::endl;

        for(int j = 0; j < allMCSamples2D[channel]->GetEntries(); ++ j)
        {
            allMCSamples2D[channel]->Write();

            std::string auxname = "channel_2D_" + std::string(channel_str) + "/" + std::string(allMCSamples2D[channel]->At(j)->GetName());
            ofoutaux << auxname << std::endl;
        }
    }

    ofoutaux.close();
    */

// TODO
/*
    draw_all(AdjustActs, AdjustActs_Err, CovMatrix, number_free_params);
*/



    ///////////////////////////////////////////////////////////////////////////
    // Single Electron Energy Histogram
    // Chisquare Test
    // Plot chisquare as a function of the xi_31 parameter
    ///////////////////////////////////////////////////////////////////////////

    
    std::cout << "ok to run chisquare test? [press enter to continue]" << std::endl;
    //std::cin.get();

//    newloglikfitter_gA_chisquaretest(minuit, AdjustActs, AdjustActs_Err);


    // TODO: this isn't working because I am calling the wrong function
//    draw_outputdiff(AdjustActs, 0.0, "afterfit", 4);
//    AdjustActs[1] = 0.0; // HSD
//    draw_outputdiff(AdjustActs, 0.0, "beforefit", 4); // TODO: makes these functions easier to use
    // by splitting the function called into several different functions


//    draw(AdjustActs, AdjustActs_Err, "afterfit_hTotalE.*");
//    AdjustActs[1] = 0.0; // HSD // TODO: should it be zero or 1? does this indicate a bug elsewhere
    // BUG: I think I have interpreted AdjustActs as a scaling parameter and I cannot do that if my
    // initial xi value is 0! because 0 * anything = 0
//    draw(AdjustActs, AdjustActs_Err, "beforefit_hTotalE.*");

}
















