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


// note: these must appear in correct order and after general includes above
#include "newLogLikFitter.h"
#include "newLogLikFitter_aux.h"
#include "newLogLikFitter_print.h"
#include "newLogLikFitter_read_parameternames_lst.h"
#include "newLogLikFitter_book1DHistograms.h"
#include "newLogLikFitter_book2DHistograms.h"
#include "newLogLikFitter_reweight.h"
#include "newLogLikFitter_loglikelihood.h"
#include "newLogLikFitter_draw.h"
#include "newLogLikFitter_draw_2D.h"
#include "newLogLikFitter_draw_outputdiff.h"
#include "newLogLikFitter_draw_all.h"
#include "newLogLikFitter_chisquaretest.h"
#include "newLogLikFitter_test.h"


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

    allDataSamples1D = new TObjArray();
    allDataSamples2D = new TObjArray();



    

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
        

        do_test_xi_31_test1(AdjustActs_copy, AdjustActs_Err_copy);
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
    TMinuit *minuit = fitBackgrounds(AdjustActs, AdjustActs_Err, CovMatrix, number_free_params, thePhase);


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
        newloglikfitter_testmyphasespace(minuit, AdjustActs, AdjustActs_Err);
        std::cout << "MPS done" << std::endl;
    }
   

    if(0)
    {
        newloglikfitter_gA_chisquaretest(minuit, AdjustActs, AdjustActs_Err);
    }


    ///////////////////////////////////////////////////////////////////////////

    std::cout << "The following adjustments (in minuit parameter units) should be made:" << std::endl;
    std::cout << "Note that gA (1) is a special parameter" << std::endl;
    // TODO: after homogenizing gA is no longer a special parameter
    //for(int i = 0; i < numberParams; i++)
    std::ofstream myFileFitResults("fit_results.txt", std::ios::out | std::ios::app);
    timestamp(myFileFitResults);
    for(int i = 0; i < numberEnabledParams; i++)
    {
        std::cout << i << " :\t" << AdjustActs[i] << " +- " << AdjustActs_Err[i] << std::endl;  
        myFileFitResults << AdjustActs[i] << " +- " << AdjustActs_Err[i] << std::endl;
    }
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
        TH1F *tmpHist = (TH1F*)allMCSamples1D[0]->At(i);
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















///////////////////////////////////////////////////////////////////////////////
// fitBackgrounds
///////////////////////////////////////////////////////////////////////////////



TMinuit* fitBackgrounds_init(double *AdjustActs, double *AdjustActs_Err)
{
//    std::cout << " 1 " << AdjustActs[1] << std::endl;
//    std::cin.get();

    std::cout << ">>>>> fitBackgrounds_init()" << std::endl;

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
    minuit->SetPrintLevel(3);

    // moved reading of parameter list file to before book histogram function
    // calls
  
    // TODO: printing sample_names array, what does this do, is it still used?
    print_paramNameMap();

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


        if(std::find(fixed_params.begin(), fixed_params.end(), i) != fixed_params.end())
        {
            // define parameter using constrained value if hard constrained

            std::cout << "minuit: fixed parameter i=" << i << std::endl;
            
            //minuit->DefineParameter(i, "_" + i_str + "_", 1.0, 0.1, 0.0, 2.0);
            //minuit->DefineParameter(i, "_" + i_str + "_FIXED", param_init_value, param_init_error, -2.0, 1.0e+06);
            //minuit->FixParameter(i);
            
            //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_FIXED", param_init_value, param_init_error);
            // TODO: change such that samples are pre-scaled by activity input value
            
            // 2020-06-17
            if(i == 1)
            {
//                std::cout << "i 1 fixed " << AdjustActs[i] << std::endl;
//                std::cin.get();
                minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_FIXED", AdjustActs[i], 0.5, 0.0, 50.0);
            }
            else
            {
                minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_FIXED", 1.0, 0.5, 0.0, 50.0);
            }
            /*
            if(i == 1)
            {
                //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_FIXED", AdjustActs[i], 0.5, 0.0, 50.0);
                minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_FIXED", 1.0, 0.5, 0.0, 50.0);
            }
            else
            {
                minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_FIXED", 1.0, 0.5, 0.0, 50.0);
            }
            */

            //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_FIXED", 1.0, 0.5, 0.0, 50.0);
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

            // 2020-06-17
            
            if(i == 1)
            {
//                std::cout << "i 1 not fixed " << AdjustActs[i] << std::endl;
//                std::cin.get();
                // TODO: fix this
                // does not work if xi_31 paramter is not number 1
                //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", xi_31_init, 0.5, 0.0, 1000.0);
                //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", AdjustActs[i], 0.5, 0.0, 1000.0);
                minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", AdjustActs[i], 0.5, -1.0, 5.0);
            }
            else
            {
                minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", 1.0, 0.5, 0.0, 1000.0);
            }
            
            /*
            if(i == 1)
            {
                // TODO: fix this
                // does not work if xi_31 paramter is not number 1
                //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", xi_31_init, 0.5, 0.0, 1000.0);
                //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", AdjustActs[i], 0.5, 0.0, 1000.0);
                minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", 1.0, 0.5, -1.0, 5.0);
            }
            else
            {
                minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", 1.0, 0.5, 0.0, 1000.0);
            }
            */


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
    std::cout << "Ready to exec fix" << std::endl;

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
    //minuit->mnexcm("SET EPS", 0.01);
    //minuit->SetEPS(1.0e-3); // TODO
    //give it the function to use
    minuit->SetFCN(logLikelihood);
    //std::cout << "calling: minuit->mnsimp()" << std::endl;
    //minuit->mnsimp();
    std::cout << "calling: minuit->Migrad()" << std::endl;


// draw 1D histograms, write chisquare values to file for range of
// different gA values
//    newloglikfitter_gA_chisquaretest(minuit, AdjustActs, AdjustActs_Err);

    /*
    Int_t npar = 0;
    Double_t fval;
    logLikelihood(npar, nullptr, fval, AdjustActs, 0);
    std::cout << "fval=" << fval << std::endl;
    draw(AdjustActs, AdjustActs_Err, CovMatrix, number_free_params);
    std::cin.get();

    AdjustActs[1] = 0.296 -0.1;
    logLikelihood(npar, nullptr, fval, AdjustActs, 0);
    std::cout << "fval=" << fval << std::endl;
    draw(AdjustActs, AdjustActs_Err, CovMatrix, number_free_params);
    std::cin.get();

    AdjustActs[1] = 0.296 + 0.1;
    logLikelihood(npar, nullptr, fval, AdjustActs, 0);
    std::cout << "fval=" << fval << std::endl;
    draw(AdjustActs, AdjustActs_Err, CovMatrix, number_free_params);
    std::cin.get();
    */


    return minuit;
}



void fitBackgrounds_setparams(TMinuit *minuit, double *AdjustActs, double *AdjustActs_Err)
{

}


void fitBackgrounds_exec(TMinuit *minuit, double *AdjustActs, double *AdjustActs_Err)
{
//                std::cout << "i 1 exec: " << AdjustActs[i] << std::endl;
//                std::cin.get();

    // TODO: re-enable
    minuit->Migrad();

    // Then get results
    //for(int i = 0; i < numberParams; i++)
    for(int i = 0; i < numberEnabledParams; i++)
    {
        minuit->GetParameter(i, AdjustActs[i], AdjustActs_Err[i]);
    }
 }



void fitBackgrounds_postexectest(TMinuit *minuit, double *AdjustActs, double *AdjustActs_Err)
{


    // TODO: remove AdjustActs, AdjustActs_Err arguments?
    // Do I want these to be copy of original values

    if(0)
    {
        newloglikfitter_100Mo_chisquaretest(minuit, AdjustActs, AdjustActs_Err);
    }

    // Some people have seen further by standing on the shoulders of giants.
    // In my case, my vision has been obscured by floundering hopeless idiots
    // doing about as much as it was possible to do to inhibit my ability to
    // see anything, either simply by being incompetent, or by actively
    // doing everything possible to throw as many spanners into the works as
    // time would allow for.

    if(0)
    {
        newloglikfitter_testmyphasespace(minuit, AdjustActs, AdjustActs_Err);
    }


}


void fitBackgrounds_getcovmatrix(TMinuit* minuit, double *&CovMatrix, int& number_free_params)
{

    // TODO: this no longer works, or does it?
    // needs to take into account the number of ENABLED free params
    // NOTE: 2020-04-16 fixed
    number_free_params = minuit->GetNumFreePars();
    CovMatrix = new Double_t[number_free_params * number_free_params];
    for(int ix{0}; ix < number_free_params * number_free_params; ++ ix)
    {
        CovMatrix[ix] = 0.;
    }
    minuit->mnemat(CovMatrix, number_free_params);  

}




//void fitBackgrounds(double *AdjustActs, double *AdjustActs_Err, double *&CovMatrix, int& number_free_params, Int_t thePhase)
TMinuit* fitBackgrounds(double *AdjustActs, double *AdjustActs_Err, double *&CovMatrix, int& number_free_params, Int_t thePhase)
{

    std::cout << ">>>>> fitBackgrounds()" << std::endl;
    
    TMinuit* minuit = fitBackgrounds_init(AdjustActs, AdjustActs_Err);
    //fitBackgrounds_exec(minuit, AdjustActs, AdjustActs_Err);
    //fitBackgrounds_postexectest(minuit, AdjustActs, AdjustActs_Err);
    //fitBackgrounds_getcovmatrix(minuit, CovMatrix, number_free_params);

    return minuit;

}



