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
#include "newLogLikFitter_read_parameternames_lst.h"
#include "newLogLikFitter_reweight.h"
#include "newLogLikFitter_loglikelihood.h"
#include "newLogLikFitter_draw.h"
#include "newLogLikFitter_draw_outputdiff.h"
#include "newLogLikFitter_chisquaretest.h"


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

void book1DHistograms(Int_t channel_counter, TString theChannel,TString thePhase, TString theHistogram);
void book2DHistograms(Int_t channel_counter, TString theChannel,TString thePhase, TString theHistogram);

//void fitBackgrounds(double *AdjustActs, double *AdjustActs_Err, double*& CovMatrix, int& number_free_params, Int_t thePhase);
TMinuit * fitBackgrounds(double *AdjustActs, double *AdjustActs_Err, double*& CovMatrix, int& number_free_params, Int_t thePhase);






void timestamp(std::ofstream& os)
{
    time_t rawtime;
    time(&rawtime);
    struct tm *info;
    info = localtime(&rawtime);
    char timestringbuf[100];
    size_t ncount = strftime(timestringbuf, sizeof(timestringbuf), "%F %T", info);
    std::string datetimestamp(timestringbuf);
    os << datetimestamp << std::endl;
}


template<typename T, typename U>
void print_map(std::map<T, U> &map, const std::string& name)
{
    std::cout << "contents of map " << name << std::endl;
    for(auto it = map.cbegin(); it != map.cend(); ++ it)
    {
        std::cout << it->first << " -> " << it->second << std::endl;
    }
    
    /*
    std::cout << "contents of map paramNameToNumberMap:" << std::endl;
    for(auto it = paramNameToNumberMap.cbegin(); it != paramNameToNumberMap.cend(); ++ it)
    {
        std::cout << it->first << " -> " << it->second << std::endl;
    }
    */
}

void draw_inputdata()
{

    // debug
    TCanvas *c_nEqNull;
    c_nEqNull = new TCanvas("c_nEqNull", "c_nEqNull"); //, 4000, 3000);
    h_nEqNull->Draw("colz");
    c_nEqNull->SaveAs("c_nEqNull.png");
    c_nEqNull->SaveAs("c_nEqNull.pdf");
    c_nEqNull->SaveAs("c_nEqNull.C");
    delete c_nEqNull;

    // debug
    TCanvas *c_nEqTwo;
    c_nEqTwo = new TCanvas("c_nEqTwo", "c_nEqTwo"); //, 4000, 3000);
    h_nEqTwo->Draw("colz");
    c_nEqTwo->SaveAs("c_nEqTwo.png");
    c_nEqTwo->SaveAs("c_nEqTwo.pdf");
    c_nEqTwo->SaveAs("c_nEqTwo.C");
    delete c_nEqTwo;
}




void draw(const double* const AdjustActs, const double* const AdjustActs_Err, const double* const CovMatrix, const int number_free_params);


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
    const Double_t xi_31_baseline{0.296}; // TODO: this is WRONG change it
    // TODO: multiple instances, move to header file
    // TODO: also change in parameter list file
    //Double_t xi_31_init = 0.8;
    Double_t xi_31_init = 2.63e-01; //1.13; //0.296; // change to baseline value for testing purposes

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


    #if 0
    draw_inputdata();
    #endif

    std::cout << std::fixed << std::endl;

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);

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
        AdjustActs_Err[ix] = 1.0;
    }
    // TODO: fix this
    AdjustActs[1] = xi_31_init;


    /*
    Int_t number_free_params = minuit->GetNumFreePars();
    Double_t *CovMatrix = new Double_t[number_free_params * number_free_params];
    for(int ix{0}; ix < number_free_params * number_free_params; ++ ix)
    {
        CovMatrix[ix] = 0.;
    }
    */
    


    // draw some gA values as output
    const int i_max = 100;
    for(int i = 0; i < i_max; ++ i)
    {

        Double_t xi_31_default = 0.296;
        Double_t xi_31_half_range = 5.0;

        Double_t xi_31_offset = 0.05;
        Double_t xi_31_min = xi_31_default - xi_31_half_range + xi_31_offset;
        Double_t xi_31_max = xi_31_default + xi_31_half_range + xi_31_offset;
        Double_t xi_31_value = ((double)i / (double)i_max) * (xi_31_max - xi_31_min) + xi_31_min;
        AdjustActs[1] = xi_31_value;

        TString xi_31_str;
        xi_31_str.Form("%f", xi_31_value);

        // TODO, put in custom directory with text file containing params
        draw(AdjustActs, AdjustActs_Err, std::string("hTotalE_") + std::string(xi_31_str) + std::string(".png"));
        draw_2D(AdjustActs, AdjustActs_Err, std::string("hHighLowEnergy_") + std::string(xi_31_str) + std::string(".png"));
        draw_outputdiff(AdjustActs, 0.296, std::string("houtputdiff_") + std::string(xi_31_str) + std::string(".png"));
    }

    std::cout << "done, check output folder for figures" << std::endl;
    return 0;



    int number_free_params = -1;
    double *CovMatrix = nullptr;

    //fitBackgrounds(AdjustActs, AdjustActs_Err, CovMatrix, number_free_params, thePhase);
    // fit of backgrounds disabled.
    // need to fit Nd150 parameter only
    // disabled in fitBackgrounds function
    
    TMinuit *minuit = fitBackgrounds(AdjustActs, AdjustActs_Err, CovMatrix, number_free_params, thePhase);


    #if 1
    if(0)
    {
        newloglikfitter_gA_chisquaretest(minuit, AdjustActs, AdjustActs_Err);
    }
    #endif


    ///////////////////////////////////////////////////////////////////////////

    std::cout << "The following adjustments (in minuit parameter units) should be made:" << std::endl;
    std::cout << "Note that gA (1) is a special parameter" << std::endl;
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
        timestamp(summary_ofstream);
        
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
        summary_ofstream.close();
    }

    
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


    draw(AdjustActs, AdjustActs_Err, CovMatrix, number_free_params);




}

void draw(const double* const AdjustActs, const double* const AdjustActs_Err, const double* const CovMatrix, const int number_free_params)
{
    // TODO: load AdjustActs, AdjustActs_Err from file

    ///////////////////////////////////////////////////////////////////////////
    // read from input file
    ///////////////////////////////////////////////////////////////////////////

    if(allDataSamples1D->GetEntries() == 0)
    {

        TFile *fin = new TFile("Nd150_2e_P" + Phase + "_fit_histograms.root", "UPDATE");
        std::cout << "reading histograms from \"\"" << std::endl;

        std::ifstream ifinaux("Nd150_2e_P" + Phase + "_fit_histograms.txt", std::ifstream::in);

        std::string auxname;

        for(int channel = 0; channel < number1DHists; ++ channel)
        {
            std::cout << "reading: 1D: channel=" << channel << std::endl;
            
            TString channel_str;
            channel_str.Form("%i", channel);
            //fout->cd("/1D");
            //fout->cd("/");
            //TDirectory *dir = gDirectory;
            //TDirectory *dir_histogram = dir->mkdir("channel_" + channel_str);
            //TDirectory *dir_histogram = dir->mkdir("channel_1D_" + channel_str);

            ifinaux >> auxname;

            TString hfullname = TString(auxname);
            //TString hfullname = "channel_1D_" + channel_str + "/" + hname;

            allDataSamples1D->Add((TH1F*)fin->Get(hfullname));
            allMCSamples1D[channel]->Add((TH1F*)fin->Get(hfullname));

        }

        for(int channel = 0; channel < number2DHists; ++ channel)
        {
            std::cout << "reading: 2D: channel=" << channel << std::endl;

            TString channel_str;
            channel_str.Form("%i", channel);
            //fout->cd("/2D");
            //fout->cd("/");
            //TDirectory *dir = gDirectory;
            //TDirectory *dir_histogram = dir->mkdir("channel_" + channel_str);
            //TDirectory *dir_histogram = dir->mkdir("channel_2D_" + channel_str);
            
            ifinaux >> auxname;

            TString hfullname = TString(auxname);
            //TString hfullname = "channel_2D_" + channel_str + "/" + hname;

            allDataSamples2D->Add((TH2F*)fin->Get(hfullname));
            allMCSamples2D[channel]->Add((TH2F*)fin->Get(hfullname));

        }

        ifinaux.close();

    }

    draw(AdjustActs, AdjustActs_Err, "hTotalE.*");
    draw_2D(AdjustActs, AdjustActs_Err, "hHighLowEnergy.*");
    draw_covariance_matrix(CovMatrix, number_free_params, "cov_matrix.*");


    draw_outputdiff(AdjustActs, 0.296, "c_outputdiff_.png");

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

                //std::string directory("scaled/hTotalE_/");
                std::string directory("scaled/" + theHistogram + "/");
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
                            //Rn222BkgFiles);//,
                            Rn222BkgFilesNew);//,
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

    std::cout << "Data" << std::endl;
    // TODO here
    // what is name in other section of code
    //std::string name(theHistogram + "data_2e");
    //std::string directory("processeddata/hTotalE_/");
    std::string directory("processeddata/" + theHistogram + "/");
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
// book2DHistograms
///////////////////////////////////////////////////////////////////////////////

void book2DHistograms_helper(
    TFile *myFile,
    Int_t channel_counter,
    TString theChannel,
    TString thePhase_arg,
    TString theHistogram,
    const int nBkgs,
    TString *BkgFiles)
    //, TH1F *tmpHist)
{
        
    TH2F *tmpHist = nullptr;

    for(int i = 0; i < nBkgs; i++)
    {
        // check if parameter is enabled
        // convert parameter string name to index

        // convert TString to std::string
        std::string mc_name = std::string(BkgFiles[i]);

        // example: "bi214_int_rot" -> "bi214_int_rot,pb214_int_rot"
        std::string search_object = MCNameToParamNameMap.at(mc_name);

        // example: "bi214_int_rot,pb214_int_rot" -> 4
        // check if parameter number exists
        // (was defined by parameter_names.lst)
        if(paramNameToNumberMap.count(search_object) > 0)
        {
            // convert from mc sample name to param number
            int param_number = paramNameToNumberMap.at(search_object);

            // check if this parameter number is enabled
            if(std::find(enabled_params.begin(), enabled_params.end(), param_number) != enabled_params.end())
            {

                // TODO:
                //std::string directory("scaled/hHighLowEnergy_/");
                std::string directory("scaled/" + theHistogram + "/");
                std::string name(theHistogram + BkgFiles[i] + "_fit_scaled");
                std::string fullname = directory + name;
                std::string new_name(theHistogram + BkgFiles[i] + "_fit");
                std::cout << "fullname=" << fullname << std::endl;

                // TODO: try catch block
                // load sample
                tmpHist = (TH2F*)myFile->Get(fullname.c_str())->Clone(new_name.c_str());

                if(tmpHist != nullptr)
                {
                    // scale by activity

                    // convert parameter number to minuit parameter number
                    //minuit_param_number = paramNumberToMinuitParamNumberMap.at(param_number);

                    // TODO: change such that samples are pre-scaled by activity input value
                    // get initial parameter values and error
                    Double_t param_init_value = 0.;
                    Double_t param_init_error = 0.; 
                    if(thePhase == 0)
                    {
                        param_init_value = paramInitValueP1Map[param_number];
                        param_init_error = paramInitErrorP1Map[param_number];
                    }
                    else if(thePhase == 1)
                    {
                        param_init_value = paramInitValueP2Map[param_number];
                        param_init_error = paramInitErrorP2Map[param_number];
                    }
                    else
                    {
                        std::cout << "ERROR: Invalid value for thePhase: thePhase=" << thePhase << std::endl;
                    }
                    Double_t scale_factor = param_init_value;

                    // account for 208 Tl branching ratio of 36 %
                    // TODO: should I move this into fit_2e code
                    // and apply using ->Fill() function call with
                    // weight = 0.36
                    if(mc_name == std::string("tl208_int_rot") ||
                       mc_name == std::string("tl208_feShield") ||
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


                    allMCSamples2D[channel_counter]->Add(tmpHist);
                    // TODO: does this work as expected for secular equlibrium samples?

                    //std::cout << tmpHist->GetName() << std::endl;

                }
                else
                {
                    std::cout << "could not find histogram in file: " << fullname << " - disabling parameter number " << param_number << std::endl;
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

            print_map(paramNameToNumberMap, "paramNameToNumberMap");
            print_map(MCNameToParamNameMap, "MCNameToParamNameMap");
        }
    }
}

// channel_counter = 0
// theChannel = "2e_"
// thePhase = "P1"
// theHistogram = "hHighLowEnergy_"
void book2DHistograms(Int_t channel_counter, TString theChannel, TString thePhase_arg, TString theHistogram) {

    std::cout << "booking 2D hists for " << theChannel << " " << thePhase_arg << std::endl;
    allMCSamples2D[channel_counter] = new TObjArray();

    TFile *aFile = TFile::Open("Nd150_" + theChannel + thePhase_arg + ".root");
    
    std::cout << "External" << std::endl;
    book2DHistograms_helper(aFile, channel_counter, theChannel,
                            thePhase_arg, theHistogram,
                            nExternalBkgs,
                            ExternalBkgFiles);

    std::cout << "Internal" << std::endl;
    book2DHistograms_helper(aFile, channel_counter, theChannel,
                            thePhase_arg, theHistogram,
                            nInternalBkgs,
                            InternalBkgFiles);

    std::cout << "Rn 222" << std::endl;
    book2DHistograms_helper(aFile, channel_counter, theChannel,
                            thePhase_arg, theHistogram,
                            nRn222Bkgs,
                            //Rn222BkgFiles);//,
                            Rn222BkgFilesNew);

    std::cout << "Rn 220" << std::endl;
    book2DHistograms_helper(aFile, channel_counter, theChannel,
                            thePhase_arg, theHistogram,
                            nRn220Bkgs,
                            Rn220BkgFiles);

    std::cout << "Neighbour" << std::endl;
    book2DHistograms_helper(aFile, channel_counter, theChannel,
                            thePhase_arg, theHistogram,
                            nNeighbours,
                            NeighbourFiles);

    std::cout << "Nd150" << std::endl;
    book2DHistograms_helper(aFile, channel_counter, theChannel,
                            thePhase_arg, theHistogram,
                            nNd150Samples,
                            Nd150Files);

    // TODO here
    // what is name in other section of code
    //std::string name(theHistogram + "data_2e");
    std::string directory("processeddata/" + theHistogram + "/");
    std::string name(theHistogram + "data_2e");
    std::string fullname = directory + name;
    std::cout << "fullname=" << fullname << std::endl;

    // load histogram from file
    // TODO: try catch block
    TH1F *tmpHist = (TH1F*)aFile->Get(fullname.c_str())->Clone();
    if(tmpHist != nullptr)
    {
        //TH1F *tmpHist = nullptr;
        // 2020-04-03: removed changing of histogram name
        //std::string hist_name("data_" + theChannel + thePhase_arg);
        //std::cout << "Get() : " << name << " from file, Clone() : " << hist_name << std::endl;
        //tmpHist = (TH1F*)gDirectory->Get(name.c_str())->Clone(hist_name.c_str());
        //tmpHist = (TH1F*)gDirectory->Get(fullname.c_str())->Clone();
        allDataSamples2D->Add((TH1F*)tmpHist);
    }
    else
    {
        std::cout << "gDirectory->GetListOfKeys() does not contain " << fullname << std::endl;
    }

    // aFile->Close();
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


    // Then get results
    //for(int i = 0; i < numberParams; i++)
    for(int i = 0; i < numberEnabledParams; i++)
    {
        minuit->GetParameter(i, AdjustActs[i], AdjustActs_Err[i]);
    }

   
#if 1
    if(0)
    {
        newloglikfitter_100Mo_chisquaretest(minuit, AdjustActs, AdjustActs_Err);
    }
#endif


    


    // Some people have seen further by standing on the shoulders of giants.
    // In my case, my vision has been obscured by floundering hopeless idiots
    // doing about as much as it was possible to do to inhibit my ability to
    // see anything, either simply by being incompetent, or by actively
    // doing everything possible to throw as many spanners into the works as
    // time would allow for.

#if 1
    if(0)
    {
            newloglikfitter_testmyphasespace(minuit, AdjustActs, AdjustActs_Err);
    }
#endif



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

  
    return minuit;

}





