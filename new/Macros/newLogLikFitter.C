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
#include "newLogLikFitter_print.h"
#include "parametergroup.h"
#include "newLogLikFitter.h"
#include "newLogLikFitter_aux.h"
#include "newLogLikFitter_printfitresult.h"
#include "newLogLikFitter_read_parameternames_lst.h"
#include "Systematics.h"
#include "newLogLikFitter_reweight.h"
#include "newLogLikFitter_reweight_apply.h"
#include "newLogLikFitter_reweight_apply_fakedata.h"
#include "newLogLikFitter_reweight_apply_data.h"
#include "newLogLikFitter_buildfakedata.h"
#include "newLogLikFitter_book1DHistograms.h"
#include "newLogLikFitter_book2DHistograms.h"
#include "newLogLikFitter_stacker_helper.h"
#include "newLogLikFitter_logpoisson.h"
#include "newLogLikFitter_drawchannel.h"
#include "newLogLikFitter_draw.h"
#include "newLogLikFitter_draw_2D.h"
#include "newLogLikFitter_draw_outputdiff.h"
#include "newLogLikFitter_draw_all.h"
#include "newLogLikFitter_rebuild_150Nd_MC.h"
#include "newLogLikFitter_rebuild_150Nd_data.h"
#include "newLogLikFitter_rebuild_150Nd_fakedata.h"
#include "MinimizeFCNAxialVector.h"
#include "newLogLikFitter_fitBackgrounds.h"
#include "newLogLikFitter_chisquaretest.h"
#include "newLogLikFitter_test.h"




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





                //
                // if n is large, then exp(-n) may fail in Poisson calc
                // n! may also fail, so use Stirling
                //
                // so calculate LL from expansion of log
                //
                //if(nMC > 10.0)
                //{
                //    // nMC > 1.0e-05 here already
                //    const double l = nMC;
                //    const double n = nData;
                //    const double stirling = n*TMath::Log(n) - n;
                //    ll_channel += -l + n*TMath::Log(l) - stirling;
                //}



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
void loadFiles(int);


//void fitBackgrounds(double *AdjustActs, double *AdjustActs_Err, double*& CovMatrix, int& number_free_params, Int_t thePhase);
TMinuit * fitBackgrounds(double *AdjustActs, double *AdjustActs_Err, double*& CovMatrix, int& number_free_params, Int_t thePhase);









bool load_from_script(
    int i,
    int &number_job_id,
    std::string &output_name,
    int &start_index,
    int &stop_index
    )
{
    int script_index = 0;
    for(;; ++ script_index)
    {
        if(script_index != i) continue;

        std::string script_fname = "script" + std::to_string(script_index) + ".txt";
        std::cout << "loading... " << script_fname << std::endl;
        std::ifstream ifs(script_fname);
        if(ifs.is_open())
        {
            for(int lineindex = 1; ifs.good(); ++ lineindex)
            {
                std::string ifs_line;
                std::getline(ifs, ifs_line);
                if(lineindex == 1)
                {
                    std::string delim = "=";
                    std::string token1 = ifs_line.substr(0, ifs_line.find(delim));
                    std::string token2 = ifs_line.substr(ifs_line.find(delim) + 1);
                    if(token1 == "NUMBER")
                    {
                        number_job_id = std::stoi(token2);
                    }
                    else
                    {
                        return false;
                    }
                }
                else if(lineindex == 2)
                {
                    std::string delim = "=";
                    std::string token1 = ifs_line.substr(0, ifs_line.find(delim));
                    std::string token2 = ifs_line.substr(ifs_line.find(delim) + 1);
                    if(token1 == "OUTPUT_NAME")
                    {
                        output_name = token2;
                    }
                    else
                    {
                        return false;
                    }
                }
                else if(lineindex == 3)
                {
                    std::string delim = "=";
                    std::string token1 = ifs_line.substr(0, ifs_line.find(delim));
                    std::string token2 = ifs_line.substr(ifs_line.find(delim) + 1);
                    if(token1 == "START_INDEX")
                    {
                        start_index = std::stoi(token2);
                    }
                    else
                    {
                        return false;
                    }
                }
                else if(lineindex == 4)
                {
                    std::string delim = "=";
                    std::string token1 = ifs_line.substr(0, ifs_line.find(delim));
                    std::string token2 = ifs_line.substr(ifs_line.find(delim) + 1);
                    if(token1 == "STOP_INDEX")
                    {
                        stop_index = std::stoi(token2);
                    }
                    else
                    {
                        return false;
                    }
                }
                else if(lineindex == 5)
                {
                    std::string delim = "=";
                    std::string token1 = ifs_line.substr(0, ifs_line.find(delim));
                    std::string token2 = ifs_line.substr(ifs_line.find(delim) + 1);
                    if(token1 == "RUNNING")
                    {
                        if(token2 == "false")
                        {
                            ifs.close();
                            std::ofstream ofs(script_fname, std::ios::out);
                            ofs << "NUMBER=" << number_job_id << std::endl;
                            ofs << "OUTPUT_NAME=" << output_name << std::endl;
                            ofs << "START_INDEX=" << start_index << std::endl;
                            ofs << "STOP_INDEX=" << stop_index << std::endl;
                            ofs << "RUNNING=true" << std::endl;
                            return true;
                        }
                        else if(token2 == "true")
                        {
                            break;
                        }
                        else
                        {
                            return false;
                        }
                    }
                    else
                    {
                        return false;
                    }
                }
            }
        }
        else
        {
            return false;
        }
    }
}





void newLogLikFitter(int i)
{
    loadFiles(i);
}


void loadFiles(int i)
{

    TH1::AddDirectory(false);
    TH2::AddDirectory(false);


    ///////////////////////////////////////////////////////////////////////////
    // parallel mode code
    ///////////////////////////////////////////////////////////////////////////

    int number_job_id;
    std::string output_name;
    int start_index;
    int stop_index;
    if(MODE_PARALLEL == 1)
    {
        bool success = load_from_script(i, number_job_id, output_name, start_index, stop_index);
        if(success == true)
        {
            std::cout << "Job Init: ID=" << number_job_id << std::endl;
            std::cout << "output_name=" << output_name << std::endl;
            std::cout << "START_INDEX=" << start_index << std::endl;
            std::cout << "STOP_INDEX=" << stop_index << std::endl;
            // do nothing = continue
        }
        else
        {
            std::cout << "fail" << std::endl;
            return;
        }
    }
    else
    {
        if(i != 0) return;
        number_job_id = 0;
        output_name = "noparallel";
        start_index = 0;
        stop_index = 51;
    }


    ///////////////////////////////////////////////////////////////////////////
    // read parameter file list
    ///////////////////////////////////////////////////////////////////////////

    // read parameter_name.lst file
    // read parameter list file
    read_parameter_list_file();
    g_pg.print();


    ///////////////////////////////////////////////////////////////////////////
    // load spectral data
    ///////////////////////////////////////////////////////////////////////////

    std::cout << std::scientific;

    // load data
    std::cout << "attempting to load spectral data from file" << std::endl;

    //bb_Q = 3.368;
    bb_Q = 3.37138;
    // ramdisk, should be faster?
    std::size_t count_G0 =
        LD_REWEIGHT_DATA_2(h_nEqNull,
                           //"/home/ecb/150Nd/150Nd-data/dG150Nd/G0/dG0.dat",
                           "/home/ecb/100Mo-150Nd/gA_theoretical_files/psf-nuclei/150Nd/0-N0/nEqNull.dat",
                           "h_nEqNull",
                           "nEqNull"//,
                           //0.0, bb_Q
                           );
    std::size_t count_G2 =
        LD_REWEIGHT_DATA_2(h_nEqTwo,
                           //"/home/ecb/150Nd/150Nd-data/dG150Nd/G2/dG2.dat",
                           "/home/ecb/100Mo-150Nd/gA_theoretical_files/psf-nuclei/150Nd/1-N2/nEqTwo.dat",
                           "h_nEqTwo",
                           "nEqTwo"//,  0.0, bb_Q
                           );
                           // TODO: change data file location!

    // phase space integrals
    // TODO: these were the wrong way around!
    // I need to re-run this code, plotting the MPS without fit
    const Double_t G0_ps_integral_MeV = 0.759721E-45; //0.744684E-45;
    const Double_t G0_ps_integral_yrinv = 0.364244E-16; //0.357034E-16;
    const Double_t G2_ps_integral_MeV = 0.429791E-45; //0.420438E-45;
    const Double_t G2_ps_integral_yrinv = 0.206061E-16; //0.201577E-16;

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
//    last_xi_31_parameter_value = 0.0;




    ///*const Double_t*/ bb_Q = 3.368;
    double count = 0;
    if(count_G0 == count_G2)
    {
        count = count_G0;
    }
    else
    {
        std::cout << "error: count_G0=" << count_G0 << ", count_G2=" << count_G2 << std::endl;
        return;
    }

    psiN0 = G0_ps_integral_MeV;
    psiN2 = G2_ps_integral_MeV; // TODO: check this is the correct option
    
    std::cout << "histogram format constructed" << std::endl;


    if(0)
    {
        draw_inputdata();
    }

    std::cout << std::fixed << std::endl;


    ///////////////////////////////////////////////////////////////////////////
    // general ROOT init
    ///////////////////////////////////////////////////////////////////////////


    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    //gStyle->SetPalette(kBird);
    //gStyle->SetPalette(kBrownCyan);
    gStyle->SetPalette(kLightTemperature);
    //gStyle->SetNumberContours(1000);

    // TODO: C++ syntax, should be TObjArray?
    // without parens
    allDataSamples1D = new TObjArray();
    allDataSamples2D = new TObjArray();
    // TODO could call build_fake_data() in this function after load?
    //allFakeDataSamples1D = nullptr;
    //allFakeDataSamples2D = nullptr;
    allFakeDataSamples1D = new TObjArray();
    allFakeDataSamples2D = new TObjArray();


    gEnablePhase1 = true;
    gEnablePhase2 = true;
    //parameter_group pg;
    //g_pg = pg;
    // working here

    


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


    gSystematics.systematic_energy_offset = -0.1;
    const int xi_31_ext_param_number = g_pg.get_xi_31_ext_param_number();
        //if(param[xi_31_ext_param_number] != g_pg.file_params.at(xi_31_ext_param_number).paramLastValue)
    //const double xi_31{param[xi_31_ext_param_number]};
    const double xi_31{g_pg.file_params.at(xi_31_ext_param_number).paramInitValue};
    std::cout << "xi_31=" << xi_31 << std::endl;
//    std::cin.get();
    // NOTE: you put whatever value of xi_31 you like in here
    // to construct fakedata with that xi_31 value
    const double xi_31_SSD = 0.296;
    //rebuild_fake_data_systematics(xi_31, xi_31_baseline);
    rebuild_fake_data_systematics(0.0, xi_31_baseline); // want to check if the fitter can fit itself to itsel
    //rebuild_fake_data_systematics(xi_31_SSD, xi_31_baseline); // want to check if the fitter can fit itself to itsel


    // 1d: Phase 1 & 2
    for(int channel = 0; channel < number1DHists; ++ channel)
    {
        //book1DHistograms(0, "2e_", "hTotalE_");
        std::cout << "book1DHistograms(" << channel << ", 2e_, " << channel_histname_1D[channel] << ")" << std::endl;
        book1DHistograms(channel, "2e_", channel_histname_1D[channel]);
    }
    //book1DHistograms(1, "2e_", "hSingleEnergy_");
    //book1DHistograms(2, "2e_", "hHighEnergy_");
    //book1DHistograms(3, "2e_", "hLowEnergy_");
    //book1DHistograms(4, "2e_", "hEnergySum_");
    //book1DHistograms(5, "2e_", "hEnergyDiff_");
    //map_1d_channel_to_phase[0] = 0;
    //map_1d_channel_to_phase[1] = 0;
    //map_1d_channel_to_phase[2] = 0;
    //map_1d_channel_to_phase[3] = 0;
    //map_1d_channel_to_phase[4] = 0;
    //map_1d_channel_to_phase[5] = 0;

    // 2d: Phase 1 & 2
    for(int channel = 0; channel < number2DHists; ++ channel)
    {
        std::cout << "book2DHistograms(" << channel << ", 2e_, " << channel_histname_2D[channel] << ")" << std::endl;
        book2DHistograms(channel, "2e_", channel_histname_2D[channel]);
    }
#if 0
    book1DHistograms( 6, "2e_", "hTotalE_");
    book1DHistograms( 7, "2e_", "hSingleEnergy_");
    book1DHistograms( 8, "2e_", "hHighEnergy_");
    book1DHistograms( 9, "2e_", "hLowEnergy_");
    book1DHistograms(10, "2e_", "hEnergySum_");
    book1DHistograms(11, "2e_", "hEnergyDiff_");
    map_1d_channel_to_phase[6] = 1;
    map_1d_channel_to_phase[7] = 1;
    map_1d_channel_to_phase[8] = 1;
    map_1d_channel_to_phase[9] = 1;
    map_1d_channel_to_phase[10] = 1;
    map_1d_channel_to_phase[11] = 1;
#endif

    std::cout << "All histograms loaded" << std::endl;



    // First, create a root file to hold all of the histograms
    //TFile *myFile = TFile::Open("Nd150_loglikResults.root", "RECREATE");

#if 0
    // 2d: Phase 1
    book2DHistograms(0, "2e_", "P2", "hHighLowEnergy_");
    map_2d_channel_to_phase[0] = 0;

    // 2d: Phase 2
    book2DHistograms(1, "2e_", "P2", "hHighLowEnergy_");
    map_2d_channel_to_phase[1] = 1;
#endif
    // Array to hold activity adjustment parameters
    //Double_t AdjustActs[numberParams];
    //Double_t AdjustActs_Err[numberParams];
    //for(Int_t ix{0}; ix < numberParams; ++ ix)
    /*
    Double_t AdjustActs[numberEnabledParams];
    Double_t AdjustActs_Err[numberEnabledParams];
    for(Int_t ix{0}; ix < numberEnabledParams; ++ ix)
    {
        AdjustActs[ix] = 1.0;
        AdjustActs_Err[ix] = 0.5; // TODO: set using parameter_names.list
    }
    */
    // TODO: fix this
    //xi_31_init_value = 0.0;
    //xi_31_init_error = 0.0;
    //get_paramInitValueError(thePhase, 1, xi_31_init_value, xi_31_init_error);
    //AdjustActs[1] = xi_31_init_value;

    /*
    std::cout << "initial value is set to " << xi_31_init_value << std::endl;

    std::vector<double> init_par;
    std::vector<double> init_err;
    for(Int_t ix{0}; ix < numberEnabledParams; ++ ix)
    {
        init_par.push_back(AdjustActs[ix]);
        init_err.push_back(AdjustActs_Err[ix]);
    }
    */

    /*
    Int_t number_free_params = minuit->GetNumFreePars();
    Double_t *CovMatrix = new Double_t[number_free_params * number_free_params];
    for(int ix{0}; ix < number_free_params * number_free_params; ++ ix)
    {
        CovMatrix[ix] = 0.;
    }
    */
    

    #if 0
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
    #endif





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




    ///////////////////////////////////////////////////////////////////////////
    // HSD fixed xi_31 = HSD fit
    ///////////////////////////////////////////////////////////////////////////

    // do not do this in parallel mode
    if(0)// || (MODE_PARALLEL == 0))
    {
        // create minimizer
        ROOT::Minuit2::MnUserParameterState theParameterStateBefore;
        ROOT::Minuit2::VariableMetricMinimizer theMinimizer;
        MinimizeFCNAxialVector theFCN;

        // initialize fit
        //fitBackgrounds_init(theParameterState, theMinimizer, AdjustActs, AdjustActs_Err);
        //const double xi_31_value = xi_31_init_value;
        //const double xi_31_error = xi_31_init_error;
        const int xi_31_param_number = g_pg.get_xi_31_ext_param_number();
        const double xi_31_value = g_pg.file_params.at(xi_31_param_number).paramInitValue;
        const double xi_31_error = g_pg.file_params.at(xi_31_param_number).paramInitError;
        std::cout << "xi_31_param_number=" << xi_31_param_number
                  << " xi_31=" << xi_31_value << " +- " << xi_31_error << std::endl;
        fitBackgrounds_init(theParameterStateBefore, theMinimizer, xi_31_value, xi_31_error);


        // fix xi_31 parameter
        TString i_str;
        i_str.Form("%i", 1);
        TString minuit_param_number_str;
        minuit_param_number_str.Form("%i", 1);
        TString minuit_param_name = "_" + i_str + "_" + minuit_param_number_str + "_";
        theParameterStateBefore.Fix(std::string(minuit_param_name));
        theParameterStateBefore.SetValue(std::string(minuit_param_name), 0.0); // HSD


        // get parameters and chi2 value before fit
        std::vector<double> params_before = theParameterStateBefore.Params();
        std::vector<double> param_errs_before = theParameterStateBefore.Errors();

        for(int i = 0; i < params_before.size(); ++ i)
        {
            std::cout << "i=" << i << " param[i]=" << params_before.at(i) << std::endl;
        }
        
        double fval_before = theFCN.operator()(params_before);
        //int ndf = theFCN.ndf - theParameterStateBefore.VariableParameters();
        int ndf = theFCN.ndf - g_pg.get_number_free_params();

        // draw before fit
        draw_input_data drawinputdata;
        drawinputdata.chi2 = fval_before;
        drawinputdata.ndf = ndf;
        drawinputdata.serial_dir = "HSD";
        drawinputdata.saveas_filename = "HSD_before";
        drawinputdata.saveas_png = true;
       
        draw(drawinputdata,
             params_before,
             param_errs_before);

        // exec fit
        // this will fit backgrounds and the 150Nd amplitude parameter
        // but xi_31 is fixed
        ROOT::Minuit2::FunctionMinimum FCN_min =
            fitBackgrounds_exec(
                theParameterStateBefore,
                theMinimizer,
                theFCN);

        // get result
        ROOT::Minuit2::MnUserParameterState theParameterStateAfter = FCN_min.UserParameters();
        std::vector<double> params_after = theParameterStateAfter.Params();
        std::vector<double> param_errs_after = theParameterStateAfter.Errors();

        for(int i = 0; i < params_after.size(); ++ i)
        {
            std::cout << "i=" << i << " param[i]=" << params_after.at(i) << std::endl;
        }

        double fval_after = theFCN.operator()(params_after);
        //ndf = theFCN.ndf - theParameterStateAfter.VariableParameters();
        ndf = theFCN.ndf - g_pg.get_number_free_params();

        // draw after fit
        drawinputdata.chi2 = fval_after;
        drawinputdata.ndf = ndf;
        drawinputdata.saveas_filename = "HSD_after";
       
        draw(drawinputdata,
             params_after,
             param_errs_after);

        theParameterStateBefore.Release(std::string(minuit_param_name));
    }



#if 1
    ///////////////////////////////////////////////////////////////////////////
    // SSD fixed xi_31 = SSD fit
    ///////////////////////////////////////////////////////////////////////////

    // do fit for SSD
    // TODO: this block doesn't really make sense, unless we fit
    // for xi_31 = SSD with xi_31 fixed
    // do not do this in parallel mode
    if(0)// || (MODE_PARALLEL == 0))
    {
        // create minimizer
        ROOT::Minuit2::MnUserParameterState theParameterStateBefore;
        ROOT::Minuit2::VariableMetricMinimizer theMinimizer;
        MinimizeFCNAxialVector theFCN;

        // initialize fit
        //fitBackgrounds_init(theParameterState, theMinimizer, AdjustActs, AdjustActs_Err);
        const int xi_31_param_number = g_pg.get_xi_31_ext_param_number();
        const double xi_31_value = g_pg.file_params.at(xi_31_param_number).paramInitValue;
        const double xi_31_error = g_pg.file_params.at(xi_31_param_number).paramInitError;
        std::cout << "xi_31_param_number=" << xi_31_param_number
                  << " xi_31=" << xi_31_value << " +- " << xi_31_error << std::endl;
        fitBackgrounds_init(theParameterStateBefore, theMinimizer, xi_31_value, xi_31_error);

        // fix xi_31 parameter
        TString i_str;
        i_str.Form("%i", 1);
        TString minuit_param_number_str;
        minuit_param_number_str.Form("%i", 1);
        TString minuit_param_name = "_" + i_str + "_" + minuit_param_number_str + "_";
        theParameterStateBefore.Fix(std::string(minuit_param_name));
        theParameterStateBefore.SetValue(std::string(minuit_param_name), 0.296); // SSD

        // get parameters and chi2 value before fit
        std::vector<double> params_before = theParameterStateBefore.Params();
        std::vector<double> param_errs_before = theParameterStateBefore.Errors();
        double fval_before = theFCN.operator()(params_before);
        //int ndf = theFCN.ndf - theParameterStateBefore.VariableParameters();
        int ndf = theFCN.ndf - g_pg.get_number_free_params();

        // draw before fit
        draw_input_data drawinputdata;
        drawinputdata.chi2 = fval_before;
        drawinputdata.ndf = ndf;
        drawinputdata.serial_dir = "SSD";
        drawinputdata.saveas_filename = "SSD_before";
        drawinputdata.saveas_filename = "SSD_systematic1_before";
        drawinputdata.saveas_png = true;
       
        draw(drawinputdata,
             params_before,
             param_errs_before);
        
        // exec fit
        // this will fit backgrounds and the 150Nd amplitude parameter
        // but xi_31 is fixed
        ROOT::Minuit2::FunctionMinimum FCN_min =
            fitBackgrounds_exec(
                theParameterStateBefore,
                theMinimizer,
                theFCN);

        // get result
        ROOT::Minuit2::MnUserParameterState theParameterStateAfter = FCN_min.UserParameters();
        std::vector<double> params_after = theParameterStateAfter.Params();
        std::vector<double> param_errs_after = theParameterStateAfter.Errors();

        double fval_after = theFCN.operator()(params_after);
        //ndf = theFCN.ndf - theParameterStateAfter.VariableParameters();
        ndf = theFCN.ndf - g_pg.get_number_free_params();

        // draw result
        drawinputdata.chi2 = fval_after;
        drawinputdata.ndf = ndf;
        drawinputdata.saveas_filename = "SSD_after";
        drawinputdata.saveas_filename = "SSD_systematic1_after";
       
        draw(drawinputdata,
             params_after,
             param_errs_after);

        theParameterStateBefore.Release(std::string(minuit_param_name));
    }
#endif

    #if 0
    std::cout << "AdjustActs: ";
    for(int i = 0; i < numberEnabledParams; ++ i)
    {
        std::cout << AdjustActs[i] << " ";
    }
    std::cout << std::endl;
    #endif 

    /*
    std::cout << "params_before: ";
    for(int i = 0; i < params_before.size(); ++ i)
    {
        std::cout << params_before.at(i) << " ";
    }
    std::cout << "fval_before=" << fval_before << std::endl;
    std::cout << std::endl;
    cin_wait();
    */

/*
    ROOT::Minuit2::FunctionMinimum FCN_min =
        fitBackgrounds(
            theParameterState,
            theMinimizer,
            theFCN,
            //AdjustActs,
            //AdjustActs_Err,
            //CovMatrix,
            //number_free_params);
*/

#if 0
    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    ///////////////////////////////////////////////////////////////////////////

#if 1
    // do not do this in parallel mode
    if(1) // || (MODE_PARALLEL == 0))
    {
        // create minimizer
        ROOT::Minuit2::MnUserParameterState theParameterStateBefore;
        ROOT::Minuit2::VariableMetricMinimizer theMinimizer;
        MinimizeFCNAxialVector theFCN;

        // initialize fit
        //fitBackgrounds_init(theParameterState, theMinimizer, AdjustActs, AdjustActs_Err);
        const int xi_31_param_number = g_pg.get_xi_31_ext_param_number();
        const double xi_31_value = g_pg.file_params.at(xi_31_param_number).paramInitValue;
        const double xi_31_error = g_pg.file_params.at(xi_31_param_number).paramInitError;
        std::cout << "xi_31_param_number=" << xi_31_param_number
                  << " xi_31=" << xi_31_value << " +- " << xi_31_error << std::endl;
        fitBackgrounds_init(theParameterStateBefore, theMinimizer, xi_31_value, xi_31_error);

        // get parameters and chi2 value before fit
        std::vector<double> params_before = theParameterStateBefore.Params();
        std::vector<double> param_errs_before = theParameterStateBefore.Errors();
        double fval_before = theFCN.operator()(params_before);
        //int ndf = theFCN.ndf - theParameterStateBefore.VariableParameters();
        int ndf = theFCN.ndf - g_pg.get_number_free_params();

        // draw before fit
        draw_input_data drawinputdata;
        drawinputdata.chi2 = fval_before;
        drawinputdata.ndf = ndf;
        drawinputdata.serial_dir = "xifree";
        drawinputdata.saveas_filename = "xifree_before";
        drawinputdata.saveas_png = true;
       
        draw(drawinputdata,
             params_before,
             param_errs_before);

        // exec fit
        // do fit with all parameters free
        ROOT::Minuit2::FunctionMinimum FCN_min =
            fitBackgrounds_exec(
                theParameterStateBefore,
                theMinimizer,
                theFCN);

        // get result
        ROOT::Minuit2::MnUserParameterState theParameterStateAfter = FCN_min.UserParameters();
        std::vector<double> params_after = theParameterStateAfter.Params();
        std::vector<double> param_errs_after = theParameterStateAfter.Errors();

        double fval_after = theFCN.operator()(params_after);
        //ndf = theFCN.ndf - theParameterStateAfter.VariableParameters();
        ndf = theFCN.ndf - g_pg.get_number_free_params();

        // draw result
        drawinputdata.chi2 = fval_after;
        drawinputdata.ndf = ndf;
        drawinputdata.saveas_filename = "xifree_after";
       
        draw(drawinputdata,
             params_after,
             param_errs_after);

        
        // minimize
        //ROOT::Minuit2::FunctionMinimum FCN_min = theMinimizer.Minimize(theFCN, init_par, init_err);
        //ROOT::Minuit2::FunctionMinimum FCN_min = theMinimizer.Minimize(theFCN, init_par, init_err);
        /*
        std::cout << "Minimization finished" << std::endl;
        std::cout << "minimum: " << FCN_min << std::endl;
        std::cout << "chi2: " << FCN_min.Fval() << std::endl;
        std::cout << "edm: " << FCN_min.Edm() << std::endl;
        */


        std::cout << "fval_after=" << fval_after << " for params_after[0]=" << params_after[0] << " params_after[1]=" << params_after[1] << std::endl;
        std::cout << "fval_before=" << fval_before << std::endl;
    }
#endif
#endif



    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit - Variable Systematic Parameter
    ///////////////////////////////////////////////////////////////////////////

#if 1
    TCanvas *results_c = nullptr;
    TCanvas *results_c_A = nullptr;
    TCanvas *results_c_chi2_before = nullptr;
    TCanvas *results_c_chi2_after = nullptr;
    // do not do this in parallel mode
    if(1) // || (MODE_PARALLEL == 0))
    {

        std::vector<double> results_x;
        std::vector<double> results_y;
        std::vector<double> results_x_A;
        std::vector<double> results_y_A;
        std::vector<double> results_x_chi2_before;
        std::vector<double> results_y_chi2_before;
        std::vector<double> results_x_chi2_after;
        std::vector<double> results_y_chi2_after;

        const int i_max = 100;
        for(int i = 0; i <= i_max; ++ i)
        {

            const double min = -0.1;
            const double max = +0.1;
            const double diff = max - min;
            const double fr = (double)i / (double)i_max;
            gSystematics.systematic_energy_offset = fr * diff + min;
            double systematic_energy_offset = gSystematics.systematic_energy_offset;
            std::cout << "seo=" << systematic_energy_offset << std::endl;
            rebuild_fake_data_systematics(0.296, xi_31_baseline); // want to check if the fitter can fit itself to itsel

            std::string name_extra = "seo_" + std::to_string(systematic_energy_offset);

            // create minimizer
            ROOT::Minuit2::MnUserParameterState theParameterStateBefore;
            ROOT::Minuit2::VariableMetricMinimizer theMinimizer;
            MinimizeFCNAxialVector theFCN;

            // initialize fit
            //fitBackgrounds_init(theParameterState, theMinimizer, AdjustActs, AdjustActs_Err);
            const int xi_31_param_number = g_pg.get_xi_31_ext_param_number();
            const double xi_31_value = g_pg.file_params.at(xi_31_param_number).paramInitValue;
            const double xi_31_error = g_pg.file_params.at(xi_31_param_number).paramInitError;
            std::cout << "xi_31_param_number=" << xi_31_param_number
                      << " xi_31=" << xi_31_value << " +- " << xi_31_error << std::endl;
            fitBackgrounds_init(theParameterStateBefore, theMinimizer, xi_31_value, xi_31_error);

            // get parameters and chi2 value before fit
            std::vector<double> params_before = theParameterStateBefore.Params();
            std::vector<double> param_errs_before = theParameterStateBefore.Errors();
            double fval_before = theFCN.operator()(params_before);
            //int ndf = theFCN.ndf - theParameterStateBefore.VariableParameters();
            int ndf = theFCN.ndf - g_pg.get_number_free_params();

            // draw before fit
            draw_input_data drawinputdata;
            drawinputdata.chi2 = fval_before;
            drawinputdata.ndf = ndf;
            drawinputdata.serial_dir = "xifree";
            drawinputdata.saveas_filename = std::string("xifree_before") + "_" + name_extra;
            drawinputdata.saveas_png = true;
           
            draw(drawinputdata,
                 params_before,
                 param_errs_before);


            // exec fit
            // do fit with all parameters free
            ROOT::Minuit2::FunctionMinimum FCN_min =
                fitBackgrounds_exec(
                    theParameterStateBefore,
                    theMinimizer,
                    theFCN);

            // get result
            ROOT::Minuit2::MnUserParameterState theParameterStateAfter = FCN_min.UserParameters();
            std::vector<double> params_after = theParameterStateAfter.Params();
            std::vector<double> param_errs_after = theParameterStateAfter.Errors();

            double fval_after = theFCN.operator()(params_after);
            //ndf = theFCN.ndf - theParameterStateAfter.VariableParameters();
            ndf = theFCN.ndf - g_pg.get_number_free_params();

            // draw result
            drawinputdata.chi2 = fval_after;
            drawinputdata.ndf = ndf;
            drawinputdata.saveas_filename = std::string("xifree_after") + "_" + name_extra;
           
            draw(drawinputdata,
                 params_after,
                 param_errs_after);

            
            // minimize
            //ROOT::Minuit2::FunctionMinimum FCN_min = theMinimizer.Minimize(theFCN, init_par, init_err);
            //ROOT::Minuit2::FunctionMinimum FCN_min = theMinimizer.Minimize(theFCN, init_par, init_err);
            /*
            std::cout << "Minimization finished" << std::endl;
            std::cout << "minimum: " << FCN_min << std::endl;
            std::cout << "chi2: " << FCN_min.Fval() << std::endl;
            std::cout << "edm: " << FCN_min.Edm() << std::endl;
            */


            std::cout << "fval_after=" << fval_after << " for params_after[0]=" << params_after[0] << " params_after[1]=" << params_after[1] << std::endl;
            std::cout << "fval_before=" << fval_before << std::endl;

            results_x.push_back(systematic_energy_offset);
            results_y.push_back(params_after.at(1));

            results_x_A.push_back(systematic_energy_offset);
            results_y_A.push_back(params_after.at(0));
            
            results_x_chi2_before.push_back(systematic_energy_offset);
            results_y_chi2_before.push_back(fval_before);

            results_x_chi2_after.push_back(systematic_energy_offset);
            results_y_chi2_after.push_back(fval_after);
        }

        TGraph *results_g = new TGraph(results_x.size(), results_x.data(), results_y.data());
        TGraph *results_g_A = new TGraph(results_x_A.size(), results_x_A.data(), results_y_A.data());
        TGraph *results_g_chi2_before = new TGraph(results_x_chi2_before.size(), results_x_chi2_before.data(), results_y_chi2_before.data());
        TGraph *results_g_chi2_after = new TGraph(results_x_chi2_after.size(), results_x_chi2_after.data(), results_y_chi2_after.data());

        results_c = new TCanvas("results", "results");
        results_g->Draw();

        results_c_A = new TCanvas("results_A", "results_A");
        results_g_A->Draw();

        results_c_chi2_before = new TCanvas("results_chi2_before", "results_chi2_before");
        results_g_chi2_before->Draw();

        results_c_chi2_after = new TCanvas("results_chi2_after", "results_chi2_after");
        results_g_chi2_after->Draw();
    }
#endif

//    draw_channel(1, params, fval, "test_channel_1.png");
//    TH1D *j1, *j2, *j3, *j4;
// TODO: re-enable
//    draw(params, param_errs, fval, j1, j2, j3, j4, "minuit_output_fake_data.*", ".", false);
    //draw(params, nullptr, "NOSAVE", fval, j1, j2, j3, j4, true);

//    std::cin.get();

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

    // reenable this one
    if(0)
    {
        gSystematics.systematic_energy_offset = 0.0;
        double systematic_energy_offset = gSystematics.systematic_energy_offset;
        std::cout << "seo=" << systematic_energy_offset << std::endl;
        rebuild_fake_data_systematics(0.296, xi_31_baseline); // want to check if the fitter can fit itself to itsel

        // create minimizer
        //ROOT::Minuit2::MnUserParameterState theParameterState;
        //ROOT::Minuit2::VariableMetricMinimizer theMinimizer;
        //MinimizeFCNAxialVector theFCN;

        // initialize fit
        //fitBackgrounds_init(theParameterState, theMinimizer, AdjustActs, AdjustActs_Err);
        //const double xi_31_value = xi_31_init_value;
        //const double xi_31_error = xi_31_init_error;
        // TODO: this function does not work as expected with the testmps function below
        //fitBackgrounds_phasespace_init(theParameterState, theMinimizer, 1.0, xi_31_value);

        std::cout << "ready to test MPS" << std::endl;
        std::cout << "output_name=" << output_name << std::endl;
        std::cout << "start_index=" << start_index << std::endl;
        std::cout << "stop_index=" << stop_index << std::endl;
        //std::cin.get();
        //newloglikfitter_testmyphasespace(theParameterState, theFCN); //, FCN_min);
        newloglikfitter_testmyphasespace_newversion(number_job_id, output_name, start_index, stop_index);
        std::cout << "MPS done" << std::endl;
    }
   
/*
    if(0)
    {
        //newloglikfitter_gA_chisquaretest(minuit, AdjustActs, AdjustActs_Err);
    }
*/
    ///////////////////////////////////////////////////////////////////////////
#if 0
    print_adjustacts(std::cout, params, param_errs);
    std::ofstream myFileFitResults("fit_results.txt", std::ios::out | std::ios::app);
    print_adjustacts(myFileFitResults, params, param_errs);
    myFileFitResults.close();


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
#endif
    
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

#if 0
    std::ofstream of_numberofeventsafterfit("of_numberofeventsafterfit.txt", std::ofstream::out | std::ofstream::app);
    timestamp(of_numberofeventsafterfit);
    for(int i = 0; i < allMCSamples1D[0]->GetEntries(); ++ i)
    {
        TH1D *tmpHist = (TH1D*)allMCSamples1D[0]->At(i);
        Double_t integral = tmpHist->Integral();
        of_numberofeventsafterfit << tmpHist->GetName() << " number of events " << integral << std::endl;
    }
    of_numberofeventsafterfit.close();
#endif

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
















