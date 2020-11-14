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
#include "newLogLikFitter_logpoisson.h"
#include "newLogLikFitter_aux.h"
#include "newLogLikFitter_printfitresult.h"
#include "newLogLikFitter_read_parameternames_lst.h"
#include "Systematics.h"
#include "newLogLikFitter_reweight.h"
#include "newLogLikFitter_reweight_apply_MC.h"
#include "newLogLikFitter_reweight_apply_fakedata.h"
#include "newLogLikFitter_reweight_apply_data.h"
#include "newLogLikFitter_buildfakedata.h"
#include "newLogLikFitter_book1DHistograms.h"
#include "newLogLikFitter_book2DHistograms.h"
#include "newLogLikFitter_stacker_helper.h"
#include "newLogLikFitter_drawchannel.h"
#include "newLogLikFitter_draw_channel_phase_with_pull.h"
#include "newLogLikFitter_draw.h"
#include "newLogLikFitter_draw_2D.h"
#include "newLogLikFitter_draw_outputdiff.h"
#include "newLogLikFitter_draw_all.h"
#include "newLogLikFitter_rebuild_150Nd_MC.h"
#include "newLogLikFitter_rebuild_150Nd_data.h"
#include "newLogLikFitter_rebuild_150Nd_fakedata.h"
#include "MinimizeFCNAxialVector.h"
#include "newLogLikFitter_fitBackgrounds.h"
#include "newLogLikFitter_mps_aux.h"
#include "newLogLikFitter_mps_draw.h"
#include "newLogLikFitter_mps_calculate.h"
//#include "newLogLikFitter_chisquaretest.h"
#include "newLogLikFitter_test.h"
//#include "newLogLikFitter_preMPSfitdriver.h"
#include "newLogLikFitter_min_point_helper.h"




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




void final_mps_draw()
{
    std::cout << "calling draw" << std::endl;
    
    int number_job_id = 0;
    std::string output_name = "noparallel";
    int start_index = 0;
    int stop_index = 301;
    if(MODE_PARALLEL == 1)
    {
        //bool success = load_from_script(i, number_job_id, output_name, start_index, stop_index);
        //if(success == true)
        //{
        //    std::cout << "Job Init: ID=" << number_job_id << std::endl;
        //    std::cout << "output_name=" << output_name << std::endl;
        //    std::cout << "START_INDEX=" << start_index << std::endl;
        //    std::cout << "STOP_INDEX=" << stop_index << std::endl;
        //    // do nothing = continue
        //}
        //else
        //{
        //    std::cout << "fail" << std::endl;
        //    return;
        //}
        output_name = "mps_resultsmatrix";
    }
    else
    {
        number_job_id = 0;
        output_name = "noparallel";
        start_index = 0;
        stop_index = 301;
    }

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    //gStyle->SetPalette(kBird);
    //gStyle->SetPalette(kBrownCyan);
    gStyle->SetPalette(kLightTemperature);
    //gStyle->SetNumberContours(1000);

    newloglikfitter_mps_draw_systematics(number_job_id, output_name, start_index, stop_index);
}



void final_mps_draw_colz()
{
    std::cout << "calling draw (colz)" << std::endl;
    
    int number_job_id = 0;
    std::string output_name = "noparallel";
    int start_index = 0;
    int stop_index = 301;
    if(MODE_PARALLEL == 1)
    {
        output_name = "mps_resultsmatrix";
    }
    else
    {
        number_job_id = 0;
        output_name = "noparallel";
        start_index = 0;
        stop_index = 301;
    }

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    //gStyle->SetPalette(kBird);
    //gStyle->SetPalette(kBrownCyan);
    gStyle->SetPalette(kLightTemperature);
    //gStyle->SetNumberContours(1000);

    newloglikfitter_mps_draw_systematics_colz(number_job_id, output_name, start_index, stop_index);
}



bool load_from_script(
    int i,
    int &number_job_id,
    std::string &output_name,
    int &start_index,
    int &stop_index
    )
{
    int script_index = 1; // NOTE: zero reserved for concatinated JID files
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
                            ifs.close();
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







void newLogLikFitter()
{
    loadFiles(0);
}

void newLogLikFitter(int i)
{
    loadFiles(i);
}










void systematic_n_init_helper
(
    std::vector<double> *systematic_n_highlow_1D_P1[][number1DHists],
    std::vector<double> *systematic_n_highlow_1D_P2[][number1DHists],
    const int systematic_n_index
)
{

    // loop over all channels
    for(int channel = 0; channel < number1DHists; ++ channel)
    {

        std::string histname = std::string(channel_histname_1D[channel]);
        std::string search_object_P1;
        std::string search_object_P2;

        search_object_P1 = histname + std::string("fakedata") + "_P1";
        search_object_P2 = histname + std::string("fakedata") + "_P2";

        TH1D *tmpDataHist1D_P1 = nullptr;
        TH1D *tmpDataHist1D_P2 = nullptr;
    
        tmpDataHist1D_P1 = (TH1D*)allFakeDataSamples1D->FindObject(search_object_P1.c_str());
        tmpDataHist1D_P2 = (TH1D*)allFakeDataSamples1D->FindObject(search_object_P2.c_str());

        // phase 1
        if(tmpDataHist1D_P1 == nullptr)
        {
            std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
            throw "problem";
        }

        // phase 2
        if(tmpDataHist1D_P2 == nullptr)
        {
            std::cout << "ERROR: Could not find object " << search_object_P2 << std::endl;
            throw "problem";
        }

        // phase 1
        for(Int_t bin_ix{1}; bin_ix <= tmpDataHist1D_P1->GetNbinsX(); ++ bin_ix)
        {
            Double_t content = tmpDataHist1D_P1->GetBinContent(bin_ix);
            //systematic_offset_high_1D_P1[channel]->push_back(content);
            systematic_n_highlow_1D_P1[systematic_n_index][channel]->push_back(content);
        }

        // phase 2
        for(Int_t bin_ix{1}; bin_ix <= tmpDataHist1D_P2->GetNbinsX(); ++ bin_ix)
        {
            Double_t content = tmpDataHist1D_P2->GetBinContent(bin_ix);
            //systematic_offset_high_1D_P2[channel]->push_back(content);
            systematic_n_highlow_1D_P2[systematic_n_index][channel]->push_back(content);
        }

    }
}











void loadFiles(int i)
{
    
    std::chrono::system_clock::time_point start_time = std::chrono::high_resolution_clock::now();
        


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
        stop_index = 301;
    }


    ///////////////////////////////////////////////////////////////////////////
    // read parameter file list
    ///////////////////////////////////////////////////////////////////////////

    // read parameter_name.lst file
    // read parameter list file
    read_parameter_list_file();
    //g_pg.print();


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

    xi_31_baseline = 0.0;

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



    std::cout << "Loading histograms" << std::endl;
    // 1d: Phase 1 & 2
    for(int channel = 0; channel < number1DHists; ++ channel)
    {
        //book1DHistograms(0, "2e_", "hTotalE_");
        book1DHistograms(channel, "2e_", channel_histname_1D[channel]);
    }

    // 2d: Phase 1 & 2
    // TODO?
    for(int channel = 0; channel < number2DHists; ++ channel)
    {
        book2DHistograms(channel, "2e_", channel_histname_2D[channel]);
    }
    std::cout << "All histograms loaded" << std::endl;





    //gSystematics.systematic_energy_offset = 0.0; //-0.1;
    //gSystematics.systematic_energy_scale = 0.0; //-0.1;
    gSystematics.reset();
    /*const*/ int xi_31_ext_param_number = g_pg.get_xi_31_ext_param_number();
        //if(param[xi_31_ext_param_number] != g_pg.file_params.at(xi_31_ext_param_number).paramLastValue)
    //const double xi_31{param[xi_31_ext_param_number]};
    double xi_31 = 0.0;
            if(g_pg.get_xi_31_int_param_number() != -1)
            {
                xi_31_ext_param_number = g_pg.get_xi_31_ext_param_number();
            }
            else
            {
    /*const double*/ xi_31 = g_pg.file_params.at(xi_31_ext_param_number).paramInitValue;
            }
    std::cout << "xi_31=" << xi_31 << std::endl;
//    std::cin.get();
    // NOTE: you put whatever value of xi_31 you like in here
    // to construct fakedata with that xi_31 value
    const double xi_31_SSD = 0.296;
    //rebuild_fake_data_systematics(xi_31, xi_31_baseline);
    //rebuild_fake_data_systematics(0.0, xi_31_baseline); // want to check if the fitter can fit itself to itsel
    //rebuild_fake_data_systematics(xi_31_SSD, xi_31_baseline); // want to check if the fitter can fit itself to itsel


    const double xi_31_systematics_reweight_value = 0.7; //xi_31_SSD;
    //const double xi_31_covariance_matrix_reweight_value = xi_31_SSD;
    //const double xi_31_covariance_matrix_reweight_value = 0.0;
    const double xi_31_covariance_matrix_reweight_value = 0.7;








    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // All Systematics OFF
    ///////////////////////////////////////////////////////////////////////////
    
    //gSystematics.systematic_energy_offset = 0.0;
    //gSystematics.systematic_energy_scale = 0.0;
    gSystematics.reset();
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);

    // loop over all channels
    for(int channel = 0; channel < number1DHists; ++ channel)
    {

        if(systematic_nominal_1D_P1[channel] == nullptr)
        {
            systematic_nominal_1D_P1[channel] = new std::vector<double>;
            systematic_nominal_1D_P2[channel] = new std::vector<double>;

            for(int i = 0; i < N_SYSTEMATICS; ++ i)
            {
                // P1
                systematic_n_low_1D_P1[i][channel] = new std::vector<double>;
                systematic_n_high_1D_P1[i][channel] = new std::vector<double>;
                systematic_n_V_MATRIX_coeff_1D_P1[i][channel] = new std::vector<double>;

                // P2
                systematic_n_low_1D_P2[i][channel] = new std::vector<double>;
                systematic_n_high_1D_P2[i][channel] = new std::vector<double>;
                systematic_n_V_MATRIX_coeff_1D_P2[i][channel] = new std::vector<double>;
            }

        }
        // check channel enabled
        //if(channel_enable_1D[channel] == 0)
        //{
        //    continue;
        //}

        std::string histname = std::string(channel_histname_1D[channel]);
        std::string search_object_P1;
        std::string search_object_P2;

        search_object_P1 = histname + std::string("fakedata") + "_P1";
        search_object_P2 = histname + std::string("fakedata") + "_P2";

        TH1D *tmpDataHist1D_P1 = nullptr;
        TH1D *tmpDataHist1D_P2 = nullptr;
    
        tmpDataHist1D_P1 = (TH1D*)allFakeDataSamples1D->FindObject(search_object_P1.c_str());
        tmpDataHist1D_P2 = (TH1D*)allFakeDataSamples1D->FindObject(search_object_P2.c_str());

        if(tmpDataHist1D_P1 == nullptr)
        {
            std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
            throw "problem";
        }
        if(tmpDataHist1D_P2 == nullptr)
        {
            std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
            throw "problem";
        }
        for(Int_t bin_ix{1}; bin_ix <= tmpDataHist1D_P1->GetNbinsX(); ++ bin_ix)
        {
            Double_t content = tmpDataHist1D_P1->GetBinContent(bin_ix);
            systematic_nominal_1D_P1[channel]->push_back(content);
        }
        for(Int_t bin_ix{1}; bin_ix <= tmpDataHist1D_P2->GetNbinsX(); ++ bin_ix)
        {
            Double_t content = tmpDataHist1D_P2->GetBinContent(bin_ix);
            systematic_nominal_1D_P2[channel]->push_back(content);
        }

    }
    gSystematics.reset();


    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Systematic: Linear Energy Shift
    // Systematic Value: +0.1 MeV
    ///////////////////////////////////////////////////////////////////////////

    gSystematics.reset();
    gSystematics.systematic_energy_offset = +0.1;
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);
    systematic_n_init_helper(systematic_n_high_1D_P1, systematic_n_high_1D_P2, 0);
    gSystematics.reset();


    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Systematic: Linear Energy Shift
    // Systematic Value: -0.1 MeV
    ///////////////////////////////////////////////////////////////////////////

    gSystematics.reset();
    gSystematics.systematic_energy_offset = -0.1;
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);
    systematic_n_init_helper(systematic_n_low_1D_P1, systematic_n_low_1D_P2, 0);
    gSystematics.reset();


    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Systematic: Energy Scale Multiplier
    // Systematic Value: -0.012 -(1.0 % + 0.2 %)
    ///////////////////////////////////////////////////////////////////////////

    gSystematics.reset();
    gSystematics.systematic_energy_scale = -0.012;
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);
    systematic_n_init_helper(systematic_n_low_1D_P1, systematic_n_low_1D_P2, 1);
    gSystematics.reset();


    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Systematic: Linear Energy Shift
    // Systematic Value: +0.012 (1.0% + 0.2%)
    ///////////////////////////////////////////////////////////////////////////

    gSystematics.reset();
    gSystematics.systematic_energy_scale = +0.012;
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);
    systematic_n_init_helper(systematic_n_high_1D_P1, systematic_n_high_1D_P2, 1);
    gSystematics.reset();



    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Systematic: Efficiency
    // Systematic Value: -5.55 %
    ///////////////////////////////////////////////////////////////////////////

    gSystematics.reset();
    gSystematics.systematic_efficiency = -5.55e-02;
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);
    systematic_n_init_helper(systematic_n_low_1D_P1, systematic_n_low_1D_P2, 2);
    gSystematics.reset();


    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Systematic: Efficiency
    // Systematic Value: +5.55 %
    ///////////////////////////////////////////////////////////////////////////

    gSystematics.reset();
    gSystematics.systematic_efficiency = +5.55e-02;
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);
    systematic_n_init_helper(systematic_n_high_1D_P1, systematic_n_high_1D_P2, 2);
    gSystematics.reset();


    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Systematic: Enrichment
    // Systematic Value: -0.5 %
    ///////////////////////////////////////////////////////////////////////////

    gSystematics.reset();
    gSystematics.systematic_enrichment = -0.5e-02;
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);
    systematic_n_init_helper(systematic_n_low_1D_P1, systematic_n_low_1D_P2, 3);
    gSystematics.reset();


    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Systematic: Enrichment
    // Systematic Value: +0.5 %
    ///////////////////////////////////////////////////////////////////////////

    gSystematics.reset();
    gSystematics.systematic_enrichment = +0.5e-02;
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);
    systematic_n_init_helper(systematic_n_high_1D_P1, systematic_n_high_1D_P2, 3);
    gSystematics.reset();


    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Systematic: Linear Energy Shift
    // Systematic Value: +3.0e-03 MeV
    ///////////////////////////////////////////////////////////////////////////

    gSystematics.reset();
    gSystematics.systematic_energy_offsetsmall = +3.0e-3;
    //gSystematics.systematic_energy_scale = 0.0;
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);
    systematic_n_init_helper(systematic_n_high_1D_P1, systematic_n_high_1D_P2, 4);
    gSystematics.reset();


    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Systematic: Linear Energy Shift
    // Systematic Value: -3.0e-03 MeV
    ///////////////////////////////////////////////////////////////////////////

    gSystematics.reset();
    gSystematics.systematic_energy_offsetsmall = -3.0e-3;
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);
    systematic_n_init_helper(systematic_n_low_1D_P1, systematic_n_low_1D_P2, 4);
    gSystematics.reset();


    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Systematic: Foil Thickness
    // Systematic Value: +
    ///////////////////////////////////////////////////////////////////////////

    gSystematics.reset();
    gSystematics.systematic_foil_thickness_virtual = +1;
    //gSystematics.systematic_energy_scale = 0.0;
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);
    systematic_n_init_helper(systematic_n_high_1D_P1, systematic_n_high_1D_P2, 5);
    gSystematics.reset();


    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Systematic: Foil Thickness
    // Systematic Value: -
    ///////////////////////////////////////////////////////////////////////////

    gSystematics.reset();
    gSystematics.systematic_foil_thickness_virtual = -1;
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);
    systematic_n_init_helper(systematic_n_low_1D_P1, systematic_n_low_1D_P2, 5);
    gSystematics.reset();


    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Systematic: dE/dX
    // Systematic Value: +
    ///////////////////////////////////////////////////////////////////////////

    gSystematics.reset();
    gSystematics.systematic_dEdX_virtual = +1;
    //gSystematics.systematic_energy_scale = 0.0;
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);
    systematic_n_init_helper(systematic_n_high_1D_P1, systematic_n_high_1D_P2, 6);
    gSystematics.reset();


    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Systematic: dE/dX
    // Systematic Value: -
    ///////////////////////////////////////////////////////////////////////////

    gSystematics.reset();
    gSystematics.systematic_dEdX_virtual = -1;
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);
    systematic_n_init_helper(systematic_n_low_1D_P1, systematic_n_low_1D_P2, 6);
    gSystematics.reset();


    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Systematic: brem
    // Systematic Value: +
    ///////////////////////////////////////////////////////////////////////////

    gSystematics.reset();
    gSystematics.systematic_brem_virtual = +1;
    //gSystematics.systematic_energy_scale = 0.0;
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);
    systematic_n_init_helper(systematic_n_high_1D_P1, systematic_n_high_1D_P2, 7);
    gSystematics.reset();


    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Systematic: brem
    // Systematic Value: -
    ///////////////////////////////////////////////////////////////////////////

    gSystematics.reset();
    gSystematics.systematic_brem_virtual = -1;
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);
    systematic_n_init_helper(systematic_n_low_1D_P1, systematic_n_low_1D_P2, 7);
    gSystematics.reset();


    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Systematic: Foil Thickness
    // Systematic Value: +
    ///////////////////////////////////////////////////////////////////////////

    gSystematics.reset();
    gSystematics.systematic_foil_thickness_nominal = +1;
    //gSystematics.systematic_energy_scale = 0.0;
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);
    systematic_n_init_helper(systematic_n_high_1D_P1, systematic_n_high_1D_P2, 8);
    gSystematics.reset();


    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Systematic: Foil Thickness
    // Systematic Value: -
    ///////////////////////////////////////////////////////////////////////////

    gSystematics.reset();
    gSystematics.systematic_foil_thickness_nominal = -1;
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);
    systematic_n_init_helper(systematic_n_low_1D_P1, systematic_n_low_1D_P2, 8);
    gSystematics.reset();


    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Systematic: dE/dX
    // Systematic Value: +
    ///////////////////////////////////////////////////////////////////////////

    gSystematics.reset();
    gSystematics.systematic_dEdX_nominal = +1;
    //gSystematics.systematic_energy_scale = 0.0;
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);
    systematic_n_init_helper(systematic_n_high_1D_P1, systematic_n_high_1D_P2, 9);
    gSystematics.reset();


    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Systematic: dE/dX
    // Systematic Value: -
    ///////////////////////////////////////////////////////////////////////////

    gSystematics.reset();
    gSystematics.systematic_dEdX_nominal = -1;
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);
    systematic_n_init_helper(systematic_n_low_1D_P1, systematic_n_low_1D_P2, 9);
    gSystematics.reset();


    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Systematic: brem
    // Systematic Value: +
    ///////////////////////////////////////////////////////////////////////////

    gSystematics.reset();
    gSystematics.systematic_brem_nominal = +1;
    //gSystematics.systematic_energy_scale = 0.0;
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);
    systematic_n_init_helper(systematic_n_high_1D_P1, systematic_n_high_1D_P2, 10);
    gSystematics.reset();


    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Systematic: brem
    // Systematic Value: -
    ///////////////////////////////////////////////////////////////////////////

    gSystematics.reset();
    gSystematics.systematic_brem_nominal = -1;
    rebuild_fake_data_systematics(xi_31_covariance_matrix_reweight_value, xi_31_baseline);
    systematic_n_init_helper(systematic_n_low_1D_P1, systematic_n_low_1D_P2, 10);
    gSystematics.reset();


    // reset systematics
    gSystematics.reset();
    //rebuild_fake_data_systematics(xi_31_SSD, xi_31_baseline);
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;
    rebuild_fake_data_systematics(xi_31_systematics_reweight_value, xi_31_baseline);


    ///////////////////////////////////////////////////////////////////////////
    // Construct Systematic Data Objects
    // Calculate Covariance Matrix Parameters
    ///////////////////////////////////////////////////////////////////////////

    for(int channel = 0; channel < number1DHists; ++ channel)
    {

        ///////////////////////////////////////////////////////////////////////
        // SYSTEMATIC N
        ///////////////////////////////////////////////////////////////////////
        
        for(int n = 0; n < N_SYSTEMATICS; ++ n)
        {

            // PHASE 1
            for(std::size_t i = 0; i < systematic_nominal_1D_P1[channel]->size(); ++ i)
            {
                double up = systematic_n_high_1D_P1[n][channel]->at(i);
                double nominal = systematic_nominal_1D_P1[channel]->at(i);
                double down = systematic_n_low_1D_P1[n][channel]->at(i);
                double value_up = up - nominal;
                double value_down = nominal - down;

                double value = 0.0;
                /*
                if(std::abs(value_up) >= std::abs(value_down))
                {
                    value = value_up;
                }
                else
                {
                    value = value_down;
                }*/

                /*
                double mp_fd_x = min_point_fake_data[0];
                double mp_fd_y = min_point_fake_data[1];
                double mp_sysn_l_x = min_point_sysn_l[n][0];
                double mp_sysn_l_y = min_point_sysn_l[n][1];
                double mp_sysn_h_x = min_point_sysn_h[n][0];
                double mp_sysn_h_y = min_point_sysn_h[n][1];

                double c_l_x = (mp_sysn_l_x - mp_fd_x);
                */

                // just do it manually
                // TODO: choose which values go here
                /*
                if((n == 0) ||
                   (n == 1) || 
                   (n == 2) || 
                   (n == 3) || 
                   (n == 4) ||
                   (n == 5) ||
                   (n == 6) ||
                   (n == 7) ||
                   (n == 8) ||
                   (n == 9) ||
                   (n == 10))
                {
                    value = value_up;
                }
                else
                {
                    value = value_down;
                }*/
                value = value_up;

                //systematic_offset_V_MATRIX_coeff_1D_P1[channel]->operator[](i) = value;
                systematic_n_V_MATRIX_coeff_1D_P1[n][channel]->push_back(value);
            }

            // PHASE 2
            for(std::size_t i = 0; i < systematic_nominal_1D_P2[channel]->size(); ++ i)
            {
                double up = systematic_n_high_1D_P2[n][channel]->at(i);
                double nominal = systematic_nominal_1D_P2[channel]->at(i);
                double down = systematic_n_low_1D_P2[n][channel]->at(i);
                double value_up = up - nominal;
                double value_down = nominal - down;

                double value = 0.0;
                /*
                if(std::abs(value_up) >= std::abs(value_down))
                {
                    value = value_up;
                }
                else
                {
                    value = value_down;
                }*/

                // just do it manually
                // TODO: choose which values go here
                /*
                if((n == 0) ||
                   (n == 1) || 
                   (n == 2) || 
                   (n == 3) || 
                   (n == 4) ||
                   (n == 5) ||
                   (n == 6) ||
                   (n == 7) ||
                   (n == 8) ||
                   (n == 9) ||
                   (n == 10))
                {
                    value = value_up;
                }
                else
                {
                    value = value_down;
                }*/
                value = value_up;

                //systematic_offset_V_MATRIX_coeff_1D_P2[channel]->operator[](i) = value;
                systematic_n_V_MATRIX_coeff_1D_P2[n][channel]->push_back(value);
            }

        }

    }





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
    







#if 0
    ///////////////////////////////////////////////////////////////////////////
    // Reproduce Summers Fit (only when xi_31 parameter is disabled)
    ///////////////////////////////////////////////////////////////////////////

    // Difference between this and HSD fit:
    // HSD fit is Summers Fit
    // This section scans the phase space of 150Nd amplitude parameter
    // (used to check the Minuit2 error reporting)

    // do not do this in parallel mode
    if(0)// || (MODE_PARALLEL == 0))
    {
        bool restore_V_ENABLE_SYSALL = V_ENABLE_SYSALL;
        V_ENABLE_SYSALL = false;

        std::cout << "Reproduction of Summers 150 Nd fit" << std::endl;
        std::cout << "Note that this only works if the xi_31 parameter is disabled in the parameter_names.lst file" << std::endl;

        TH1D *mps1D = new TH1D("mps1D", "mps1D", 50, 0.84, 1.02);
        TH1D *mps1D_before = new TH1D("mps1D", "mps1D", 50, 0.84, 1.02);
        double min_found = std::numeric_limits<double>::infinity();
        for(Int_t index = 0; index < mps1D->GetNbinsX(); ++ index)
        {

            Double_t amplitude = mps1D->GetBinCenter(index);

            // create minimizer
            ROOT::Minuit2::MnUserParameterState theParameterStateBefore;
            ROOT::Minuit2::VariableMetricMinimizer theMinimizer;
            MinimizeFCNAxialVector theFCN;

            // initialize fit
            //fitBackgrounds_init(theParameterState, theMinimizer, AdjustActs, AdjustActs_Err);
            //const double xi_31_value = xi_31_init_value;
            //const double xi_31_error = xi_31_init_error;
            const int xi_31_param_number = g_pg.get_xi_31_ext_param_number(); // don't know what this will do when this param is disabled
            const double xi_31_value = g_pg.file_params.at(xi_31_param_number).paramInitValue;
            const double xi_31_error = g_pg.file_params.at(xi_31_param_number).paramInitError;
            std::cout << "xi_31_param_number=" << xi_31_param_number
                      << " xi_31=" << xi_31_value << " +- " << xi_31_error << std::endl;
            fitBackgrounds_init(theParameterStateBefore, theMinimizer, xi_31_value, xi_31_error);

            TString i_str;
            i_str.Form("%i", 0);
            TString minuit_param_number_str;
            minuit_param_number_str.Form("%i", 0);
            TString minuit_param_name = "_" + i_str + "_" + minuit_param_number_str + "_";
            theParameterStateBefore.Fix(std::string(minuit_param_name));
            theParameterStateBefore.SetValue(std::string(minuit_param_name), amplitude);

            // fix xi_31 parameter
            // NOTE: do not do this, should be zero
            /*
            TString i_str;
            i_str.Form("%i", 1);
            TString minuit_param_number_str;
            minuit_param_number_str.Form("%i", 1);
            TString minuit_param_name = "_" + i_str + "_" + minuit_param_number_str + "_";
            theParameterStateBefore.Fix(std::string(minuit_param_name));
            theParameterStateBefore.SetValue(std::string(minuit_param_name), 0.0); // HSD
            */

            // get parameters and chi2 value before fit
            std::vector<double> params_before = theParameterStateBefore.Params();
            std::vector<double> param_errs_before = theParameterStateBefore.Errors();

            for(int i = 0; i < params_before.size(); ++ i)
            {
                std::cout << "i=" << i << " param[i]=" << params_before.at(i) << " +- " << param_errs_before.at(i) << std::endl;
            }
            
            double fval_before = theFCN.operator()(params_before);
            //int ndf = theFCN.ndf - theParameterStateBefore.VariableParameters();
            int nch = theFCN.nch;
            int nfp = g_pg.get_number_free_params();
            int ndf = nch - nfp;

            mps1D_before->SetBinContent(index, fval_before);

            // draw before fit
            draw_input_data drawinputdata;
            drawinputdata.chi2 = fval_before;
            drawinputdata.nch = nch;
            drawinputdata.nfp = nfp;
            drawinputdata.serial_dir = "ND150fit";
            drawinputdata.saveas_filename = "ND150_before";
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
                std::cout << "i=" << i << " param[i]=" << params_after.at(i) << " +- " << param_errs_after.at(i) << std::endl;
            }

            double fval_after = theFCN.operator()(params_after);
            //ndf = theFCN.ndf - theParameterStateAfter.VariableParameters();
            nch = theFCN.nch;
            nfp = g_pg.get_number_free_params();
            ndf = nch - nfp;

            mps1D->SetBinContent(index, fval_after);

            // draw after fit
            drawinputdata.chi2 = fval_after;
            drawinputdata.nch = nch;
            drawinputdata.nfp = nfp;
            drawinputdata.saveas_filename = "ND150_after";
           
            draw(drawinputdata,
                 params_after,
                 param_errs_after);

            if(fval_before < min_found)
            {
                min_found = fval_before;
            }
            if(fval_after < min_found)
            {
                min_found = fval_after;
            }

            //theParameterStateBefore.Release(std::string(minuit_param_name));
        
        }

        TCanvas *c_mps1D = new TCanvas("c_mps1D", "c_mps1D");
        //TCanvas *c_mps1D_before = new TCanvas("c_mps1D_before", "c_mps1D_before");
        mps1D_before->SetLineColor(kBlue);
        mps1D_before->GetYaxis()->SetRangeUser(min_found - 0.2, min_found + 1.5);
        mps1D_before->Draw("hist");
        mps1D->SetLineColor(kRed);
        mps1D->Draw("histsame");
        c_mps1D->SaveAs("c_mps1D.png");
        c_mps1D->SaveAs("c_mps1D.pdf");


        V_ENABLE_SYSALL = restore_V_ENABLE_SYSALL;

        // rest of analysis does not make sense
        return 0;
    }
#endif


#if 0
    ///////////////////////////////////////////////////////////////////////////
    // HSD fixed xi_31 = HSD fit
    ///////////////////////////////////////////////////////////////////////////

    // do not do this in parallel mode
    if(1)// || (MODE_PARALLEL == 0))
    {
        bool restore_V_ENABLE_SYSALL = V_ENABLE_SYSALL;
        V_ENABLE_SYSALL = false;
        //V_ENABLE_SYSALL = true;

        bool restore_g_mode_fake_data = g_mode_fake_data;
        g_mode_fake_data = false;

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
        int nch = theFCN.nch;
        int nfp = g_pg.get_number_free_params();
        int ndf = nch - nfp;

        // draw before fit
        draw_input_data drawinputdata;
        drawinputdata.chi2 = fval_before;
        drawinputdata.nch = nch;
        drawinputdata.nfp = nfp;
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
        nch = theFCN.nch;
        nfp = g_pg.get_number_free_params();
        ndf = nch - nfp;

        // draw after fit
        drawinputdata.chi2 = fval_after;
        drawinputdata.nch = nch;
        drawinputdata.nfp = nfp;
        drawinputdata.saveas_filename = "HSD_after";
       
        draw(drawinputdata,
             params_after,
             param_errs_after);

        theParameterStateBefore.Release(std::string(minuit_param_name));

        std::cout << "HSD Fit: NEMO3 Data (same as Summers fit)" << std::endl;
        std::cout << "SYSTEMATICS: CONSTANT OFFSET DISABLED: " << gSystematics.systematic_energy_offset << " MeV" << std::endl;
        std::cout << "SYSTEMATICS: CONSTANT SCALE DISABLED: " << gSystematics.systematic_energy_scale << " MeV" << std::endl;
        std::cout << "SYSTEMATICS: EFFICIENCY DISABLED: " << gSystematics.systematic_efficiency << "" << std::endl;
        std::cout << "SYSTEMATICS: ENRICHMENT DISABLED: " << gSystematics.systematic_enrichment << "" << std::endl;
        std::cout << "SYSTEMATICS: CONSTANT OFFSETSMALL DISABLED: " << gSystematics.systematic_energy_offsetsmall << " MeV" << std::endl;
        std::cout << "Result: " << std::endl;
        std::cout << "fval_before=" << fval_before << std::endl;
        std::cout << "fval_after=" << fval_after
                  << " for params_after[0]=" << params_after[0]
                  << " +- " << param_errs_after[0]
                  << " params_after[1]=" << params_after[1]
                  << " +- " << param_errs_after[1]
                  << std::endl;
        std::cout << std::endl;

        V_ENABLE_SYSALL = restore_V_ENABLE_SYSALL;

        g_mode_fake_data = restore_g_mode_fake_data;

        return 0;
    }
#endif



    // argumnets:
    // CHANNELMODE, SYS_EN, MODE_FAKE, N_FREE, FIX_XI, XI_VAL,
    // min_point, _err, _fval
    // min_point_fname, draw_fname, draw_odir

    ///////////////////////////////////////////////////////////////////////////
    // HSD fixed xi_31 = HSD fit
    // CHANNEL 0 VERSION
    ///////////////////////////////////////////////////////////////////////////
    min_point_helper("CH0", false, false, 1, true, 0.0,
        min_point_data_HSD_CH0, min_point_data_HSD_CH0_err, min_point_data_HSD_CH0_fval,
        "min_point_data_HSD_CH0", "HSD_CH0_data", "HSD_CH0",
        xi_31_systematics_reweight_value);

    ///////////////////////////////////////////////////////////////////////////
    // HSD fixed xi_31 = HSD fit
    ///////////////////////////////////////////////////////////////////////////
    min_point_helper("CH1", false, false, 1, true, 0.0,
        min_point_data_HSD, min_point_data_HSD_err, min_point_data_HSD_fval,
        "min_point_data_HSD", "HSD_data", "HSD",
        xi_31_systematics_reweight_value);

    ///////////////////////////////////////////////////////////////////////////
    // SSD fixed xi_31 = SSD fit
    ///////////////////////////////////////////////////////////////////////////
    min_point_helper("CH1", false, false, 1, true, 0.296,
        min_point_data_SSD, min_point_data_SSD_err, min_point_data_SSD_fval,
        "min_point_data_SSD", "SSD_data", "SSD",
        xi_31_systematics_reweight_value);

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Data
    // Systematics disabled
    ///////////////////////////////////////////////////////////////////////////
    min_point_helper("CH1", false, false, 2, false, -1000.0,
        min_point_data, min_point_data_err, min_point_data_fval,
        "min_point_data", "xifree_data", "xifree",
        xi_31_systematics_reweight_value);

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Data
    // Systematics enabled
    ///////////////////////////////////////////////////////////////////////////
    min_point_helper("CH1", true, false, 2, false, -1000.0,
        min_point_data_SYSALL, min_point_data_SYSALL_err, min_point_data_SYSALL_fval,
        "min_point_data_SYSALL", "xifree_data_SYSALL", "xifree",
        xi_31_systematics_reweight_value);

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // Systematics disabled
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake, min_point_fake_err, min_point_fake_fval,
        "min_point_fake", "xifree_fake", "xifree",
        xi_31_systematics_reweight_value);

    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // Systematics enabled
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", true, true, 2, false, -1000.0,
        min_point_fake_SYSALL, min_point_fake_SYSALL_err, min_point_fake_SYSALL_fval,
        "min_point_fake_SYSALL", "xifree_fake_SYSALL", "xifree",
        xi_31_systematics_reweight_value);

    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // Constant Energy Offset
    // Low Systematic
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    gSystematics.systematic_energy_offset = -0.1;
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake_sysn_l[0], min_point_fake_sysn_l_err[0], min_point_fake_sysn_l_fval[0],
        "min_point_fake_SYS0L", "xifree_fake_SYS0L", "xifree",
        xi_31_systematics_reweight_value);

    ND150_FAKEDATA_SCALE_FACTOR = 1.0;


    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // Constant Energy Offset
    // High Systematic
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    gSystematics.systematic_energy_offset = +0.1;
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake_sysn_h[0], min_point_fake_sysn_h_err[0], min_point_fake_sysn_h_fval[0],
        "min_point_fake_SYS0H", "xifree_fake_SYS0H", "xifree",
        xi_31_systematics_reweight_value);

    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // Energy Scale Multiplier
    // Low Systematic
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    gSystematics.systematic_energy_scale = -0.012;
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake_sysn_l[1], min_point_fake_sysn_l_err[1], min_point_fake_sysn_l_fval[1],
        "min_point_fake_SYS1L", "xifree_fake_SYS1L", "xifree",
        xi_31_systematics_reweight_value);

    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // Energy Scale Multiplier
    // High Systematic
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    gSystematics.systematic_energy_scale = +0.012;
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake_sysn_h[1], min_point_fake_sysn_h_err[1], min_point_fake_sysn_h_fval[1],
        "min_point_fake_SYS1H", "xifree_fake_SYS1H", "xifree",
        xi_31_systematics_reweight_value);

    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // Efficiency
    // Low Systematic
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    gSystematics.systematic_efficiency = -5.55e-02;
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake_sysn_l[2], min_point_fake_sysn_l_err[2], min_point_fake_sysn_l_fval[2],
        "min_point_fake_SYS2L", "xifree_fake_SYS2L", "xifree",
        xi_31_systematics_reweight_value);

    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // Efficiency
    // High Systematic
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    gSystematics.systematic_efficiency = +5.55e-02;
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake_sysn_h[2], min_point_fake_sysn_h_err[2], min_point_fake_sysn_h_fval[2],
        "min_point_fake_SYS2H", "xifree_fake_SYS2H", "xifree",
        xi_31_systematics_reweight_value);

    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // Enrichment
    // Low Systematic
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    gSystematics.systematic_enrichment = -0.5e-02;
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake_sysn_l[3], min_point_fake_sysn_l_err[3], min_point_fake_sysn_l_fval[3],
        "min_point_fake_SYS3L", "xifree_fake_SYS3L", "xifree",
        xi_31_systematics_reweight_value);

    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // Enrichment
    // High Systematic
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    gSystematics.systematic_enrichment = +0.5e-02;
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake_sysn_h[3], min_point_fake_sysn_h_err[3], min_point_fake_sysn_h_fval[3],
        "min_point_fake_SYS3H", "xifree_fake_SYS3H", "xifree",
        xi_31_systematics_reweight_value);
    
    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // Constant Energy Offset Small
    // Low Systematic
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    gSystematics.systematic_energy_offsetsmall = -3.0e-3;
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake_sysn_l[4], min_point_fake_sysn_l_err[4], min_point_fake_sysn_l_fval[4],
        "min_point_fake_SYS4L", "xifree_fake_SYS4L", "xifree",
        xi_31_systematics_reweight_value);

    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // Constant Energy Offset Small
    // High Systematic
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    gSystematics.systematic_energy_offsetsmall = +3.0e-3;
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake_sysn_h[4], min_point_fake_sysn_h_err[4], min_point_fake_sysn_h_fval[4],
        "min_point_fake_SYS4H", "xifree_fake_SYS4H", "xifree",
        xi_31_systematics_reweight_value);

    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // Foil Thickness
    // Low Systematic
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    gSystematics.systematic_foil_thickness_virtual = -1;
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake_sysn_l[5], min_point_fake_sysn_l_err[5], min_point_fake_sysn_l_fval[5],
        "min_point_fake_SYS5L", "xifree_fake_SYS5L", "xifree",
        xi_31_systematics_reweight_value);

    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // Foil Thickness
    // High Systematic
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    gSystematics.systematic_foil_thickness_virtual = +1;
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake_sysn_h[5], min_point_fake_sysn_h_err[5], min_point_fake_sysn_h_fval[5],
        "min_point_fake_SYS5H", "xifree_fake_SYS5H", "xifree",
        xi_31_systematics_reweight_value);

    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // dEdX
    // Low Systematic
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    gSystematics.systematic_dEdX_virtual = -1;
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake_sysn_l[6], min_point_fake_sysn_l_err[6], min_point_fake_sysn_l_fval[6],
        "min_point_fake_SYS6L", "xifree_fake_SYS6L", "xifree",
        xi_31_systematics_reweight_value);

    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // dEdX
    // High Systematic
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    gSystematics.systematic_dEdX_virtual = +1;
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake_sysn_h[6], min_point_fake_sysn_h_err[6], min_point_fake_sysn_h_fval[6],
        "min_point_fake_SYS6H", "xifree_fake_SYS6H", "xifree",
        xi_31_systematics_reweight_value);
    
    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // Bremsstrahlung
    // Low Systematic
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    gSystematics.systematic_brem_virtual = -1;
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake_sysn_l[7], min_point_fake_sysn_l_err[7], min_point_fake_sysn_l_fval[7],
        "min_point_fake_SYS7L", "xifree_fake_SYS7L", "xifree",
        xi_31_systematics_reweight_value);

    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // Bremsstrahlung
    // High Systematic
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    gSystematics.systematic_brem_virtual = +1;
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake_sysn_h[7], min_point_fake_sysn_h_err[7], min_point_fake_sysn_h_fval[7],
        "min_point_fake_SYS7H", "xifree_fake_SYS7H", "xifree",
        xi_31_systematics_reweight_value);

    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // Foil Thickness (nominal)
    // Low Systematic
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    gSystematics.systematic_foil_thickness_nominal = -1;
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake_sysn_l[8], min_point_fake_sysn_l_err[8], min_point_fake_sysn_l_fval[8],
        "min_point_fake_SYS8L", "xifree_fake_SYS8L", "xifree",
        xi_31_systematics_reweight_value);

    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // Foil Thickness (nominal)
    // High Systematic
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    gSystematics.systematic_foil_thickness_nominal = +1;
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake_sysn_h[8], min_point_fake_sysn_h_err[8], min_point_fake_sysn_h_fval[8],
        "min_point_fake_SYS8H", "xifree_fake_SYS8H", "xifree",
        xi_31_systematics_reweight_value);

    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // dEdX (nominal)
    // Low Systematic
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    gSystematics.systematic_dEdX_nominal = -1;
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake_sysn_l[9], min_point_fake_sysn_l_err[9], min_point_fake_sysn_l_fval[9],
        "min_point_fake_SYS9L", "xifree_fake_SYS9L", "xifree",
        xi_31_systematics_reweight_value);

    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // dEdX (nominal)
    // High Systematic
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    gSystematics.systematic_dEdX_nominal = +1;
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake_sysn_h[9], min_point_fake_sysn_h_err[9], min_point_fake_sysn_h_fval[9],
        "min_point_fake_SYS9H", "xifree_fake_SYS9H", "xifree",
        xi_31_systematics_reweight_value);

    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // Bremsstrahlung (nominal)
    // Low Systematic
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    gSystematics.systematic_brem_nominal = -1;
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake_sysn_l[10], min_point_fake_sysn_l_err[10], min_point_fake_sysn_l_fval[10],
        "min_point_fake_SYS10L", "xifree_fake_SYS10L", "xifree",
        xi_31_systematics_reweight_value);

    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // Bremsstrahlung (nominal)
    // High Systematic
    ///////////////////////////////////////////////////////////////////////////
    gSystematics.reset();
    gSystematics.systematic_brem_nominal = +1;
    ND150_FAKEDATA_SCALE_FACTOR = 1.15;

    min_point_helper("CH1", false, true, 2, false, -1000.0,
        min_point_fake_sysn_h[10], min_point_fake_sysn_h_err[10], min_point_fake_sysn_h_fval[10],
        "min_point_fake_SYS10H", "xifree_fake_SYS10H", "xifree",
        xi_31_systematics_reweight_value);

    ND150_FAKEDATA_SCALE_FACTOR = 1.0;

    gSystematics.reset();
    rebuild_fake_data_systematics(xi_31_systematics_reweight_value, xi_31_baseline);




    // 2020-10-05: NOTE
    //
    // Fit to data, set min_point used to be here
    // Moved below, as this function requires the state machine to be in mode
    // V_ENABLE_SYSx = FALSE
    // (ie; systematics disabled)
    //
    // This means that the first call to MinimizeFCNAxialVector::operator()
    // has V_ENABLE_SYSx = false, and therefore the V_PHYS_MATHMORE objects
    // are initilized with systematics DISABLED.
    //
    // Now this call block has been moved to the end (before MPS block)
    // the V_PHYS_MATHMORE objects will contain systematics. I do not know
    // if this will affect the result of the calculations of the min_point
    //
    // Probably it will not?
    //
    // Further note: This still does not resolve the problem, as the first
    // call to MinimizeFCNAxialVector::operator() has the following settings:
    //
    // V_ENABLE_SYSALL = false;
    // V_ENABLE_SYS1 = false;
    // V_ENABLE_SYS2 = false;
    // V_ENABLE_SYS3 = false;
    // V_ENABLE_SYS4 = false;
    //
    //











    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
















    return 0;






    // TODO: ll_walk is saved from the previous fit (call to exec)
    // which in this case is the high systematic energy offset fit




    ll_walk_save.clear();

    // reenable this one
    if(1)
    {
        
        //DRAW_V_PHYS_MATRIX = false;

        //std::cout << "recalculate_V_PHYS_MATHMORE=" << recalculate_V_PHYS_MATHMORE << std::endl;
        //std::cout << "recalculate_V_PHYS_SYS=" << recalculate_V_PHYS_SYS << std::endl;
        //recalculate_V_PHYS_MATHMORE = true;
        recalculate_V_PHYS_SYS = true;
        //std::cout << "recalculate_V_PHYS_MATHMORE=" << recalculate_V_PHYS_MATHMORE << std::endl;
        //std::cout << "recalculate_V_PHYS_SYS=" << recalculate_V_PHYS_SYS << std::endl;

        std::cout << "V_ENABLE_SYSALL=" << V_ENABLE_SYSALL << std::endl;
        for(int i = 0; i < N_SYSTEMATICS; ++ i)
        {
            std::cout << "V_ENABLE_SYS" << i << "=" << V_ENABLE_SYSn[i] << std::endl;
        }
        std::cout << "g_mode_fake_data=" << g_mode_fake_data << std::endl;

        gSystematics.reset();
        gSystematics_print();

        ND150_FAKEDATA_SCALE_FACTOR = 1.15;
        rebuild_fake_data_systematics(xi_31_systematics_reweight_value, xi_31_baseline);
        ND150_FAKEDATA_SCALE_FACTOR = 1.0;

        std::cout << "ready to test MPS" << std::endl;
        std::cout << "output_name=" << output_name << std::endl;
        std::cout << "start_index=" << start_index << std::endl;
        std::cout << "stop_index=" << stop_index << std::endl;
        std::cout << "calling calculate" << std::endl;
        newloglikfitter_mps_calculate(number_job_id, output_name, start_index, stop_index);

        std::cout << "MPS done" << std::endl;
    }
   


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

    
    //std::cout << "ok to run chisquare test? [press enter to continue]" << std::endl;
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

    std::chrono::system_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> runtime_sec = end_time - start_time;
    std::cout << "execution time: " <<  runtime_sec.count() / 3600.0 << " hours" << std::endl;
}























    ///////////////////////////////////////////////////////////////////////////
    // All Parameter Fit
    // Fake Data
    // Variable Systematic Parameter
    ///////////////////////////////////////////////////////////////////////////

#if 0
    
    bool g_mode_fake_data_restore_value = g_mode_fake_data;
    g_mode_fake_data = true;
    
    TCanvas *results_c_xi_31 = nullptr;
    TCanvas *results_c_A = nullptr;
    TCanvas *results_c_chi2_before = nullptr;
    TCanvas *results_c_chi2_after = nullptr;

    std::vector<TCanvas*> results_c_allparams;

    // do not do this in parallel mode
    if(0) // || (MODE_PARALLEL == 0))
    {
        // TODO: enable/disable SYSX?

        std::vector<double> results_x_xi_31;
        std::vector<double> results_y_xi_31;
        std::vector<double> results_x_A;
        std::vector<double> results_y_A;
        std::vector<double> results_x_chi2_before;
        std::vector<double> results_y_chi2_before;
        std::vector<double> results_x_chi2_after;
        std::vector<double> results_y_chi2_after;

        std::vector<std::vector<double>> results_y_allparams;

        const int i_max = 0; //10; //0;
        for(int i = 0; i <= i_max; ++ i)
        {

            const double min = -0.01;
            const double max = +0.01;
            const double diff = max - min;
            const double fr = (double)i / (double)i_max;
            //gSystematics.systematic_energy_offset = fr * diff + min;
            gSystematics.systematic_energy_offset = 0.0;
            gSystematics.systematic_energy_scale = 0.0;
            double systematic_energy_offset = gSystematics.systematic_energy_offset;
            double systematic_energy_scale = gSystematics.systematic_energy_scale;
            std::cout << "seo=" << systematic_energy_offset << std::endl;
            std::cout << "sem=" << systematic_energy_scale << std::endl;
            //rebuild_fake_data_systematics(0.296, xi_31_baseline); // want to check if the fitter can fit itself to itself
            rebuild_fake_data_systematics(0.0, xi_31_baseline); // want to check if the fitter can fit itself to itself
            // just check the output looks sensible

            std::string name_extra = "seo_" + std::to_string(systematic_energy_offset)
                                   + "_sem_" + std::to_string(systematic_energy_scale);

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
            int nch = theFCN.nch;
            int nfp = g_pg.get_number_free_params();
            int ndf = nch - nfp;

            // initilize output vectors
            if(i == 0)
            {
                for(int i = 0; i < g_pg.file_params.size(); ++ i)
                {
                    std::vector<double> blankvector;
                    results_y_allparams.push_back(blankvector);
                }
            }

            // draw before fit
            draw_input_data drawinputdata;
            drawinputdata.chi2 = fval_before;
            drawinputdata.nch = nch;
            drawinputdata.nfp = nfp;
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
            nch = theFCN.nch;
            nfp = g_pg.get_number_free_params();
            ndf = nch - nfp;

            // draw result
            drawinputdata.chi2 = fval_after;
            drawinputdata.nch = nch;
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

            results_x_xi_31.push_back(systematic_energy_offset);
            results_y_xi_31.push_back(params_after.at(1));

            results_x_A.push_back(systematic_energy_offset);
            results_y_A.push_back(params_after.at(0));
            
            results_x_chi2_before.push_back(systematic_energy_offset);
            results_y_chi2_before.push_back(fval_before);

            results_x_chi2_after.push_back(systematic_energy_offset);
            results_y_chi2_after.push_back(fval_after);

            for(int i = 0; i < g_pg.file_params.size(); ++ i)
            {
                results_y_allparams.at(i).push_back(params_after.at(i));
            }
        }

        TGraph *results_g_xi_31 = new TGraph(results_x_xi_31.size(), results_x_xi_31.data(), results_y_xi_31.data());
        TGraph *results_g_A = new TGraph(results_x_A.size(), results_x_A.data(), results_y_A.data());
        TGraph *results_g_chi2_before = new TGraph(results_x_chi2_before.size(), results_x_chi2_before.data(), results_y_chi2_before.data());
        TGraph *results_g_chi2_after = new TGraph(results_x_chi2_after.size(), results_x_chi2_after.data(), results_y_chi2_after.data());

        
        for(int i = 0; i < g_pg.file_params.size(); ++ i)
        {
            std::cout << "i=" << i << std::endl;

            TGraph *results_g_allparams =
                new TGraph(results_x_xi_31.size(),
                           results_x_xi_31.data(),
                           results_y_allparams.at(i).data());

            TString cname;
            cname.Form("results_P%d", i);
            std::cout << "cname=" << cname << std::endl;
            results_c_allparams.push_back(new TCanvas(cname, cname));
            results_g_allparams->Draw();
        }


        results_c_xi_31 = new TCanvas("results_xi_31", "results_xi_31");
        results_g_xi_31->Draw();

        results_c_A = new TCanvas("results_A", "results_A");
        results_g_A->Draw();

        results_c_chi2_before = new TCanvas("results_chi2_before", "results_chi2_before");
        results_g_chi2_before->Draw();

        results_c_chi2_after = new TCanvas("results_chi2_after", "results_chi2_after");
        results_g_chi2_after->Draw();
    }

    g_mode_fake_data = g_mode_fake_data_restore_value;
#endif

