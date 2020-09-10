#ifndef MINIMIZEFCNAXIALVECTOR_H
#define MINIMIZEFCNAXIALVECTOR_H



// Root headers, minuit
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"






//Double_t getNumberMC1D(const Int_t channel, const Int_t bin_ix, const std::vector<double> &p);
//Double_t getNumberMC2D(const Int_t channel, const Int_t bin_ix, const Int_t bin_iy, const std::vector<double> &p);



// TODO: may need to change precision to float not double using
// MnMachinePrecision::setPrecision

class MinimizeFCNAxialVector : public ROOT::Minuit2::FCNBase
{


    public:


    MinimizeFCNAxialVector()
        : error_def{1.0} // LL // TODO: check as I am actually minimizing chisquare
    {
    }

    virtual
    ~MinimizeFCNAxialVector()
    {
    }

    virtual
    double Up() const
    {
        return error_def;
    }

    virtual
    double operator()(const std::vector<double> &param) const
    {


        const int MODE_LOGPOISSON = 0;
        const int MODE_CHI2 = 1;
        //const int MODE_METRIC = MODE_LOGPOISSON;
        const int MODE_METRIC = MODE_CHI2;

        const int debuglevel = 1;
        


        // TODO: subtract numfreepars
        ndf = 0;
        int ndf_P1 = 0;
        int ndf_P2 = 0;

        /*
        std::cout << "operator()-> param: ";
        for(int pix = 0; pix < param.size(); ++ pix)
        {
            std::cout << param[pix] << " ";
        }
        std::cout << std::endl;
        */
       
        Int_t nPar = param.size();
        Double_t fval = 0.0;

        // TODO: need to re-enable to use MODE fake data
        /*
        if(g_mode_fake_data == true)
        {
            if(allFakeDataSamples1D == nullptr)
            {
                if(debuglevel > 0)
                {
                    std::cout << "calling build_fake_data()" << std::endl;
                }
                std::cin.get();
                build_fake_data();
            }
        }
        */
        // TODO: manually call this from main function (loadFiles)


        // save the walk
        ll_walk.push_back(std::make_pair(param.at(1), param.at(0)));


        // draw the output
        //TString fname;
        //fname.Form("lliter_%d", counter);
        //draw_channel(1, p, std::string(fname));



        // error mode
        // 1 = data
        // 2 = MC
        // 3 = both in quadrature
        const int EMODE = 1; // only use errors from data:
        // reason:
        // what do errors represent? uncertainty on a sample
        // in the case of real or fake data, the uncertainty is
        // sqrt(N) for each bin, because this is the statistical uncertainty
        // in a counting experiment (not really, it's actually poisson but
        // gaussian approximation is good enough)
        // and there is no uncertainty on the MC because we assume infinite
        // statistics for these samples (not quite true but number of events
        // is much higher than number of data events) - we rescale the MC
        // to the required activity for the livetime of the experiment but
        // regardless of this sqrt(N) would be rescaled accordingly as well
        // and N for MC >> N for data, adding in the rescaling factor we find
        // that the MC uncertainty is even smaller


        const int xi_31_ext_param_number = g_pg.get_xi_31_ext_param_number(); 


        if(debuglevel >= 2)
        {
            std::cout << std::scientific;
            std::cout << "logLikelihood"
                      << " param[0]=" << param[0]
                      << " param[" << xi_31_ext_param_number << "]="
                      << param[xi_31_ext_param_number] << std::endl;
        }


        /*
        if(systematic_energy_offset != systematic_energy_offset_last)
        {
            systematic_energy_offset_last = systematic_energy_offset;
            
            if(debuglevel >= 2)
            {
                std::cout << __func__ << " rebuild_150Nd_data() -> systematic_energy_offset=" << systematic_energy_offset << std::endl;
            }

            rebuild_150Nd_data();
        }
        */


        // TODO: will not work if parameter number changes
        // 2020-06-17
        //if(p[1] != last_xi_31_parameter_value)
        //if(p[1] != paramLastValueMap[1])
        //if(param[axial_vector_parameter_0_param_number] != paramLastValueMap[axial_vector_parameter_0_param_number])
        if(param[xi_31_ext_param_number] != g_pg.file_params.at(xi_31_ext_param_number).paramLastValue)
        {

            // TODO: rebuild nd150 xi_31 paramter histogram here
            //std::cout << "rebuilding 150 Nd MC" << std::endl;

            const double xi_31{param[xi_31_ext_param_number]};
            if(debuglevel >= 3)
            {
                std::cout << "rebuild 150Nd MC" << std::endl;
                std::cout << "param[" << xi_31_ext_param_number << "]=" << param[xi_31_ext_param_number] << " != " << g_pg.file_params.at(xi_31_ext_param_number).paramLastValue << std::endl;
            }
            rebuild_150Nd_MC(xi_31, xi_31_baseline);


    /*
            TCanvas *ctmp = new TCanvas("ctmp", "ctmp");
            hSingleEnergy->Draw();
            TString fname;
            fname.Form("ctmp_%d.png", counter);
            ctmp->SaveAs(fname);
            ++ counter;
    */


            ///////////////////////////////////////////////////////////////////////
            // set last value
            ///////////////////////////////////////////////////////////////////////
            

            //last_xi_31_parameter_value = xi_31;
            // TODO: check that I set .paramLastValue in g_pg.file_params later at end of this function
            // if so can remove this statement
            g_pg.file_params.at(xi_31_ext_param_number).paramLastValue = xi_31;
        }
        else
        {
            if(debuglevel >= 4)
            {
                std::cout << "dont rebuild 150Nd MC" << std::endl;
            }
        }

        // TODO: add check here to see if any disabled parameters are accessed


        bool new_method = true;
        //bool new_method = false;
        if(new_method)
        {

            ///////////////////////////////////////////////////////////////////////
            // loglikelihood/chisquare function
            // new method using V Matrix
            ///////////////////////////////////////////////////////////////////////


            ///////////////////////////////////////////////////////////////////////
            // allocation

            // number of bins: dictated by
            // 1: number of channels
            // 2: number of phases (always 2)
            // 3: number of bins in histograms
            //    (always 50 for 1D, and 50 * 50 for 2D)
            // therefore total number of bins is:
            // number of 1D channels * 2 phases * 50 bins per hist
            // + number of 2D channels * 2 phases * 50*50 bins per hist
            //
            // How to calculate the index, assuming the current local bin
            // index is bin_x
            //
            // for 1D:
            // for Phase 1:
            // super_index = channel * 2 * 50 + bin_x;
            // for Phase 2:
            // super_index = channel * 2 * 50 + 50 + bin_ix
            // for 2D:
            // for Phase 1:
            // super_index = number1DChannels * 2 * 50 + channel * 2 * 50 * 50 + 50 * bin_x + bin_x
            // for Phase 2:
            // super_index = number1DChannels * 2 * 50 + channel * 2 * 50 * 50 + 50 * 50 + 50 * bin_y + bin_x
            const Int_t NUM_BINS_XY = 50 * 2 * number1DHists + 50 * 50 * 2 * number2DHists;
            if(V_CHEN == nullptr)
            {
                std::cout << "Alloc V_CHEN" << std::endl;
                V_CHEN = new TH2D("V_CHEN", "V_CHEN",
                                  NUM_BINS_XY, 0.0, NUM_BINS_XY,
                                  NUM_BINS_XY, 0.0, NUM_BINS_XY);
            }
            if(V_PHYS == nullptr)
            {
                std::cout << "Alloc V_PHYS" << std::endl;
                V_PHYS = new TH2D("V_PHYS", "V_PHYS",
                                  NUM_BINS_XY, 0.0, NUM_BINS_XY,
                                  NUM_BINS_XY, 0.0, NUM_BINS_XY);
            }
            if(V_PHYS_STAT == nullptr)
            {
                std::cout << "Alloc V_PHYS_STAT" << std::endl;
                V_PHYS_STAT = new TH2D("V_PHYS_STAT", "V_PHYS_STAT",
                                  NUM_BINS_XY, 0.0, NUM_BINS_XY,
                                  NUM_BINS_XY, 0.0, NUM_BINS_XY);
            }
            if(V_PHYS_SYS == nullptr)
            {
                std::cout << "Alloc V_PHYS_SYS" << std::endl;
                V_PHYS_SYS = new TH2D("V_PHYS_SYS", "V_PHYS_SYS",
                                  NUM_BINS_XY, 0.0, NUM_BINS_XY,
                                  NUM_BINS_XY, 0.0, NUM_BINS_XY);
            }
            if(D == nullptr)
            {
                std::cout << "Alloc D" << std::endl;
                D = new TH2D("D", "D",
                                  NUM_BINS_XY, 0.0, NUM_BINS_XY,
                                  1, 0.0, 1.0);

                // clear the contents of D
                for(Int_t bin_i{1}; bin_i <= D->GetNbinsX(); ++ bin_i)
                {
                    const Double_t zero = 0.0;
                    D->SetBinContent(bin_i, 1, zero);
                }
            }
            if(M == nullptr)
            {
                std::cout << "Alloc M" << std::endl;
                M = new TH2D("M", "M",
                                  NUM_BINS_XY, 0.0, NUM_BINS_XY,
                                  1, 0.0, 1.0);

                // clear the contents of M
                for(Int_t bin_i{1}; bin_i <= M->GetNbinsX(); ++ bin_i)
                {
                    const Double_t zero = 0.0;
                    M->SetBinContent(bin_i, 1, zero);
                }
            }
            if(D_minus_M == nullptr)
            {
                std::cout << "Alloc D_minus_M" << std::endl;
                D_minus_M = new TH2D("D_minus_M", "D_minus_M",
                                  NUM_BINS_XY, 0.0, NUM_BINS_XY,
                                  1, 0.0, 1.0);

                // clear the contents of M
                for(Int_t bin_i{1}; bin_i <= D_minus_M->GetNbinsX(); ++ bin_i)
                {
                    const Double_t zero = 0.0;
                    D_minus_M->SetBinContent(bin_i, 1, zero);
                }
            }
            if(V_SUPER == nullptr)
            {
                std::cout << "Alloc V_SUPER" << std::endl;
                V_SUPER = new TH2D("V_SUPER", "V_SUPER",
                                  NUM_BINS_XY, 0.0, NUM_BINS_XY,
                                  NUM_BINS_XY, 0.0, NUM_BINS_XY);
            }


            ///////////////////////////////////////////////////////////////////////
            // set values for M, D and D_minus_M

            // can D change?
            // yes - if systematic parameters change and we are in fake data mode
            // just assume it always changes for now
            
            // set the contents of D

            // loop over all channels
            for(int channel = 0; channel < number1DHists; ++ channel)
            {
                if(debuglevel >= 4)
                {
                    std::cout << "channel=" << channel << std::endl;
                }


                // TODO: this no longer works because there are multiple phases
                // in the tmpData1D arrays
                // copy code from below to find correct object
                // TODO: phase 1 and phase 2 data objects
                //Double_t *tmpData1D_P1 = nullptr; //(TH1D*)allDataSamples1D->At(channel);
                //Double_t *tmpData1D_P2 = nullptr;
                //Double_t *tmpFakeData1D_P1 = nullptr; //(TH1D*)allFakeDataSamples1D->At(channel);
                //Double_t *tmpFakeData1D_P2 = nullptr;

                std::string histname = std::string(channel_histname_1D[channel]);
                std::string search_object_P1;
                std::string search_object_P2;
                if(g_mode_fake_data == false)
                {
                    search_object_P1 = histname + std::string(DataFile) + "_P1";
                    search_object_P2 = histname + std::string(DataFile) + "_P2";
                }
                else if(g_mode_fake_data == true)
                {
                    search_object_P1 = histname + std::string("fakedata") + "_P1";
                    search_object_P2 = histname + std::string("fakedata") + "_P2";
                }
                TH1D *tmpDataHist1D_P1 = nullptr;
                TH1D *tmpDataHist1D_P2 = nullptr;
                
                if(debuglevel >= 6)
                {
                    std::cout << "search_object_P1=" << search_object_P1
                              << " search_object_P2=" << search_object_P2 << std::endl;
                }

                if(g_mode_fake_data == false)
                {
                    tmpDataHist1D_P1 = (TH1D*)allDataSamples1D->FindObject(search_object_P1.c_str());
                    tmpDataHist1D_P2 = (TH1D*)allDataSamples1D->FindObject(search_object_P2.c_str());
                }
                else if(g_mode_fake_data == true)
                {
                    tmpDataHist1D_P1 = (TH1D*)allFakeDataSamples1D->FindObject(search_object_P1.c_str());
                    tmpDataHist1D_P2 = (TH1D*)allFakeDataSamples1D->FindObject(search_object_P2.c_str());
                }

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

                //tmpData1D_P1 = new Double_t[tmpDataHist1D_P1->GetNbinsX()]; //tmpHistP1->Clone("tmpData1D_P1");
                //tmpData1D_P2 = new Double_t[tmpDataHist1D_P2->GetNbinsX()]; //tmpHistP2->Clone("tmpData1D_P2");
                for(Int_t bin_x{0}; bin_x < tmpDataHist1D_P1->GetNbinsX(); ++ bin_x)
                {
                    //tmpData1D_P1[bin_x] = tmpDataHist1D_P1->GetBinContent(bin_x);
                    Int_t super_index = channel * 2 * 50 + bin_x;
                    Double_t content_input = 0.0;
                    Double_t content_add = tmpDataHist1D_P1->GetBinContent(bin_x + 1);
                    Double_t content_output = content_input + content_add;
                    D->SetBinContent(super_index + 1, 1, content_output);
                }
                
                //std::cout << "LIST OF P2 DATA channel=" << channel << std::endl;
                //std::cin.get();
                for(Int_t bin_x{0}; bin_x < tmpDataHist1D_P2->GetNbinsX(); ++ bin_x)
                {
                    //tmpData1D_P2[bin_x] = tmpDataHist1D_P2->GetBinContent(bin_x);
                    Int_t super_index = channel * 2 * 50 + 50 + bin_x;
                    Double_t content_input = 0.0;
                    Double_t content_add = tmpDataHist1D_P2->GetBinContent(bin_x + 1);
                    Double_t content_output = content_input + content_add;
                    D->SetBinContent(super_index + 1, 1, content_output);
                    //std::cout << "bin_x=" << bin_x + 1 << " content=" << content_output << std::endl;
                }

            }


            // can M change?
            // yes - if xi_31 parameter changes, which is detectable using g_pg
            // as shown above
            // can it change another way?
            // it depends on all "amplitude" parameters, so need to rebuild if any
            // of these parameters changes
            // just assume it always changes for now


            // loop over all channels
            for(int channel = 0; channel < number1DHists; ++ channel)
            {
                if(debuglevel >= 4)
                {
                    std::cout << "channel=" << channel << std::endl;
                }

                /*
                Double_t *tmpTotalMC1D_P1 = new Double_t[tmpDataHist1D_P1->GetNbinsX()];
                Double_t *tmpTotalMC1D_P2 = new Double_t[tmpDataHist1D_P2->GetNbinsX()];
                for(Int_t bin_x{1}; bin_x <= tmpDataHist1D_P1->GetNbinsX(); ++ bin_x)
                {
                    tmpTotalMC1D_P1[bin_x] = 0.0;
                }
                for(Int_t bin_x{1}; bin_x <= tmpDataHist1D_P2->GetNbinsX(); ++ bin_x)
                {
                    tmpTotalMC1D_P2[bin_x] = 0.0;
                }
                */

                // loop over all the parameters
                std::map<int, file_parameter>::iterator it{g_pg.file_params.begin()};
                for(; it != g_pg.file_params.end(); ++ it)
                {
                    int paramNumberInt = -1;

                    int paramNumber = it->second.paramNumber;
                    std::cout << "paramNumber=" << paramNumber << std::endl;
                    bool paramEnabled = it->second.paramEnabled;
                    bool paramEnabledP1 = it->second.paramEnabledP1;
                    bool paramEnabledP2 = it->second.paramEnabledP2;
                    double paramInitValue = it->second.paramInitValue;
                    double paramInitError = it->second.paramInitError;
                    int paramConstraintMode = it->second.paramConstraintMode;

                    if(debuglevel >= 5)
                    {
                        std::cout << "paramNumber=" << paramNumber << std::endl;
                    }

                    // NOTE: notes on the V MATRIX method and the enable/disable
                    // of parameters
                    //
                    // 1: how can parameters be enabled/disabled
                    //
                    // a:
                    // There are two sources, the first is the
                    // parameter_names.lst file, which has options to
                    // enable/disable the parameter totally, or enable/disable
                    // it for P1 and P2
                    // 
                    // b:
                    // The other source is the channel enable/disable switch
                    // in the header file
                    // This affects the V_CHEN matrix because it switches
                    // on/off groups of BINS in the chisquare calculation
                    //
                    // By contrast, (a) switches on and off parameters
                    // which does not affect the number of bins only the bin
                    // content

                    if(paramEnabled == false)
                    {
                        if(debuglevel > 0)
                        {
                            std::cout << __func__ << " channel [channel TODO] is disabled (overall) skip" << std::endl;
                            std::cin.get();
                        }
                        continue;
                    }
                    if((paramEnabledP1 || paramEnabledP2) == false)
                    {
                        // if both are false
                        if(debuglevel > 0)
                        {
                            std::cout << __func__ << " channel [channel TODO] is disabled for P1 and P2, skip" << std::endl;
                            std::cin.get();
                        }
                        continue;
                    }
                    /*
                    bool ok = false;
                    if(paramEnabled == true)
                    {
                        if(gEnablePhase1 == true)
                        {
                            if(paramEnabledP1 == true)
                            {
                                ok = true;
                            }
                        }

                        if(gEnablePhase2 == true)
                        {
                            if(paramEnabledP2 == true)
                            {
                                ok = true;
                            }
                        }
                    }
                    else
                    {
                        ok = true;
                        continue;
                    }
                    if(ok == false)
                    {
                        if(debuglevel > 0)
                        {
                            std::cout << __func__ << " ok == false" << std::endl;
                            std::cin.get();
                        }
                        continue;
                    }
                    */
                    // TODO: this correctly ignores any parameter which is disabled such that
                    // paramEnabled == false
                    // however, it also ignores parameters which are set as disabled for P1
                    // and P2, and when these are irrelevent due to the value of gEnablePhaseX
                    // so... the internal and external index will not match
                    // need to add some code to fix this when the parameters are read from
                    // file, (probably)
                    // unless I just ignore that here... perhaps paramEnabled dictates
                    // whether parameter is drawn and the phase1/phase2 enable flag
                    // is to decide whether minuit does the fit or not (in which case
                    // the param may still contribute to chisquare but may not be minimized
                    // by minuit)


                    // looping over parameters...
                    // get a paramNumber
                    // corresponds to a paramNumberInternal (minuit param number)
                    // which is the index of the vector param where the current
                    // activity can be found
                    // also get a paramName, which is something like ac228_bi212_tl208_int
                    // which is not very useful
                    // but get a list of mc names, which for this example will be
                    // ac228_int_rot, bi212_int_rot, tl208_int_rot
                    // the corresponding MC samples names will be
                    // theHistogram + BkgFiles[i] + "_Px_fit"
                    // aka:
                    // hTotalE_ + bi212_int_rot + _P1_fit
                    // hTotalE_ + bi212_int_rot + _P2_fit
                    // now, the paramNumber can be used to look up the parameter in
                    // g_pg, from which we can obtain paramEnabled, paramEnabledP1
                    // and paramEnabledP2, as well as the constraint mode, constraint
                    // value and error
                    // so, have to strip off the "_fit" ending, and strip off the
                    // histogram (channel) name from the front, aka remove "hTotalE_"
                    // then left with something like
                    // bi212_int_rot_P1 and a channel number which comes from the loop
                    // over allDataSamples1D->GetEntries()
                    // EXCEPT this will now be wrong, because we have data for P1 and P2
                    // so what we want to loop over is an index channel which runs from
                    // 0 to number1DHists
                    // so we have a channel index which is used to select the channel
                    // which selects the index of allMCSamples1D[channel]
                    // we then have to somehow "get" from this TObject
                    // the histogram with the name bi212_int_rot_P1
                    // and bi212_int_rot_P2
                    // once we have done that we switch on paramEnabled, and paramEnabledP1
                    // and paramEnabledP2 to decide whether or not we skip (in the case
                    // where paramEnabled == false), apply a scaling with param[index]
                    // for P1 and P2 if they are enabled, or simply add a value
                    // where the scaling factor is given by paramInitValue for P1 and P2
                    // if they are disabled
                    // TODO: does this last part make sense? no it doesn't. skip if
                    // P1/P2 not enabled

                    std::vector<std::string>::iterator mc_name_it{it->second.MCNameList.begin()};
                    for(; mc_name_it != it->second.MCNameList.end(); ++ mc_name_it)
                    {
                        std::string mc_name = *mc_name_it;
                        std::string histname = std::string(channel_histname_1D[channel]);
                        std::string search_object_P1 = histname + mc_name + "_P1_fit";
                        std::string search_object_P2 = histname + mc_name + "_P2_fit";
                        TH1D *tmpHist1D_P1 = nullptr;
                        TH1D *tmpHist1D_P2 = nullptr;
        
                        if(debuglevel >= 6)
                        {
                            std::cout << "search_object_P1=" << search_object_P1
                                      << " search_object_P2=" << search_object_P2 << std::endl;
                        }

                        paramNumberInt = g_pg.ExtToIntParamNumberMap.at(paramNumber);
                        if(debuglevel >= 5)
                        {
                            std::cout << "paramNumber=" << paramNumber << " -> " << paramNumberInt << std::endl;
                        }

                        // phase 1
                        if(paramEnabledP1 == true)
                        {
                            tmpHist1D_P1 = (TH1D*)allMCSamples1D[channel]->FindObject(search_object_P1.c_str());

                            if(tmpHist1D_P1 == nullptr)
                            {
                                std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
                                throw "problem";
                            }

                            Double_t scale_factor_P1 = param.at(paramNumberInt);
                            if(debuglevel >= 6)
                            {
                                std::cout << "enabled P1: scale factor P1: " << scale_factor_P1 << std::endl;
                            }

                            for(Int_t bin_x{0}; bin_x < tmpHist1D_P1->GetNbinsX(); ++ bin_x)
                            {
                                Int_t super_index = channel * 2 * 50 + bin_x;
                                Double_t content_input = M->GetBinContent(super_index + 1, 1);
                                Double_t content_add = scale_factor_P1 * tmpHist1D_P1->GetBinContent(bin_x + 1);
                                Double_t content_output = content_input + content_add;
                                M->SetBinContent(super_index + 1, 1, content_output);
                                //std::cout << "content_input=" << content_input << " content_output=" << content_output << " content_add=" << content_add << std::endl;
                            }
                        }
                        else
                        {
                            if(debuglevel >= 4)
                            {
                                std::cout << "disabled P1" << std::endl;
                            }
                        }

                        // phase 2
                        if(paramEnabledP2 == true)
                        {
                            tmpHist1D_P2 = (TH1D*)allMCSamples1D[channel]->FindObject(search_object_P2.c_str());
                        
                            if(tmpHist1D_P2 == nullptr)
                            {
                                std::cout << "ERROR: Could not find object " << search_object_P2 << std::endl;
                                throw "problem";
                            }

                            Double_t scale_factor_P2 = param.at(paramNumberInt);
                            if(debuglevel >= 4)
                            {
                                std::cout << "enabled P2: scale factor P2: " << scale_factor_P2 << std::endl;
                            }

                            //std::cout << "super_index (start) : " << channel * 2 * 50 + 50 + 0 << std::endl;
                            for(Int_t bin_x{0}; bin_x < tmpHist1D_P2->GetNbinsX(); ++ bin_x)
                            {
                                Int_t super_index = channel * 2 * 50 + 50 + bin_x;
                                Double_t content_input = M->GetBinContent(super_index + 1, 1);
                                Double_t content_add = scale_factor_P2 * tmpHist1D_P2->GetBinContent(bin_x + 1);
                                Double_t content_output = content_input + content_add;
                                M->SetBinContent(super_index + 1, 1, content_output);
                                //std::cout << "debug: " << "super_index=" << super_index << " content_input=" << content_input << " content_add=" << content_add << " content_output=" << content_output << " M:" << M->GetBinContent(super_index + 1, 1) << std::endl;

                                if(channel == 1)
                                {
                                    //if(mc_name.find("tl208") != std::string::npos)
                                    //{
                                        if(bin_x + 1 == 13)
                                        {
                                            std::cout << "mc_name=" << mc_name << " content_add, bin 13=" << content_add << std::endl;
                                        }
                                    //}
                                }
                            }
                        }
                        else
                        {
                            if(debuglevel >= 4)
                            {
                                std::cout << "disabled P2" << std::endl;
                            }
                        }




                        // TODO: what is parameter is NOT enabled

                    } // mc sample name iterator

                } // file_param iterator

                std::cout << "LIST OF P2 DATA (construting M here) channel=" << channel << std::endl;
                std::cout << "super_index (start, below) : " << channel * 2 * 50 + 50 + 0 << std::endl;
                for(Int_t bin_x{0}; bin_x < 50; ++ bin_x)
                {
                    //tmpData1D_P2[bin_x] = tmpDataHist1D_P2->GetBinContent(bin_x);
                    Int_t super_index = channel * 2 * 50 + 50 + bin_x;
                    Double_t content_add = M->GetBinContent(super_index + 1, 1);
                    //std::cout << "channel=" << channel << " super_index=" << super_index << " bin_x=" << bin_x + 1 << " content=" << content_add << std::endl;
                }
                //std::cin.get();

            } // channel

            // draw
            /*
            TCanvas *c_D = new TCanvas("c_D", "c_D", 10000, 10000);
            D->Draw("colz");
            c_D->SaveAs("debug_c_D.png");
            TCanvas *c_M = new TCanvas("c_M", "c_M", 10000, 10000);
            M->Draw("colz");
            c_M->SaveAs("debug_c_M.png");
            */

            // M is set, D is set

            // set D_minus_M
            for(Int_t binx{1}; binx <= D->GetNbinsX(); ++ binx)
            {
                Double_t content_D = D->GetBinContent(binx, 1);
                Double_t content_M = M->GetBinContent(binx, 1);
                if(content_D != 0.0)
                {
                    std::cout << "binx=" << binx << " ~> " << content_D << " " << content_M << std::endl;
                }
                Double_t content_D_minus_M = content_D - content_M;
                D_minus_M->SetBinContent(binx, 1, content_D_minus_M);
            }

            // draw
            
            TCanvas *c_D_minus_M = new TCanvas("c_D_minus_M", "c_D_minus_M", 10000, 10000);
            D_minus_M->Draw("colz");
            c_D_minus_M->SaveAs("debug_c_D_minus_M.png");
            
            std::cin.get();
            

            // set V_PHYS_STAT
            for(Int_t binx{1}; binx <= M->GetNbinsX(); ++ binx)
            {
                Double_t content = M->GetBinContent(binx);
            }

        }
        else // if not new_method
        {

            ///////////////////////////////////////////////////////////////////////////
            // loglikelihood, 1D channels
            ///////////////////////////////////////////////////////////////////////////


            double loglik = 0.0;

        //   std::cout << "getting 1D histograms" << std::endl;

            // std::cout << allDataSamples1D->GetEntries()  << std::endl;

            // loop over all channels
            for(int channel = 0; channel < number1DHists; ++ channel)
            {
                if(debuglevel >= 4)
                {
                    std::cout << "channel=" << channel << std::endl;
                }


                // check channel enabled
                if(channel_enable_1D[channel] == 0)
                {
                    if(debuglevel >= 5)
                    {
                        std::cout << "1D: channel " << channel << " disabled, skip" << std::endl;
                    }
                    continue;
                }

                // LL value for this channel
                double ll_channel_P1 = 0.0;
                double ll_channel_P2 = 0.0;

                // TODO: what is this for and is it required?
                TString channel_str;
                channel_str.Form("%i", channel);

                // TODO: this no longer works because there are multiple phases
                // in the tmpData1D arrays
                // copy code from below to find correct object
                // TODO: phase 1 and phase 2 data objects
                Double_t *tmpData1D_P1 = nullptr; //(TH1D*)allDataSamples1D->At(channel);
                Double_t *tmpData1D_P2 = nullptr;
                //Double_t *tmpFakeData1D_P1 = nullptr; //(TH1D*)allFakeDataSamples1D->At(channel);
                //Double_t *tmpFakeData1D_P2 = nullptr;

                std::string histname = std::string(channel_histname_1D[channel]);
                std::string search_object_P1;
                std::string search_object_P2;
                if(g_mode_fake_data == false)
                {
                    search_object_P1 = histname + std::string(DataFile) + "_P1";
                    search_object_P2 = histname + std::string(DataFile) + "_P2";
                }
                else if(g_mode_fake_data == true)
                {
                    search_object_P1 = histname + std::string("fakedata") + "_P1";
                    search_object_P2 = histname + std::string("fakedata") + "_P2";
                }
                TH1D *tmpDataHist1D_P1 = nullptr;
                TH1D *tmpDataHist1D_P2 = nullptr;
                
                if(debuglevel >= 6)
                {
                    std::cout << "search_object_P1=" << search_object_P1
                              << " search_object_P2=" << search_object_P2 << std::endl;
                }

                if(g_mode_fake_data == false)
                {
                    tmpDataHist1D_P1 = (TH1D*)allDataSamples1D->FindObject(search_object_P1.c_str());
                    tmpDataHist1D_P2 = (TH1D*)allDataSamples1D->FindObject(search_object_P2.c_str());
                }
                else if(g_mode_fake_data == true)
                {
                    tmpDataHist1D_P1 = (TH1D*)allFakeDataSamples1D->FindObject(search_object_P1.c_str());
                    tmpDataHist1D_P2 = (TH1D*)allFakeDataSamples1D->FindObject(search_object_P2.c_str());
                }

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

                tmpData1D_P1 = new Double_t[tmpDataHist1D_P1->GetNbinsX()]; //tmpHistP1->Clone("tmpData1D_P1");
                tmpData1D_P2 = new Double_t[tmpDataHist1D_P2->GetNbinsX()]; //tmpHistP2->Clone("tmpData1D_P2");
                for(Int_t bin_x{1}; bin_x <= tmpDataHist1D_P1->GetNbinsX(); ++ bin_x)
                {
                    tmpData1D_P1[bin_x] = tmpDataHist1D_P1->GetBinContent(bin_x);
                }
                for(Int_t bin_x{1}; bin_x <= tmpDataHist1D_P2->GetNbinsX(); ++ bin_x)
                {
                    tmpData1D_P2[bin_x] = tmpDataHist1D_P2->GetBinContent(bin_x);
                }

                // TODO: problem with this method is that it continually creates
                // new histograms due to clone, which is very slow?
                // might be better to organize data into vectors, and use those
                // to store the bin contents
                // NOTE: DONE

                //TH1D *tmpTotalMC1D_P1 = nullptr;
                //TH1D *tmpTotalMC1D_P2 = nullptr;
                Double_t *tmpTotalMC1D_P1 = new Double_t[tmpDataHist1D_P1->GetNbinsX()];
                Double_t *tmpTotalMC1D_P2 = new Double_t[tmpDataHist1D_P2->GetNbinsX()];
                for(Int_t bin_x{1}; bin_x <= tmpDataHist1D_P1->GetNbinsX(); ++ bin_x)
                {
                    tmpTotalMC1D_P1[bin_x] = 0.0;
                }
                for(Int_t bin_x{1}; bin_x <= tmpDataHist1D_P2->GetNbinsX(); ++ bin_x)
                {
                    tmpTotalMC1D_P2[bin_x] = 0.0;
                }

                // loop over all the parameters
                std::map<int, file_parameter>::iterator it{g_pg.file_params.begin()};
                for(; it != g_pg.file_params.end(); ++ it)
                {
                    int paramNumberInt = -1;

                    int paramNumber = it->second.paramNumber;
                    bool paramEnabled = it->second.paramEnabled;
                    bool paramEnabledP1 = it->second.paramEnabledP1;
                    bool paramEnabledP2 = it->second.paramEnabledP2;
                    double paramInitValue = it->second.paramInitValue;
                    double paramInitError = it->second.paramInitError;
                    int paramConstraintMode = it->second.paramConstraintMode;

                    if(debuglevel >= 5)
                    {
                        std::cout << "paramNumber=" << paramNumber << std::endl;
                    }

                    bool ok = false;
                    if(paramEnabled == true)
                    {
                        if(gEnablePhase1 == true)
                        {
                            if(paramEnabledP1 == true)
                            {
                                ok = true;
                            }
                        }

                        if(gEnablePhase2 == true)
                        {
                            if(paramEnabledP2 == true)
                            {
                                ok = true;
                            }
                        }
                    }
                    else
                    {
                        ok = true;
                        continue;
                    }
                    if(ok == false)
                    {
                        if(debuglevel > 0)
                        {
                            std::cout << __func__ << " ok == false" << std::endl;
                            std::cin.get();
                        }
                        continue;
                    }
                    // TODO: this correctly ignores any parameter which is disabled such that
                    // paramEnabled == false
                    // however, it also ignores parameters which are set as disabled for P1
                    // and P2, and when these are irrelevent due to the value of gEnablePhaseX
                    // so... the internal and external index will not match
                    // need to add some code to fix this when the parameters are read from
                    // file, (probably)
                    // unless I just ignore that here... perhaps paramEnabled dictates
                    // whether parameter is drawn and the phase1/phase2 enable flag
                    // is to decide whether minuit does the fit or not (in which case
                    // the param may still contribute to chisquare but may not be minimized
                    // by minuit)


                    // looping over parameters...
                    // get a paramNumber
                    // corresponds to a paramNumberInternal (minuit param number)
                    // which is the index of the vector param where the current
                    // activity can be found
                    // also get a paramName, which is something like ac228_bi212_tl208_int
                    // which is not very useful
                    // but get a list of mc names, which for this example will be
                    // ac228_int_rot, bi212_int_rot, tl208_int_rot
                    // the corresponding MC samples names will be
                    // theHistogram + BkgFiles[i] + "_Px_fit"
                    // aka:
                    // hTotalE_ + bi212_int_rot + _P1_fit
                    // hTotalE_ + bi212_int_rot + _P2_fit
                    // now, the paramNumber can be used to look up the parameter in
                    // g_pg, from which we can obtain paramEnabled, paramEnabledP1
                    // and paramEnabledP2, as well as the constraint mode, constraint
                    // value and error
                    // so, have to strip off the "_fit" ending, and strip off the
                    // histogram (channel) name from the front, aka remove "hTotalE_"
                    // then left with something like
                    // bi212_int_rot_P1 and a channel number which comes from the loop
                    // over allDataSamples1D->GetEntries()
                    // EXCEPT this will now be wrong, because we have data for P1 and P2
                    // so what we want to loop over is an index channel which runs from
                    // 0 to number1DHists
                    // so we have a channel index which is used to select the channel
                    // which selects the index of allMCSamples1D[channel]
                    // we then have to somehow "get" from this TObject
                    // the histogram with the name bi212_int_rot_P1
                    // and bi212_int_rot_P2
                    // once we have done that we switch on paramEnabled, and paramEnabledP1
                    // and paramEnabledP2 to decide whether or not we skip (in the case
                    // where paramEnabled == false), apply a scaling with param[index]
                    // for P1 and P2 if they are enabled, or simply add a value
                    // where the scaling factor is given by paramInitValue for P1 and P2
                    // if they are disabled
                    // TODO: does this last part make sense? no it doesn't. skip if
                    // P1/P2 not enabled

                    std::vector<std::string>::iterator mc_name_it{it->second.MCNameList.begin()};
                    for(; mc_name_it != it->second.MCNameList.end(); ++ mc_name_it)
                    {
                        std::string mc_name = *mc_name_it;
                        std::string histname = std::string(channel_histname_1D[channel]);
                        std::string search_object_P1 = histname + mc_name + "_P1_fit";
                        std::string search_object_P2 = histname + mc_name + "_P2_fit";
                        TH1D *tmpHist1D_P1 = nullptr;
                        TH1D *tmpHist1D_P2 = nullptr;
        
                        if(debuglevel >= 6)
                        {
                            std::cout << "search_object_P1=" << search_object_P1
                                      << " search_object_P2=" << search_object_P2 << std::endl;
                        }

                        paramNumberInt = g_pg.ExtToIntParamNumberMap.at(paramNumber);
                        if(debuglevel >= 5)
                        {
                            std::cout << "paramNumber=" << paramNumber << " -> " << paramNumberInt << std::endl;
                        }

                        if(paramEnabledP1 == true)
                        {

                            tmpHist1D_P1 = (TH1D*)allMCSamples1D[channel]->FindObject(search_object_P1.c_str());

                            if(tmpHist1D_P1 == nullptr)
                            {
                                std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
                                throw "problem";
                            }

                            Double_t scale_factor_P1 = param.at(paramNumberInt);
                            if(debuglevel >= 6)
                            {
                                std::cout << "enabled P1: scale factor P1: " << scale_factor_P1 << std::endl;
                            }
                            //tmpHist1D_P1->Scale(scale_factor_P1);

                            for(Int_t bin_x{1}; bin_x <= tmpHist1D_P1->GetNbinsX(); ++ bin_x)
                            {
                                tmpTotalMC1D_P1[bin_x] += scale_factor_P1 * tmpHist1D_P1->GetBinContent(bin_x);
                            }

                            // found the mc samples, they are stored in tmpHist1D_P1 and tmpHist1D_P2
                            /*
                            if(tmpTotalMC1D_P1 == nullptr)
                            {
                                tmpTotalMC1D_P1 = tmpHist1D_P1->Clone("tmpTotalMC1D_P1");
                            }
                            else
                            {
                                tmpTotalMC1D_P1->Add(tmpHist1D_P1);
                            }
                            */
                        }
                        else
                        {
                            if(debuglevel >= 4)
                            {
                                std::cout << "disabled P1" << std::endl;
                            }
                        }

                        if(paramEnabledP2 == true)
                        {

                            tmpHist1D_P2 = (TH1D*)allMCSamples1D[channel]->FindObject(search_object_P2.c_str());
                        
                            if(tmpHist1D_P2 == nullptr)
                            {
                                std::cout << "ERROR: Could not find object " << search_object_P2 << std::endl;
                                throw "problem";
                            }

                            Double_t scale_factor_P2 = param.at(paramNumberInt);
                            if(debuglevel >= 4)
                            {
                                std::cout << "enabled P2: scale factor P2: " << scale_factor_P2 << std::endl;
                            }
                            //tmpHist1D_P2->Scale(scale_factor_P2);

                            for(Int_t bin_x{1}; bin_x <= tmpHist1D_P2->GetNbinsX(); ++ bin_x)
                            {
                                tmpTotalMC1D_P2[bin_x] += scale_factor_P2 * tmpHist1D_P2->GetBinContent(bin_x);
                            }
                            
                            // found the mc samples, they are stored in tmpHist1D_P1 and tmpHist1D_P2
                            /*
                            if(tmpTotalMC1D_P2 == nullptr)
                            {
                                tmpTotalMC1D_P2 = tmpHist1D_P1->Clone("tmpTotalMC1D_P2");
                            }
                            else
                            {
                                tmpTotalMC1D_P2->Add(tmpHist1D_P2);
                            }
                            */
                        }
                        else
                        {
                            if(debuglevel >= 4)
                            {
                                std::cout << "disabled P2" << std::endl;
                            }
                        }




                        // TODO: what is parameter is NOT enabled

                    } // mc sample name iterator

                } // file_param iterator


                ///////////////////////////////////////////////////////////////////
                // loop over all bins and calculate LL (Phase 1 & 2)
                ///////////////////////////////////////////////////////////////////

                // NOTE: this block assumes that tmpData1D_P1, tmpData1D_P2 are
                // not nullptr
                // which cannot happen unless all params are disabled for P1 or P2
                // NOTE: assume nbinsx same for both phases

                for(Int_t bin_x = 1; bin_x <= tmpDataHist1D_P1->GetNbinsX(); ++ bin_x)
                {
                
                    ///////////////////////////////////////////////////////////////
                    // PHASE 1
                    ///////////////////////////////////////////////////////////////

                    Double_t nData_P1 = 0.0; // ->GetBinContent(bin_ix);
                    nData_P1 = tmpData1D_P1[bin_x];
                    Double_t nMC_P1 = tmpTotalMC1D_P1[bin_x];

                    // poisson is a probability, so values are between 0 and 1 which means that
                    // log of this value is always negative
                    // NOTE: Log(poisson) is always negative! so fval goes UP NOT DOWN
                    if(nMC_P1 >= 0.0)
                    {
                        if(nMC_P1 == 0.0)
                        {
                            continue;
                        }
                        
                        double lp_P1 = 0.0;
                        if(MODE_METRIC == MODE_LOGPOISSON)
                        {
                            lp_P1 = logpoisson(nData_P1, nMC_P1);
                            //double lp = logpoisson(nFakeData, nMC);
                            //double lp = logpoisson_sterling(nFakeData, nMC);
                        }
                        else if(MODE_METRIC == MODE_CHI2)
                        {
                            // chi2 using sqrt(N), N = number of MC
                            // this should be the same as the poisson method
                            // except for events where nData > 0, nMC <= 0
                            //double chi2_P1 = std::pow(nData_P1 - nMC_P1, 2.0) / std::pow(std::sqrt(nMC_P1), 2.0);

                            //if(nData_P1 == 0.0) continue;
                            //double chi2_P1 = std::pow(nData_P1 - nMC_P1, 2.0) / std::pow(std::sqrt(nData_P1), 2.0);
                            if(nMC_P1 == 0.0) continue;
                            double chi2_P1 = std::pow(nData_P1 - nMC_P1, 2.0) / std::pow(std::sqrt(nMC_P1), 2.0);
                            
                            //double chi2_P1 = std::pow((nData_P1 - nMC_P1) / nData_P1, 2.0);
                            // TODO: presumably I need the if(nMC_P1 == 0.0) rejection
                            lp_P1 = chi2_P1;
                        }

                        ll_channel_P1 += lp_P1;
                        if(nData_P1 != 0.0 && nMC_P1 != 0.0)
                        {
                            ++ ndf_P1;
                        }

                        if(debuglevel >= 6)
                        {
                            std::cout << "bin_x=" << bin_x << " lp_P1=" << lp_P1 << " nData_P1=" << nData_P1 << " nMC_P1=" << nMC_P1 << " (1)" << std::endl;
                        }
                    }
                    // not sure we are dealing with zero bins correctly, should
                    // ignore?
                    // this appears to happen a lot
                    else
                    {
                    //TODO cleanup
                        std::cout << "MC WENT LOWER THAN ZERO - 1D P1" << std::endl;
                        std::cout << "nMC_P1=" << nMC_P1 << std::endl;

                        std::cout << "ignore this bin" << std::endl;
                        continue; // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                        double lp_P1 = 0.0;
                        if(MODE_METRIC == MODE_LOGPOISSON)
                        {
                            lp_P1 = logpoisson(nData_P1, 0.0);
                            //double lp_P1 = logpoisson(nData_P1, nMC_P1);
                            //double lp = logpoisson(nFakeData, nMC);
                            //double lp = logpoisson_sterling(nFakeData, nMC);
                        }
                        else if(MODE_METRIC == MODE_CHI2)
                        {
                            // chi2 using sqrt(N), N = number of MC
                            // this should be the same as the poisson method
                            // except for events where nData > 0, nMC <= 0
                            double chi2_P1 = std::pow(nData_P1 - 0.0, 2.0) / std::pow(std::sqrt(0.0), 2.0);
                            // TODO: presumably I need the if(nMC_P1 == 0.0) rejection
                            lp_P1 = chi2_P1;
                        }
                        //double lp = logpoisson(nFakeData, 1.0e-05);
                        //double lp = logpoisson_sterling(nFakeData, 1.0e-05);
                        ll_channel_P1 += lp_P1;
                        if(nData_P1 != 0.0 && nMC_P1 != 0.0)
                        {
                            ++ ndf_P1;
                        }

                        if(debuglevel >= 3)
                        {
                            std::cout << "bin_x=" << bin_x << " lp_P1=" << lp_P1 << " nData_P1=" << nData_P1 << " nMC_P1=" << nMC_P1 << " (2)" << std::endl;
                        }
                    }


                    ///////////////////////////////////////////////////////////////
                    // PHASE 2
                    ///////////////////////////////////////////////////////////////

                    Double_t nData_P2 = 0.0; // ->GetBinContent(bin_ix);
                    nData_P2 = tmpData1D_P2[bin_x];
                    Double_t nMC_P2 = tmpTotalMC1D_P2[bin_x];

                    // poisson is a probability, so values are between 0 and 1 which means that
                    // log of this value is always negative
                    // NOTE: Log(poisson) is always negative! so fval goes UP NOT DOWN
                    if(nMC_P2 >= 0.0)
                    {
                        if(nMC_P2 == 0.0)
                        {
                            continue;
                        }
                        //double lp_P2 = logpoisson(nData_P2, nMC_P2);
                        double lp_P2 = 0.0;
                        if(MODE_METRIC == MODE_LOGPOISSON)
                        {
                            lp_P2 = logpoisson(nData_P2, nMC_P2);
                            //double lp = logpoisson(nFakeData, nMC);
                            //double lp = logpoisson_sterling(nFakeData, nMC);
                        }
                        else if(MODE_METRIC == MODE_CHI2)
                        {
                            // chi2 using sqrt(N), N = number of MC
                            // this should be the same as the poisson method
                            // except for events where nData > 0, nMC <= 0
                            //double chi2_P2 = std::pow(nData_P2 - nMC_P2, 2.0) / std::pow(std::sqrt(nMC_P2), 2.0);

                            //if(nData_P2 == 0.0) continue;
                            //double chi2_P2 = std::pow(nData_P2 - nMC_P2, 2.0) / std::pow(std::sqrt(nData_P2), 2.0);
                            
                            if(nMC_P2 == 0.0) continue;
                            double chi2_P2 = std::pow(nData_P2 - nMC_P2, 2.0) / std::pow(std::sqrt(nMC_P2), 2.0);

                            if(nData_P2 != nMC_P2)
                            {
                                std::cout << "nData_P2=" << nData_P2 << " nMC_P2=" << nMC_P2 << std::endl;
                            }
                            
                            // TODO: presumably I need the if(nMC_P1 == 0.0) rejection
                            lp_P2 = chi2_P2;
                        }
                        //double lp = logpoisson(nFakeData, nMC);
                        //double lp = logpoisson_sterling(nFakeData, nMC);
                        ll_channel_P2 += lp_P2;
                        if(nData_P1 != 0.0 && nMC_P1 != 0.0)
                        {
                            ++ ndf_P2;
                        }

                        if(debuglevel >= 6)
                        {
                            std::cout << "bin_x=" << bin_x << " lp_P2=" << lp_P2 << " nData_P2=" << nData_P2 << " nMC_P2=" << nMC_P2 << " (1)" << std::endl;
                        }
                    }
                    // not sure we are dealing with zero bins correctly, should
                    // ignore?
                    // this appears to happen a lot
                    else
                    {
                        std::cout << "MC WENT LOWER THAN ZERO - 1D P2" << std::endl;
                        std::cout << "nMC_P2=" << nMC_P2 << std::endl;

                        std::cout << "ignore this bin" << std::endl;
                        continue;

                        double lp_P2 = logpoisson(nData_P2, 0.0);
                        //double lp = logpoisson(nFakeData, 1.0e-05);
                        //double lp = logpoisson_sterling(nFakeData, 1.0e-05);
                        ll_channel_P2 += lp_P2;
                        if(nData_P1 != 0.0 && nMC_P1 != 0.0)
                        {
                            ++ ndf_P2;
                        }

                        if(debuglevel >= 3)
                        {
                            std::cout << "bin_x=" << bin_x << " lp_P2=" << lp_P2 << " nData_P2=" << nData_P2 << " nMC_P2=" << nMC_P2 << " (2)" << std::endl;
                        }
                    }

                } // bin_ix


                delete [] tmpData1D_P1;
                delete [] tmpData1D_P2;
                tmpData1D_P1 = nullptr;
                tmpData1D_P2 = nullptr;

                delete [] tmpTotalMC1D_P1;
                delete [] tmpTotalMC1D_P2;
                tmpTotalMC1D_P1 = nullptr;
                tmpTotalMC1D_P2 = nullptr;


                if(debuglevel >= 3)
                {
                    std::cout << "1D: channel " << channel << " enabled, ll_chanel_P1=" << ll_channel_P1 << " ll_channel_P2=" << ll_channel_P2 << std::endl;
                }
                loglik += ll_channel_P1;
                loglik += ll_channel_P2;

            } // channel

            ndf = ndf_P1 + ndf_P2;
            if(debuglevel >= 4)
            {
                std::cout << "ndf=" << ndf << " ndf_P1=" << ndf_P2 << " ndf_P2=" << ndf_P2 << std::endl;
            }
            ndf_P1 = 0;
            ndf_P2 = 0;

            
            ///////////////////////////////////////////////////////////////////////
            // loglikelihood, 2D channels
            ///////////////////////////////////////////////////////////////////////

            // loop over all channels
            for(int channel = 0; channel < number2DHists; ++ channel)
            {
                if(debuglevel >= 4)
                {
                    std::cout << "channel=" << channel << " (2D)" << std::endl;
                }

                // check channel enabled
                if(channel_enable_2D[channel] == 0)
                {
                    if(debuglevel >= 5)
                    {
                        std::cout << "2D: channel " << channel << " disabled, skip" << std::endl;
                    }
                    continue;
                }

                // LL value for this channel
                double ll_channel_P1 = 0.0;
                double ll_channel_P2 = 0.0;

                // TODO: what is this for and is it required?
                TString channel_str;
                channel_str.Form("%i", channel);

                Double_t *tmpData2D_P1 = nullptr; //(TH2D*)allDataSamples2D->At(channel);
                Double_t *tmpData2D_P2 = nullptr;
                //Double_t *tmpFakeData2D_P1 = nullptr; //(TH2D*)allFakeDataSamples2D->At(channel);
                //Double_t *tmpFakeData2D_P2 = nullptr;
                std::string histname = std::string(channel_histname_2D[channel]);
                std::string search_object_P1;
                std::string search_object_P2;
                if(g_mode_fake_data == false)
                {
                    search_object_P1 = histname + std::string(DataFile) + "_P1";
                    search_object_P2 = histname + std::string(DataFile) + "_P2";
                }
                else if(g_mode_fake_data == true)
                {
                    search_object_P1 = histname + std::string("fakedata") + "_P1";
                    search_object_P2 = histname + std::string("fakedata") + "_P2";
                }

                if(debuglevel >= 6)
                {
                    std::cout << "search_object_P1=" << search_object_P1
                              << " search_object_P2=" << search_object_P2 << std::endl;
                }

                TH2D *tmpDataHistP1 = nullptr;
                TH2D *tmpDataHistP2 = nullptr;
                
                if(g_mode_fake_data == false)
                {
                    tmpDataHistP1 = (TH2D*)allDataSamples2D->FindObject(search_object_P1.c_str());
                    tmpDataHistP2 = (TH2D*)allDataSamples2D->FindObject(search_object_P2.c_str());
                }
                else if(g_mode_fake_data == true)
                {
                    tmpDataHistP1 = (TH2D*)allFakeDataSamples2D->FindObject(search_object_P1.c_str());
                    tmpDataHistP2 = (TH2D*)allFakeDataSamples2D->FindObject(search_object_P2.c_str());
                }

                tmpData2D_P1 = new Double_t[tmpDataHistP1->GetNbinsX() * tmpDataHistP1->GetNbinsY()];
                tmpData2D_P2 = new Double_t[tmpDataHistP2->GetNbinsX() * tmpDataHistP2->GetNbinsY()];
                for(Int_t bin_x{1}; bin_x <= tmpDataHistP1->GetNbinsX(); ++ bin_x)
                {
                    for(Int_t bin_y{1}; bin_y <= tmpDataHistP1->GetNbinsY(); ++ bin_y)
                    {
                        tmpData2D_P1[bin_x + bin_y * tmpDataHistP1->GetNbinsX()] = tmpDataHistP1->GetBinContent(bin_x, bin_y);
                    }
                }
                for(Int_t bin_x{1}; bin_x <= tmpDataHistP2->GetNbinsX(); ++ bin_x)
                {
                    for(Int_t bin_y{1}; bin_y <= tmpDataHistP2->GetNbinsY(); ++ bin_y)
                    {
                        tmpData2D_P2[bin_x + bin_y * tmpDataHistP2->GetNbinsX()] = tmpDataHistP2->GetBinContent(bin_x, bin_y);
                    }
                }

                // TODO: problem with this method is that it continually creates
                // new histograms due to clone, which is very slow?
                // might be better to organize data into vectors, and use those
                // to store the bin contents

                //TH2D *tmpTotalMC2D_P1 = nullptr;
                //TH2D *tmpTotalMC2D_P2 = nullptr;
                Double_t *tmpTotalMC2D_P1 = new Double_t[tmpDataHistP1->GetNbinsX() * tmpDataHistP1->GetNbinsY()];
                Double_t *tmpTotalMC2D_P2 = new Double_t[tmpDataHistP2->GetNbinsX() * tmpDataHistP2->GetNbinsY()];
                for(Int_t bin_x{1}; bin_x <= tmpDataHistP1->GetNbinsX(); ++ bin_x)
                {
                    for(Int_t bin_y{1}; bin_y <= tmpDataHistP1->GetNbinsY(); ++ bin_y)
                    {
                        tmpTotalMC2D_P1[bin_x + bin_y * tmpDataHistP2->GetNbinsY()] = 0.0;
                    }
                }
                for(Int_t bin_x{1}; bin_x <= tmpDataHistP2->GetNbinsX(); ++ bin_x)
                {
                    for(Int_t bin_y{1}; bin_y <= tmpDataHistP2->GetNbinsY(); ++ bin_y)
                    {
                        tmpTotalMC2D_P2[bin_x + bin_y * tmpDataHistP2->GetNbinsX()] = 0.0;
                    }
                }

                // loop over all the parameters
                std::map<int, file_parameter>::iterator it{g_pg.file_params.begin()};
                for(; it != g_pg.file_params.end(); ++ it)
                {
                    int paramNumberInt = -1;

                    int paramNumber = it->second.paramNumber;
                    bool paramEnabled = it->second.paramEnabled;
                    bool paramEnabledP1 = it->second.paramEnabledP1;
                    bool paramEnabledP2 = it->second.paramEnabledP2;
                    double paramInitValue = it->second.paramInitValue;
                    double paramInitError = it->second.paramInitError;
                    int paramConstraintMode = it->second.paramConstraintMode;

                    if(debuglevel >= 4)
                    {
                        std::cout << "paramNumber=" << paramNumber << std::endl;
                    }

                    bool ok = false;
                    if(paramEnabled == true)
                    {
                        if(gEnablePhase1 == true)
                        {
                            if(paramEnabledP1 == true)
                            {
                                ok = true;
                            }
                        }

                        if(gEnablePhase2 == true)
                        {
                            if(paramEnabledP2 == true)
                            {
                                ok = true;
                            }
                        }
                    }
                    else
                    {
                        ok = true;
                        continue;
                    }
                    if(ok == false)
                    {
                        if(debuglevel >= 0)
                        {
                            std::cout << __func__ << " ok == false" << std::endl;
                            std::cin.get();
                        }
                        continue;
                    }
                    // TODO: this correctly ignores any parameter which is disabled such that
                    // paramEnabled == false
                    // however, it also ignores parameters which are set as disabled for P1
                    // and P2, and when these are irrelevent due to the value of gEnablePhaseX
                    // so... the internal and external index will not match
                    // need to add some code to fix this when the parameters are read from
                    // file, (probably)
                    // unless I just ignore that here... perhaps paramEnabled dictates
                    // whether parameter is drawn and the phase1/phase2 enable flag
                    // is to decide whether minuit does the fit or not (in which case
                    // the param may still contribute to chisquare but may not be minimized
                    // by minuit)

                    std::vector<std::string>::iterator mc_name_it{it->second.MCNameList.begin()};
                    for(; mc_name_it != it->second.MCNameList.end(); ++ mc_name_it)
                    {
                        std::string mc_name = *mc_name_it;
                        std::string histname = std::string(channel_histname_2D[channel]);
                        std::string search_object_P1 = histname + mc_name + "_P1_fit";
                        std::string search_object_P2 = histname + mc_name + "_P2_fit";
                        TH2D *tmpHist2D_P1 = nullptr;
                        TH2D *tmpHist2D_P2 = nullptr;

                        paramNumberInt = g_pg.ExtToIntParamNumberMap.at(paramNumber);
                        if(debuglevel >= 5)
                        {
                            std::cout << "paramNumber=" << paramNumber << " -> " << paramNumberInt << std::endl;
                        }

                        if(paramEnabledP1 == true)
                        {
                            if(debuglevel >= 6)
                            {
                                std::cout << "enabled P1" << std::endl;
                            }

                            tmpHist2D_P1 = (TH2D*)allMCSamples2D[channel]->FindObject(search_object_P1.c_str());

                            if(tmpHist2D_P1 == nullptr)
                            {
                                std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
                                throw "problem";
                            }

                            Double_t scale_factor_P1 = param.at(paramNumberInt);
                            if(debuglevel >= 6)
                            {
                                std::cout << "scale_factor_P1=" << scale_factor_P1 << std::endl;
                            }
                            //tmpHist2D_P1->Scale(scale_factor_P1);

                            for(Int_t bin_x{1}; bin_x <= tmpHist2D_P1->GetNbinsX(); ++ bin_x)
                            {
                                tmpTotalMC2D_P1[bin_x] += scale_factor_P1 * tmpHist2D_P1->GetBinContent(bin_x);
                            }
                        }

                        if(paramEnabledP2 == true)
                        {
                            if(debuglevel >= 6)
                            {
                                std::cout << "enabled P2" << std::endl;
                            }

                            tmpHist2D_P2 = (TH2D*)allMCSamples2D[channel]->FindObject(search_object_P2.c_str());
                        
                            if(tmpHist2D_P2 == nullptr)
                            {
                                std::cout << "ERROR: Could not find object " << search_object_P2 << std::endl;
                                throw "problem";
                            }

                            Double_t scale_factor_P2 = param.at(paramNumberInt);
                            if(debuglevel >= 6)
                            {
                                std::cout << "scale_factor_P2=" << scale_factor_P2 << std::endl;
                            }
                            //tmpHist2D_P2->Scale(scale_factor_P2);

                            for(Int_t bin_x{1}; bin_x <= tmpHist2D_P2->GetNbinsX(); ++ bin_x)
                            {
                                tmpTotalMC2D_P2[bin_x] += scale_factor_P2 * tmpHist2D_P2->GetBinContent(bin_x);
                            }
                        }




                        // TODO: what is parameter is NOT enabled

                    } // mc sample name iterator

                } // file_param iterator


                ///////////////////////////////////////////////////////////////////
                // loop over all bins and calculate LL (Phase 1 & 2)
                ///////////////////////////////////////////////////////////////////

                // NOTE: this block assumes that tmpData1D_P1, tmpData1D_P2 are
                // not nullptr
                // which cannot happen unless all params are disabled for P1 or P2
                // NOTE: assume nbinsx same for both phases

                for(Int_t bin_x = 1; bin_x <= tmpDataHistP1->GetNbinsX(); ++ bin_x)
                {
                    for(Int_t bin_y = 1; bin_y <= tmpDataHistP1->GetNbinsY(); ++ bin_y)
                    {
                        Double_t nData_P1 = 0.0; // ->GetBinContent(bin_ix);
                        Double_t nData_P2 = 0.0; // ->GetBinContent(bin_ix);

                        nData_P1 = tmpData2D_P1[bin_x + bin_y * tmpDataHistP1->GetNbinsX()];
                        nData_P2 = tmpData2D_P2[bin_x + bin_y * tmpDataHistP2->GetNbinsX()];

                        Double_t nMC_P1 = tmpTotalMC2D_P1[bin_x + bin_y * tmpDataHistP1->GetNbinsX()];
                        Double_t nMC_P2 = tmpTotalMC2D_P2[bin_x + bin_y * tmpDataHistP2->GetNbinsX()];

                        
                        ///////////////////////////////////////////////////////////////
                        // PHASE 1
                        ///////////////////////////////////////////////////////////////

                        // poisson is a probability, so values are between 0 and 1 which means that
                        // log of this value is always negative
                        // NOTE: Log(poisson) is always negative! so fval goes UP NOT DOWN
                        if(nMC_P1 >= 0.0)
                        {
                            if(nMC_P1 == 0.0)
                            {
                                continue;
                            }
                            double lp_P1 = 0.0;
                            if(MODE_METRIC == MODE_LOGPOISSON)
                            {
                                lp_P1 = logpoisson(nData_P1, nMC_P1);
                                //double lp = logpoisson(nFakeData, nMC);
                                //double lp = logpoisson_sterling(nFakeData, nMC);
                            }
                            else if(MODE_METRIC == MODE_CHI2)
                            {
                                // chi2 using sqrt(N), N = number of MC
                                // this should be the same as the poisson method
                                // except for events where nData > 0, nMC <= 0
                                double chi2_P1 = std::pow(nData_P1 - nMC_P1, 2.0) / std::pow(std::sqrt(nMC_P1), 2.0);
                                //double chi2_P1 = std::pow(nData_P1 - nMC_P1, 2.0) / std::pow(std::sqrt(nData_P1), 2.0);
                                // TODO: presumably I need the if(nMC_P1 == 0.0) rejection
                                lp_P1 = chi2_P1;
                            }
                            //double lp_P1 = logpoisson(nData_P1, nMC_P1);
                            //double lp = logpoisson(nFakeData, nMC);
                            //double lp = logpoisson_sterling(nFakeData, nMC);
                            ll_channel_P1 += lp_P1;
                            if(nData_P1 != 0.0 && nMC_P1 != 0.0)
                            {
                                ++ ndf_P1;
                            }

                            if(debuglevel >= 6)
                            {
                                std::cout << "bin_x=" << bin_x << " bin_y=" << bin_y << " lp_P1=" << lp_P1 << " nData_P1=" << nData_P1 << " nMC_P1=" << nMC_P1 << " (1)" << std::endl;
                            }
                        }
                        // not sure we are dealing with zero bins correctly, should
                        // ignore?
                        // this appears to happen a lot
                        else
                        {
                            std::cout << "MC WENT LOWER THAN ZERO - 2D P1" << std::endl;
                            std::cout << "nMC_P1=" << nMC_P1 << std::endl;

                            std::cout << "ignore this bin" << std::endl;
                            continue;

                            double lp_P1 = logpoisson(nData_P1, 0.0);
                            //double lp = logpoisson(nFakeData, 1.0e-05);
                            //double lp = logpoisson_sterling(nFakeData, 1.0e-05);
                            ll_channel_P1 += lp_P1;
                            if(nData_P1 != 0.0 && nMC_P1 != 0.0)
                            {
                                ++ ndf_P1;
                            }

                            if(debuglevel >= 3)
                            {
                                std::cout << "bin_x=" << bin_x << " bin_y=" << bin_y << " lp_P1=" << lp_P1 << " nData_P1=" << nData_P1 << " nMC_P1=" << nMC_P1 << " (2)" << std::endl;
                            }
                        }


                        ///////////////////////////////////////////////////////////////
                        // PHASE 2
                        ///////////////////////////////////////////////////////////////

                        // poisson is a probability, so values are between 0 and 1 which means that
                        // log of this value is always negative
                        // NOTE: Log(poisson) is always negative! so fval goes UP NOT DOWN
                        if(nMC_P2 >= 0.0)
                        {
                            if(nMC_P2 == 0.0)
                            {
                                continue;
                            }
                            double lp_P2 = 0.0;
                            if(MODE_METRIC == MODE_LOGPOISSON)
                            {
                                lp_P2 = logpoisson(nData_P2, nMC_P2);
                                //double lp = logpoisson(nFakeData, nMC);
                                //double lp = logpoisson_sterling(nFakeData, nMC);
                            }
                            else if(MODE_METRIC == MODE_CHI2)
                            {
                                // chi2 using sqrt(N), N = number of MC
                                // this should be the same as the poisson method
                                // except for events where nData > 0, nMC <= 0
                                double chi2_P2 = std::pow(nData_P2 - nMC_P2, 2.0) / std::pow(std::sqrt(nMC_P2), 2.0);
                                //double chi2_P2 = std::pow(nData_P2 - nMC_P2, 2.0) / std::pow(std::sqrt(nData_P2), 2.0);
                                // TODO: presumably I need the if(nMC_P2 == 0.0) rejection
                                lp_P2 = chi2_P2;
                            }
                            //double lp_P2 = logpoisson(nData_P2, nMC_P2);
                            //double lp = logpoisson(nFakeData, nMC);
                            //double lp = logpoisson_sterling(nFakeData, nMC);
                            ll_channel_P2 += lp_P2;
                            if(nData_P1 != 0.0 && nMC_P1 != 0.0)
                            {
                                ++ ndf_P2;
                            }

                            if(debuglevel >= 6)
                            {
                                std::cout << "bin_x=" << bin_x << " bin_y=" << bin_y << " lp_P2=" << lp_P2 << " nData_P2=" << nData_P2 << " nMC_P2=" << nMC_P2 << " (1)" << std::endl;
                            }
                        }
                        // not sure we are dealing with zero bins correctly, should
                        // ignore?
                        // this appears to happen a lot
                        else
                        {
                            std::cout << "MC WENT LOWER THAN ZERO - 2D P2" << std::endl;
                            std::cout << "nMC_P2=" << nMC_P2 << std::endl;

                            std::cout << "ignore this bin" << std::endl;
                            continue;

                            double lp_P2 = logpoisson(nData_P2, 0.0);
                            //double lp = logpoisson(nFakeData, 1.0e-05);
                            //double lp = logpoisson_sterling(nFakeData, 1.0e-05);
                            ll_channel_P2 += lp_P2;
                            if(nData_P1 != 0.0 && nMC_P1 != 0.0)
                            {
                                ++ ndf_P2;
                            }

                            if(debuglevel >= 3)
                            {
                                std::cout << "bin_x=" << bin_x << " bin_y=" << bin_y << " lp_P2=" << lp_P2 << " nData_P2=" << nData_P2 << " nMC_P2=" << nMC_P2 << " (2)" << std::endl;
                            }
                        }
                    } // bin_y
                } // bin_x


                delete [] tmpData2D_P1;
                delete [] tmpData2D_P2;
                tmpData2D_P1 = nullptr;
                tmpData2D_P2 = nullptr;

                delete [] tmpTotalMC2D_P1;
                delete [] tmpTotalMC2D_P2;
                tmpTotalMC2D_P1 = nullptr;
                tmpTotalMC2D_P2 = nullptr;


                if(debuglevel >= 5)
                {
                    std::cout << "2D: channel " << channel << " enabled, ll_channel_P1=" << ll_channel_P1 << " ll_channel_P2=" << ll_channel_P2 << std::endl;
                }
                loglik += ll_channel_P1;
                loglik += ll_channel_P2;

            } // channel
            
            ndf += ndf_P1 + ndf_P2;
            if(debuglevel >= 6)
            {
                std::cout << "ndf=" << ndf << std::endl;
            }


            ///////////////////////////////////////////////////////////////////////////
            // loglikelihood, 2D channels
            ///////////////////////////////////////////////////////////////////////////

            

            int mode = MODE_CONSTRAINT_UNDEFINED;

            // penalty terms section
            double penalty_term = 0.0;

            // loop over all the parameters
            std::map<int, file_parameter>::iterator it{g_pg.file_params.begin()};
            for(; it != g_pg.file_params.end(); ++ it)
            {
                int paramNumberInt = -1;

                int paramNumber = it->second.paramNumber;
                bool paramEnabled = it->second.paramEnabled;
                bool paramEnabledP1 = it->second.paramEnabledP1;
                bool paramEnabledP2 = it->second.paramEnabledP2;
                double paramInitValue = it->second.paramInitValue;
                double paramInitError = it->second.paramInitError;
                double paramConstraintValue = it->second.paramConstraintValue;
                double paramConstraintError = it->second.paramConstraintError;
                int paramConstraintMode = it->second.paramConstraintMode;

                // stop it crashing
                // (below)
                if(paramEnabled == false)
                {
                    continue;
                }

                paramNumberInt = g_pg.ExtToIntParamNumberMap.at(paramNumber);

                if(debuglevel >= 5)
                {
                    std::cout << "paramNumber=" << paramNumber << " -> " << paramNumberInt << std::endl;
                }

                if(paramEnabled == true)
                {
                    if((paramEnabledP1 == true) || (paramEnabledP2 == true))
                    {
                        // do nothing
                    }
                    else
                    {
                        // not enabled
                        continue;
                    }
                }
                else
                {
                    // not enabled
                    continue;
                }

                if(paramConstraintMode == MODE_CONSTRAINT_SOFT)
                {
                    if(debuglevel >= 6)
                    {
                        std::cout << "soft" << std::endl;
                    }
                    // do nothing
                }
                else if(paramConstraintMode == MODE_CONSTRAINT_FREE)
                {
                    //-- ndf;
                    if(debuglevel >= 6)
                    {
                        std::cout << "free" << std::endl;
                    }
                    continue;
                }
                else if(paramConstraintMode == MODE_CONSTRAINT_HARD)
                {
                    if(debuglevel >= 6)
                    {
                        std::cout << "hard" << std::endl;
                    }
                    continue;
                }
                else
                {
                    std::cout << "ERROR: Invalid value for paramNumber=" << paramNumber << ", paramConstraintMode=" << paramConstraintMode << std::endl;
                }

                // this parameter is from minuit internal and thus is in minuit
                // internal units (not Bq)
                // have to convert to Bq units
                double param_value = param.at(paramNumberInt); // TODO this might break if we try to fit P1 seperatly and P2 seperatly
                double activity_value_Bq = paramInitValue;
            
                // convert to Bq
                // multiply by the initial value
                // activity_value_Bq should really be called param_initial_value
                //double activity_value_Bq = 0.0;
                //double tmp_err;
                //get_paramInitValueError(thePhase, i, activity_value_Bq, tmp_err);

                double value = param_value * activity_value_Bq;
                //double penalty = std::pow((param_value - constraint) / error, 2.0);
                double penalty = 0.0;

                // TODO
                if(EMODE == 1)
                {
                    // data
                }
                else if(EMODE == 2)
                {
                    // MC
                }   
                else if(EMODE == 3)
                {
                    // quadrature
                }
                
                penalty = std::pow((value - paramConstraintValue) / paramConstraintError, 2.0);
                if(debuglevel >= 5)
                {
                    //std::cout << "j=" << j << std::endl;
                    std::cout << "paramNumber=" << paramNumber
                              << " value=" << value
                              << " param_value=" << param_value
                              << " paramConstraintValue=" << paramConstraintValue
                              << " paramConstraintError=" << paramConstraintError
                              << " penalty=" << penalty
                              << std::endl;
                }


                // TODO: is this the correct error term?
                // error on constraint rather than error on current fit value?

                penalty_term += penalty;
            }

    #if 0
            // TODO: I don't like this - should loop over the enabled params
            // however, this should still work as it is
            // NOTE: this also loops over the xi parameter, but since it is
            // free and not soft, nothing happens
            for(int i = 0; i < g_pg.numberParams(); ++ i)
            {

                //if(std::find(enabled_params.begin(), enabled_params.end(), i) == enabled_params.end())
                //{
                    //std::cout << "parameter number " << param_number << " is disabled" << std::endl;
                    //std::cout << "ERROR: i=" << i << " - parameter is DISABLED" << std::endl;
                    //std::cin.get();
                //    continue;
                //}


                if(g_pg.file_params.at(i).paramEnabled == true)
                {
                    int constraint_mode = g_pg.file_params.at(i).paramConstraintMode;

                    if(constraint_mode == MODE_CONSTRAINT_SOFT)
                    {

                    }
                    else if(constraint_mode == MODE_CONSTRAINT_HARD)
                    {
                        continue;
                    }
                    else if(constraint_mode == MODE_CONSTRAINT_FREE)
                    {
                        continue;
                    }
                    else
                    {
                        std::cout << "ERROR: Invalid value for parameter_number i=" << i << ", constraint_mode=" << constraint_mode << std::endl;
                    }

                }

                // soft constraint is applied here
                double constraint = 0.;
                double error = 0.;
                //get_paramConstraintValueError(thePhase, i, constraint, error);
                constraint = g_pg.file_params.at(i).paramConstraintValue;
                error = g_pg.file_params.at(i).paramConstraintError;
                // NOTE: these values read from parameter list file and thus are in
                // units of activity (Bq)

                // TODO: can optimize this code

                if(error < 0.0)
                {
                    std::cout << "ERROR: Invalid error value: error=" << error << std::endl;
                }
                /*
                if(thePhase == 0)
                {
                    mode = paramConstrainModeP1Map[i];

                    if(mode == MODE_CONSTRAINT_SOFT)
                    {
                        // do nothing, soft constraint will be applied below
                    }
                    else if(mode == MODE_CONSTRAINT_HARD)
                    {
                        // parameter fixed by minuit, continue to next param
                        
                        // NOTE: changed to ignore HARD
                        continue;
                    }
                    else if(mode == MODE_CONSTRAINT_FREE)
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

                    if(mode == MODE_CONSTRAINT_SOFT)
                    {
                        // do nothing, soft constraint will be applied below
                    }
                    else if(mode == MODE_CONSTRAINT_HARD)
                    {
                        // parameter fixed by minuit, continue to next param
                        
                        // NOTE: changed to ignore HARD
                        continue;
                    }
                    else if(mode == MODE_CONSTRAINT_FREE)
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
                */
                
                // soft constraint is applied here
                double constraint = 0.;
                double error = 0.;
                //get_paramConstraintValueError(thePhase, i, constraint, error);
                constraint = g_pg.file_params.at(i).paramConstraintValue;
                error = g_pg.file_params.at(i).paramConstraintError;
                // NOTE: these values read from parameter list file and thus are in
                // units of activity (Bq)

                // TODO: can optimize this code

                if(error < 0.0)
                {
                    std::cout << "ERROR: Invalid error value: error=" << error << std::endl;
                }

                // check if hard parameter and error == 0.0
                /*
                if(mode == MODE_CONSTRAINT_HARD)
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
                */

                //double param_value = p[i];
                // i is index of the parameter number (external / non minuit)
                // convert to internal parameter number (minuit param number)
                /*
                int j = paramNumberToMinuitParamNumberMap.at(i);
                if(j < numberEnabledParams)
                {
                    // ok
                }
                else
                {
                    throw std::runtime_error("error: invalid value of j (internal param number)");
                }
                double param_value = param[j];
                */
                // this parameter is from minuit internal and thus is in minuit
                // internal units (not Bq)
                // have to convert to Bq units
                double param_value = param[i]; // TODO this might break if we try to fit P1 seperatly and P2 seperatly
                double activity_value_Bq = g_pg.file_params.at(i).paramInitValue;
            
                // convert to Bq
                // multiply by the initial value
                // activity_value_Bq should really be called param_initial_value
                //double activity_value_Bq = 0.0;
                //double tmp_err;
                //get_paramInitValueError(thePhase, i, activity_value_Bq, tmp_err);

                double value = param_value * activity_value_Bq;
                //double penalty = std::pow((param_value - constraint) / error, 2.0);
                double penalty = 0.0;

                // TODO
                if(EMODE == 1)
                {
                    // data
                }
                else if(EMODE == 2)
                {
                    // MC
                }   
                else if(EMODE == 3)
                {
                    // quadrature
                }
                
                penalty = std::pow((value - constraint) / error, 2.0);
                if(debuglevel >= 3)
                {
                    //std::cout << "j=" << j << std::endl;
                    std::cout << "i=" << i
                              << " value=" << value
                              << " constraint=" << constraint
                              << " error=" << error << std::endl;
                }


                // TODO: is this the correct error term?
                // error on constraint rather than error on current fit value?

                penalty_term += penalty;
            }
    #endif

            if(debuglevel >= 4)
            {
                std::cout << "penalty_term=" << penalty_term << std::endl;
            }
            /*
            if(penalty_term != 0.0)
            {
                std::cout << "penalty_term=" << penalty_term << std::endl;
            }
            */
          
            if(MODE_METRIC == MODE_LOGPOISSON)
            {
                fval = -2.0 * loglik + penalty_term;
            }
            else if(MODE_METRIC == MODE_CHI2)
            {
                fval = loglik + penalty_term;
            }

            //tmpData->Delete();

            // hook
            //global_chisquare = fval;


            // set last parameter values
            // could also loop over nPar?
            for(int i = 0; i < param.size(); ++ i)
            {
                // TODO: is the index the same?

                double value = param.at(i);
                //paramLastValueMap[i] = value;
                //minuitParamCurrent[i] = value;
                //minuitParamLast[i] = value;
                int paramNumberExt = g_pg.IntToExtParamNumberMap.at(i);
                g_pg.file_params.at(paramNumberExt).paramLastValue = value;
            }

            

            if(debuglevel >= 2)
            {
                std::cout << "fval=" << fval << std::endl;
            }
            //std::cin.get();
            if(debuglevel >= 10)
            {
                std::cout << "waiting" << std::endl;
                std::cin.get();
            }
            return fval;

        } // old method
    }


    protected:

    double error_def;

    public:

    mutable int ndf;

};




///////////////////////////////////////////////////////////////////////////////
// get number MC 1D
///////////////////////////////////////////////////////////////////////////////
#if 0
Double_t getNumberMC1D(
    const Int_t channel,
    const Int_t bin_ix,
    const std::vector<double> &p
    //const Double_t *const p
    )
{


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
    TH1D *tmpHist;
    //int which_param;

    //std::cout << "getting number of MC... "  << channel << std::endl;

    //std::cout << "channel=" << channel << std::endl;
    //std::cout << allMCSamples1D[channel]->GetEntries() << std::endl;


    // copied from above
    //for(int k = 0; k < allMCSamples1D[channel]->GetEntries(); k++)
    for(int j = 0; j < allMCSamples1D[channel]->GetEntries(); j++)
    {

        // disabled params do not appear in this array

        tmpHist = (TH1D*)allMCSamples1D[channel]->At(j);

        TString tmpName = tmpHist->GetName();

        int which_param = -1;
        bool found_param = false;

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
                int param_number = -1;
                bool success = g_pg.convert_MC_name_to_param_number(tmp_sample_name, param_number);
                //if(MCNameToParamNumberMap.count(tmp_sample_name) > 0)
                if(success == true)
                {
                    //int paramNumber = MCNameToParamNumberMap.at(tmp_sample_name);
                    // TODO: removed std::string, change tmpName type to be std::string from TString
                
                    //std::cout << "paramNumber=" << paramNumber << " -> tmp_sample_name=" << tmp_sample_name << " ~> tmpName=" << tmpName << std::endl;                    
                    //which_param = paramNumber;
                    //which_param = paramNumberToMinuitParamNumberMap.at(paramNumber);
                    // TODO: might still need the conversion from param number to minuit param number
                    found_param = true;

                    //std::cout << tmp_sample_name << " -> " << paramNumber << std::endl;

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

            // check here to see if param is disabled
            // TODO
            //
            // since which_param must be an enabled parameter, no longer need this
            // leave for now as a check?
            //if(std::find(enabled_params.begin(), enabled_params.end(), which_param) == enabled_params.end())
            /*
            if(std::find(enabled_params.begin(), enabled_params.end(), minuitParamNumberToParamNumberMap.at(which_param)) == enabled_params.end())
            {
                std::cout << "ERROR: which_param=" << which_param << " - parameter is DISABLED" << std::endl;
                std::cin.get();
            }
            */
            // TODO: enabled checks

            if(which_param == 1)
            {
                std::cout << "Error: which_param == 1 !" << std::endl;
                std::cin.get();
            }

            nMC += p[which_param] * (double)tmpHist->GetBinContent(bin_ix);

        }
        else
        {
            std::cout << "error could not find histogram: tmpName=" << tmpName << std::endl;
        } 

    }

    //std::cin.get();

    return nMC;

}
#endif


#if 0
Double_t getNumberMC2D(
    const Int_t channel,
    const Int_t bin_ix,
    const Int_t bin_iy,
    //const Double_t *const p
    const std::vector<double> &p
    )
{

    double nMC = 0.;

    // (1) grab a hist from the sample list of this channel
    // (2) figure out which parameter it corresponds to
    TH1D *tmpHist;
    //int which_param;

    //std::cout << "getting number of MC... "  << channel << std::endl;
    //std::cout << allMCSamples[channel]->GetEntries() << std::endl;

    // copied from above
    //for(int k = 0; k < allMCSamples1D[channel]->GetEntries(); k++)
    for(int j = 0; j < allMCSamples2D[channel]->GetEntries(); j++)
    {

        // get histogram name
        tmpHist = (TH1D*)allMCSamples2D[channel]->At(j);
        TString tmpName = tmpHist->GetName();

        // set paramter to default value
        int which_param = -1;
        bool found_param = false;

        {
            std::string tmp_hist_name(tmpName);
            auto i_start = tmp_hist_name.find('_') + 1;
            auto i_end = tmp_hist_name.rfind('_');
            if(i_end - i_start > 0)
            {
                // strip off WHAT - TODO comment
                std::string tmp_sample_name = tmp_hist_name.substr(i_start, i_end - i_start);
                //std::cout << "tmp_sample_name=" << tmp_sample_name << std::endl;
                //std::cin.get(); // check still working

                if(MCNameToParamNumberMap.count(tmp_sample_name) > 0)
                {
                    int paramNumber = MCNameToParamNumberMap.at(tmp_sample_name);
                
                    // convert to minuit (internal) parameter number
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

        if(found_param == true)
        {

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

            if(which_param == 1)
            {
                std::cout << "Error: which_param == 1 !" << std::endl;
                std::cin.get();
            }

            // add to number of MC
            nMC += p[which_param] * (double)tmpHist->GetBinContent(bin_ix, bin_iy);
        }
        else
        {
            std::cout << "error could not find histogram: tmpName=" << tmpName << std::endl;
        }

    } // j

    return nMC;


}
#endif





#endif // MINIMIZEFCNAXIALVECTOR_H
