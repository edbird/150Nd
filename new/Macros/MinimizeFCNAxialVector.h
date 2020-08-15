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
        : error_def{0.5} // LL // TODO: check as I am actually minimizing chisquare
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
        if(allFakeDataSamples1D == nullptr)
        {
            build_fake_data();
        }
        */


        // save the walk
        ll_walk.push_back(std::make_pair(param.at(1), param.at(0)));


        bool debugprint = false;
        bool mode_fake_data = g_mode_fake_data; //false; //true;


        // draw the output
        //TString fname;
        //fname.Form("lliter_%d", counter);
        //draw_channel(1, p, std::string(fname));



        // error mode
        // 1 = data
        // 2 = MC
        // 3 = both in quadrature
        const int EMODE = 2;


        int xi_31_ext_param_number = g_pg.get_xi_31_ext_param_number(); 


        if(debugprint)
        {
            std::cout << std::scientific;
            std::cout << "logLikelihood" << std::endl;
            std::cout << "param[0]=" << param[0] << " param[" << xi_31_ext_param_number << "]="
                      << param[xi_31_ext_param_number] << std::endl;
        }


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
            if(debugprint)
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
            if(debugprint)
            {
                std::cout << "dont rebuild 150Nd MC" << std::endl;
            }
        }

        // TODO: add check here to see if any disabled parameters are accessed


        ///////////////////////////////////////////////////////////////////////////
        // loglikelihood, 1D channels
        ///////////////////////////////////////////////////////////////////////////


        double loglik = 0.0;

    //   std::cout << "getting 1D histograms" << std::endl;

        // std::cout << allDataSamples1D->GetEntries()  << std::endl;

        // loop over all channels
        for(int channel = 0; channel < number1DHists; ++ channel)
        {

            // check channel enabled
            if(channel_enable_1D[channel] == 0)
            {
                if(debugprint)
                {
                    std::cout << "1D: channel " << channel << " disabled, skip" << std::endl;
                }
                continue;
            }

            // LL value for this channel
            double ll_channel = 0.0;

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
            std::string search_object_P1 = histname + std::string(DataFile) + "_P1";
            std::string search_object_P2 = histname + std::string(DataFile) + "_P2";
            TH1D *tmpDataHistP1 = nullptr;
            TH1D *tmpDataHistP2 = nullptr;
            
            if(mode_fake_data == false)
            {
                tmpDataHistP1 = (TH1D*)allDataSamples1D->FindObject(search_object_P1.c_str());
                tmpDataHistP2 = (TH1D*)allDataSamples1D->FindObject(search_object_P2.c_str());
            }
            if(mode_fake_data == true)
            {
                tmpDataHistP1 = (TH1D*)allFakeDataSamples1D->FindObject(search_object_P1.c_str());
                tmpDataHistP2 = (TH1D*)allFakeDataSamples1D->FindObject(search_object_P2.c_str());
            }

            tmpData1D_P1 = new Double_t[tmpDataHistP1->GetNbinsX()]; //tmpHistP1->Clone("tmpData1D_P1");
            tmpData1D_P2 = new Double_t[tmpDataHistP2->GetNbinsX()]; //tmpHistP2->Clone("tmpData1D_P2");
            for(Int_t bin_x{1}; bin_x <= tmpDataHistP1->GetNbinsX(); ++ bin_x)
            {
                tmpData1D_P1[bin_x] = tmpDataHistP1->GetBinContent(bin_x);
            }
            for(Int_t bin_x{1}; bin_x <= tmpDataHistP2->GetNbinsX(); ++ bin_x)
            {
                tmpData1D_P2[bin_x] = tmpDataHistP2->GetBinContent(bin_x);
            }

            // TODO: problem with this method is that it continually creates
            // new histograms due to clone, which is very slow?
            // might be better to organize data into vectors, and use those
            // to store the bin contents

            //TH1D *tmpTotalMC1D_P1 = nullptr;
            //TH1D *tmpTotalMC1D_P2 = nullptr;
            Double_t *tmpTotalMC1D_P1 = new Double_t[tmpDataHistP1->GetNbinsX()];
            Double_t *tmpTotalMC1D_P2 = new Double_t[tmpDataHistP2->GetNbinsX()];
            for(Int_t bin_x{1}; bin_x <= tmpDataHistP1->GetNbinsX(); ++ bin_x)
            {
                tmpTotalMC1D_P1[bin_x] = 0.0;
            }
            for(Int_t bin_x{1}; bin_x <= tmpDataHistP2->GetNbinsX(); ++ bin_x)
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

                if(debugprint)
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
                if(ok == false)
                {
                    if(debugprint)
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

                    paramNumberInt = g_pg.ExtToIntParamNumberMap.at(paramNumber);
                    if(debugprint)
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

                    if(paramEnabledP2 == true)
                    {
                        tmpHist1D_P2 = (TH1D*)allMCSamples1D[channel]->FindObject(search_object_P2.c_str());
                    
                        if(tmpHist1D_P2 == nullptr)
                        {
                            std::cout << "ERROR: Could not find object " << search_object_P2 << std::endl;
                            throw "problem";
                        }

                        Double_t scale_factor_P2 = param.at(paramNumberInt);
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
                    double lp_P1 = logpoisson(nData_P1, nMC_P1);
                    //double lp = logpoisson(nFakeData, nMC);
                    //double lp = logpoisson_sterling(nFakeData, nMC);
                    ll_channel += lp_P1;

                    if(debugprint)
                    {
                        std::cout << "bin_x=" << bin_x << " lp_P1=" << lp_P1 << " nData_P1=" << nData_P1 << " nMC_P1=" << nMC_P1 << " (1)" << std::endl;
                    }
                }
                // not sure we are dealing with zero bins correctly, should
                // ignore?
                // this appears to happen a lot
                else
                {
                    std::cout << "MC WENT LOWER THAN ZERO" << std::endl;

                    double lp_P1 = logpoisson(nData_P1, 0.0);
                    //double lp = logpoisson(nFakeData, 1.0e-05);
                    //double lp = logpoisson_sterling(nFakeData, 1.0e-05);
                    ll_channel += lp_P1;

                    if(debugprint)
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
                    double lp_P2 = logpoisson(nData_P2, nMC_P2);
                    //double lp = logpoisson(nFakeData, nMC);
                    //double lp = logpoisson_sterling(nFakeData, nMC);
                    ll_channel += lp_P2;

                    if(debugprint)
                    {
                        std::cout << "bin_x=" << bin_x << " lp_P2=" << lp_P2 << " nData_P2=" << nData_P2 << " nMC_P2=" << nMC_P2 << " (1)" << std::endl;
                    }
                }
                // not sure we are dealing with zero bins correctly, should
                // ignore?
                // this appears to happen a lot
                else
                {
                    std::cout << "MC WENT LOWER THAN ZERO" << std::endl;

                    double lp_P2 = logpoisson(nData_P2, 0.0);
                    //double lp = logpoisson(nFakeData, 1.0e-05);
                    //double lp = logpoisson_sterling(nFakeData, 1.0e-05);
                    ll_channel += lp_P2;

                    if(debugprint)
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


            if(debugprint)
            {
                std::cout << "1D: channel " << channel << " enabled, ll=" << ll_channel << std::endl;
            }
            loglik += ll_channel;

        } // channel


        
        ///////////////////////////////////////////////////////////////////////
        // loglikelihood, 2D channels
        ///////////////////////////////////////////////////////////////////////

        // loop over all channels
        for(int channel = 0; channel < number2DHists; ++ channel)
        {

            // check channel enabled
            if(channel_enable_2D[channel] == 0)
            {
                if(debugprint)
                {
                    std::cout << "2D: channel " << channel << " disabled, skip" << std::endl;
                }
                continue;
            }

            // LL value for this channel
            double ll_channel = 0.0;

            // TODO: what is this for and is it required?
            TString channel_str;
            channel_str.Form("%i", channel);

            Double_t *tmpData2D_P1 = nullptr; //(TH2D*)allDataSamples2D->At(channel);
            Double_t *tmpData2D_P2 = nullptr;
            //Double_t *tmpFakeData2D_P1 = nullptr; //(TH2D*)allFakeDataSamples2D->At(channel);
            //Double_t *tmpFakeData2D_P2 = nullptr;

            std::string histname = std::string(channel_histname_2D[channel]);
            std::string search_object_P1 = histname + std::string(DataFile) + "_P1";
            std::string search_object_P2 = histname + std::string(DataFile) + "_P2";
            TH2D *tmpDataHistP1 = nullptr;
            TH2D *tmpDataHistP2 = nullptr;
            
            if(mode_fake_data == false)
            {
                tmpDataHistP1 = (TH2D*)allDataSamples2D->FindObject(search_object_P1.c_str());
                tmpDataHistP2 = (TH2D*)allDataSamples2D->FindObject(search_object_P2.c_str());
            }
            if(mode_fake_data == true)
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

                if(debugprint)
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
                if(ok == false)
                {
                    if(debugprint)
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
                    if(debugprint)
                    {
                        std::cout << "paramNumber=" << paramNumber << " -> " << paramNumberInt << std::endl;
                    }

                    if(paramEnabledP1 == true)
                    {
                        tmpHist2D_P1 = (TH2D*)allMCSamples2D[channel]->FindObject(search_object_P1.c_str());

                        if(tmpHist2D_P1 == nullptr)
                        {
                            std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
                            throw "problem";
                        }

                        Double_t scale_factor_P1 = param.at(paramNumberInt);
                        //tmpHist2D_P1->Scale(scale_factor_P1);

                        for(Int_t bin_x{1}; bin_x <= tmpHist2D_P1->GetNbinsX(); ++ bin_x)
                        {
                            tmpTotalMC2D_P1[bin_x] += scale_factor_P1 * tmpHist2D_P1->GetBinContent(bin_x);
                        }
                    }

                    if(paramEnabledP2 == true)
                    {
                        tmpHist2D_P2 = (TH2D*)allMCSamples2D[channel]->FindObject(search_object_P2.c_str());
                    
                        if(tmpHist2D_P2 == nullptr)
                        {
                            std::cout << "ERROR: Could not find object " << search_object_P2 << std::endl;
                            throw "problem";
                        }

                        Double_t scale_factor_P2 = param.at(paramNumberInt);
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
                        double lp_P1 = logpoisson(nData_P1, nMC_P1);
                        //double lp = logpoisson(nFakeData, nMC);
                        //double lp = logpoisson_sterling(nFakeData, nMC);
                        ll_channel += lp_P1;

                        if(debugprint)
                        {
                            std::cout << "bin_x=" << bin_x << " bin_y=" << bin_y << " lp_P1=" << lp_P1 << " nData_P1=" << nData_P1 << " nMC_P1=" << nMC_P1 << " (1)" << std::endl;
                        }
                    }
                    // not sure we are dealing with zero bins correctly, should
                    // ignore?
                    // this appears to happen a lot
                    else
                    {
                        std::cout << "MC WENT LOWER THAN ZERO" << std::endl;

                        double lp_P1 = logpoisson(nData_P1, 0.0);
                        //double lp = logpoisson(nFakeData, 1.0e-05);
                        //double lp = logpoisson_sterling(nFakeData, 1.0e-05);
                        ll_channel += lp_P1;

                        if(debugprint)
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
                        double lp_P2 = logpoisson(nData_P2, nMC_P2);
                        //double lp = logpoisson(nFakeData, nMC);
                        //double lp = logpoisson_sterling(nFakeData, nMC);
                        ll_channel += lp_P2;

                        if(debugprint)
                        {
                            std::cout << "bin_x=" << bin_x << " bin_y=" << bin_y << " lp_P2=" << lp_P2 << " nData_P2=" << nData_P2 << " nMC_P2=" << nMC_P2 << " (1)" << std::endl;
                        }
                    }
                    // not sure we are dealing with zero bins correctly, should
                    // ignore?
                    // this appears to happen a lot
                    else
                    {
                        std::cout << "MC WENT LOWER THAN ZERO" << std::endl;

                        double lp_P2 = logpoisson(nData_P2, 0.0);
                        //double lp = logpoisson(nFakeData, 1.0e-05);
                        //double lp = logpoisson_sterling(nFakeData, 1.0e-05);
                        ll_channel += lp_P2;

                        if(debugprint)
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


            if(debugprint)
            {
                std::cout << "2D: channel " << channel << " enabled, ll=" << ll_channel << std::endl;
            }
            loglik += ll_channel;

        } // channel


#if 0
        // there are i samples for each channel
        for(int channel = 0; channel < allDataSamples1D->GetEntries(); ++ channel)
        {
            // currently the P1 and P2 samples are mixed together here
            // TODO: seperate them

            if(channel_enable_1D[channel] == 0)
            {
                if(debugprint)
                {
                    std::cout << "1D: channel " << channel << " disabled, skip" << std::endl;
                }
                continue;
            }
            

            double ll_channel = 0.0;

            TString i_str;
            i_str.Form("%i", channel);
            //std::cout << i << std::endl;

            tmpData1D = (TH1D*)allDataSamples1D->At(channel);
            tmpFakeData1D = (TH1D*)allFakeDataSamples1D->At(channel);

            // std::cout << tmpData1D->Integral() << std::endl;

            int nBinsX = tmpData1D->GetNbinsX();
            for(int bin_ix = 1; bin_ix <= nBinsX; ++ bin_ix)
            {
                //Int_t nData = 0;
                Double_t nData = 0;
                if(mode_fake_data == false)
                {
                    nData = tmpData1D->GetBinContent(bin_ix);
                }
                if(mode_fake_data == true)
                {
                    nData = tmpFakeData1D->GetBinContent(bin_ix);
                }
                // i is the index of the sample ?
                // ix is the bin index
                // p is a pointer to an array of parameter values

                // fake data
                //double nFakeData = getNumberMC1D(channel, bin_ix, paramInitValueMap);
                //for(int k = 0; k < numberEnabledParams; ++ k)
                //{
                //    std::cout << "minuitParamInit[" << k << "]=" << minuitParamInit[k] << " p[" << k << "]=" << p[k] << std::endl;
                //}
                //double nFakeData = getNumberMC1D(channel, bin_ix, minuitParamInit);
                // TODO: double or int
                // TODO: I don't know if this works. param index might be different
                // internal vs external issue
                // the index for paramInitValueMap must be an external param index
                // this is because paramInitValueP1Map and paramInitValueP2Map
                // are indexed using an external param index and therefore
                // paramInitValueMap is also
                // The function getNumberMC1D uses internal index format

                double nMC = getNumberMC1D(channel, bin_ix, param);

                //std::cout << "bin_ix=" << bin_ix << " nFakeData=" << nFakeData << " nMC=" << nMC << std::endl;
                //std::cout << "for bin_ix=" << bin_ix << " nMC=" << nMC << " nData=" << nData << std::endl;
                

                // poisson is a probability, so values are between 0 and 1 which means that
                // log of this value is always negative
                // NOTE: Log(poisson) is always negative! so fval goes UP NOT DOWN
                if(nMC >= 0.0)
                {
                    double lp = logpoisson(nData, nMC);
                    //double lp = logpoisson(nFakeData, nMC);
                    //double lp = logpoisson_sterling(nFakeData, nMC);
                    ll_channel += lp;

                    if(debugprint)
                    {
                        std::cout << "bin_ix=" << bin_ix << " lp=" << lp << " nData=" << nData << " nMC=" << nMC << " (1)" << std::endl;
                    }
                }
                // not sure we are dealing with zero bins correctly, should
                // ignore?
                // this appears to happen a lot
                else
                {
                    std::cout << "MC WENT LOWER THAN ZERO" << std::endl;

                    double lp = logpoisson(nData, 0.0);
                    //double lp = logpoisson(nFakeData, 1.0e-05);
                    //double lp = logpoisson_sterling(nFakeData, 1.0e-05);
                    ll_channel += lp;

                    if(debugprint)
                    {
                        std::cout << "bin_ix=" << bin_ix << " lp=" << lp << " nData=" << nData << " nMC=" << nMC << " (2)" << std::endl;
                    }
                }
            
            } //~bins

            if(debugprint)
            {
                std::cout << "1D: channel " << channel << " enabled, ll=" << ll_channel << std::endl;
            }
            loglik += ll_channel;

        } //channels
#endif        

        // add contribution from g_A related histograms (single electron energy)
        // either as single electron
        // or as high+low energy 1D
        // or as high+low energy 2D


        ///////////////////////////////////////////////////////////////////////////
        // loglikelihood, 2D channels
        ///////////////////////////////////////////////////////////////////////////

        
        //  std::cout << "getting 2D histograms" << std::endl;
#if 0
        TH2D *tmpData2D;
        TH2D *tmpFakeData2D;
        // std::cout << allDataSamples2D->GetEntries()  << std::endl;

        // there are i samples for each channel
        for(int channel = 0; channel < allDataSamples2D->GetEntries(); ++ channel)
        {

            if(channel_enable_2D[channel] == 0)
            {
                if(debugprint)
                {
                    std::cout << "2D: channel " << channel << " disabled, skip" << std::endl;
                }
                continue;
            }

            double ll_channel = 0.0;

            TString i_str;
            i_str.Form("%i", channel);
            //std::cout << i << std::endl;

            // TODO: can I remove this Clone() call safely to improve speed?
            //tmpData1D = (TH1D*)allDataSamples1D->At(i)->Clone("tmpData1D" + i_str + "_");
            tmpData2D = (TH2D*)allDataSamples2D->At(channel);
            tmpFakeData2D = (TH2D*)allFakeDataSamples2D->At(channel);

            // std::cout << tmpData2D->Integral() << std::endl;

            int nBinsX = tmpData2D->GetNbinsX();
            int nBinsY = tmpData2D->GetNbinsY();
            for(int bin_ix = 1; bin_ix <= nBinsX; ++ bin_ix)
            {
                for(int bin_iy = 1; bin_iy <= nBinsY; ++ bin_iy)
                {
                    //Int_t nData = (Int_t)tmpData2D->GetBinContent(bin_ix, bin_iy);
                    // i is the index of the sample (actually histogram type/channel)
                    // ix is the bin index
                    // p is a pointer to an array of parameter values
                    double nData = 0.0;
                    if(mode_fake_data == false)
                    {
                        nData = tmpData2D->GetBinContent(bin_ix, bin_iy);
                    }
                    if(mode_fake_data == true)
                    {
                        nData = tmpFakeData2D->GetBinContent(bin_ix, bin_iy);
                    }

                    
                    double nMC = getNumberMC2D(channel, bin_ix, bin_iy, param);

                    if(nMC >= 0.0)
                    {
                        double lp = logpoisson(nData, nMC);
                        ll_channel += lp;
                    }
                    else
                    {
                        std::cout << "2D MC WENT LOWER THAN ZERO" << std::endl;

                        double lp = logpoisson(nData, 0.0);
                        ll_channel += lp;
                    }

                } // binsY
            } // binX

            if(debugprint)
            {
                std::cout << "2D: channel " << channel << " enabled, ll=" << ll_channel << std::endl;
            }
            loglik += ll_channel;

        } // channels
#endif

     
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

            paramNumberInt = g_pg.ExtToIntParamNumberMap.at(paramNumber);

            if(debugprint)
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
                // do nothing
            }
            else if(paramConstraintMode == MODE_CONSTRAINT_FREE)
            {
                continue;
            }
            else if(paramConstraintMode == MODE_CONSTRAINT_HARD)
            {
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
            if(debugprint)
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
            if(debugprint)
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

        if(debugprint)
        {
            std::cout << "penalty_term=" << penalty_term << std::endl;
        }
        /*
        if(penalty_term != 0.0)
        {
            std::cout << "penalty_term=" << penalty_term << std::endl;
        }
        */
      
        fval = -2.0 * loglik + penalty_term;

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

        return fval;
    }



    protected:

    double error_def;


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
