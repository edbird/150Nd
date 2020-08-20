#ifndef NEWLOGLIKFITTER_BUILDFAKEDATA_H
#define NEWLOGLIKFITTER_BUILDFAKEDATA_H


///////////////////////////////////////////////////////////////////////////////
// function declarations
///////////////////////////////////////////////////////////////////////////////



// changing resolution has moved fit point
// does changing x maximum move the contours? (it used to)
// does not appear to do so by changing x upper limit
// or x lower limit
// or y axis upper limit
// disable other backgrounds and check result


void rebuild_150Nd_MC(const double, const double);


// NOTE:
// this function does not work if reweight is called before this function is
// called
// for obvious reasons
void build_fake_data()
{

    //bool debugprint = true;
    int debuglevel = 1;

    allFakeDataSamples1D = new TObjArray();
    allFakeDataSamples2D = new TObjArray();

    for(int channel = 0; channel < number1DHists; ++ channel)
    {
        TH1D *tmpTotalMC1D_P1 = nullptr;
        TH1D *tmpTotalMC1D_P2 = nullptr;

        if(debuglevel >= 2)
        {
            std::cout << "channel=" << channel << std::endl;
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

            if(debuglevel >= 3)
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
                if(debuglevel >= 1)
                {
                    std::cout << __func__ << " ok == false" << std::endl;
                    std::cin.get();
                }
                continue;
            }


            // iterate over list of MC
            std::vector<std::string>::iterator mc_name_it{it->second.MCNameList.begin()};
            for(; mc_name_it != it->second.MCNameList.end(); ++ mc_name_it)
            {
                std::string mc_name = *mc_name_it;
                std::string histname = std::string(channel_histname_1D[channel]);
                std::string search_object_P1 = histname + mc_name + "_P1_fit";
                std::string search_object_P2 = histname + mc_name + "_P2_fit";
                std::string hname_P1 = histname + std::string(DataFile) + "_P1";
                std::string hname_P2 = histname + std::string(DataFile) + "_P2";
                //std::string hname_P1 = histname + "fakedata_2e" + "_P1";
                //std::string hname_P2 = histname + "fakedata_2e" + "_P2";
                TH1D *tmpHist1D_P1 = nullptr;
                TH1D *tmpHist1D_P2 = nullptr;
                std::cout << "hname_P1=" << hname_P1 << std::endl;
                std::cout << "hname_P2=" << hname_P2 << std::endl;
                //std::cin.get();

                if(debuglevel >= 3)
                {
                    std::cout << "search_object_P1=" << search_object_P1
                              << " search_object_P2=" << search_object_P2 << std::endl;
                }

                paramNumberInt = g_pg.ExtToIntParamNumberMap.at(paramNumber);
                if(debuglevel >= 3)
                {
                    std::cout << "paramNumber=" << paramNumber << " -> " << paramNumberInt << std::endl;
                }

                if(paramEnabledP1 == true)
                {

                    if(debuglevel >= 2)
                    {
                        std::cout << "enabled P1" << std::endl;
                    }

                    tmpHist1D_P1 = (TH1D*)allMCSamples1D[channel]->FindObject(search_object_P1.c_str())->Clone();

                    if(tmpHist1D_P1 == nullptr)
                    {
                        std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
                        throw "problem";
                    }

                    Double_t scale_factor_P1 = paramInitValue;
                    if(debuglevel >= 2)
                    {
                        std::cout << "scale factor P1: " << scale_factor_P1 << std::endl;
                    }
                    tmpHist1D_P1->Scale(scale_factor_P1);

                    // found the mc samples, they are stored in tmpHist1D_P1 and tmpHist1D_P2
                    if(tmpTotalMC1D_P1 == nullptr)
                    {
                        if(debuglevel >= 3)
                        {
                            std::cout << "Clone()" << std::endl;
                        }
                        tmpTotalMC1D_P1 = (TH1D*)tmpHist1D_P1->Clone(TString(hname_P1));
                    }
                    else
                    {
                        if(debuglevel >= 3)
                        {
                            std::cout << "Add()" << std::endl;
                        }
                        tmpTotalMC1D_P1->Add(tmpHist1D_P1);
                    }
                }
                else
                {
                    if(debuglevel >= 2)
                    {
                        std::cout << "disabled P1" << std::endl;
                    }
                }

                if(paramEnabledP2 == true)
                {
                    if(debuglevel >= 2)
                    {
                        std::cout << "enabled P2" << std::endl;
                    }

                    tmpHist1D_P2 = (TH1D*)allMCSamples1D[channel]->FindObject(search_object_P2.c_str())->Clone();
                
                    if(tmpHist1D_P2 == nullptr)
                    {
                        std::cout << "ERROR: Could not find object " << search_object_P2 << std::endl;
                        throw "problem";
                    }

                    Double_t scale_factor_P2 = paramInitValue;
                    if(debuglevel >= 2)
                    {
                        std::cout << "scale factor P2: " << scale_factor_P2 << std::endl;
                    }
                    tmpHist1D_P2->Scale(scale_factor_P2);

                    // found the mc samples, they are stored in tmpHist1D_P1 and tmpHist1D_P2
                    if(tmpTotalMC1D_P2 == nullptr)
                    {
                        if(debuglevel >= 3)
                        {
                            std::cout << "Clone()" << std::endl;
                        }
                        tmpTotalMC1D_P2 = (TH1D*)tmpHist1D_P2->Clone(TString(hname_P2));
                    }
                    else
                    {
                        if(debuglevel >= 3)
                        {
                            std::cout << "Add()" << std::endl;
                        }
                        tmpTotalMC1D_P2->Add(tmpHist1D_P2);
                    }
                }
                else
                {
                    if(debuglevel >= 2)
                    {
                        std::cout << "disabled P2" << std::endl;
                    }
                }




                // TODO: what is parameter is NOT enabled

            } // mc sample name iterator

        } // file_param iterator

        allFakeDataSamples1D->Add(tmpTotalMC1D_P1);
        allFakeDataSamples1D->Add(tmpTotalMC1D_P2);

    } // channel


    for(int channel = 0; channel < number2DHists; ++ channel)
    {
        TH2D *tmpTotalMC2D_P1 = nullptr;
        TH2D *tmpTotalMC2D_P2 = nullptr;

        if(debuglevel >= 2)
        {
            std::cout << "channel=" << channel << " (2D)" << std::endl;
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

            if(debuglevel >= 3)
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
                if(debuglevel >= 1)
                {
                    std::cout << __func__ << " ok == false" << std::endl;
                    std::cin.get();
                }
                continue;
            }

            // 
            std::vector<std::string>::iterator mc_name_it{it->second.MCNameList.begin()};
            for(; mc_name_it != it->second.MCNameList.end(); ++ mc_name_it)
            {
                std::string mc_name = *mc_name_it;
                std::string histname = std::string(channel_histname_2D[channel]);
                std::string search_object_P1 = histname + mc_name + "_P1_fit";
                std::string search_object_P2 = histname + mc_name + "_P2_fit";
                TH2D *tmpHist2D_P1 = nullptr;
                TH2D *tmpHist2D_P2 = nullptr;
                std::string hname_P1 = histname + std::string(DataFile) + "_P1";
                std::string hname_P2 = histname + std::string(DataFile) + "_P2";
                //std::string hname_P1 = histname + "fakedata_2e" + "_P1";
                //std::string hname_P2 = histname + "fakedata_2e" + "_P2";
                std::cout << "hname_P1=" << hname_P1 << std::endl;
                std::cout << "hname_P2=" << hname_P2 << std::endl;
                //std::cin.get();

                paramNumberInt = g_pg.ExtToIntParamNumberMap.at(paramNumber);
                if(debuglevel >= 3)
                {
                    std::cout << "paramNumber=" << paramNumber << " -> " << paramNumberInt << std::endl;
                }

                if(paramEnabledP1 == true)
                {
                    if(debuglevel >= 2)
                    {
                        std::cout << "enabled P1" << std::endl;
                    }

                    tmpHist2D_P1 = (TH2D*)allMCSamples2D[channel]->FindObject(search_object_P1.c_str())->Clone();

                    if(tmpHist2D_P1 == nullptr)
                    {
                        std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
                        throw "problem";
                    }

                    Double_t scale_factor_P1 = paramInitValue;
                    if(debuglevel >= 2)
                    {
                        std::cout << "scale_factor_P1=" << scale_factor_P1 << std::endl;
                    }
                    
                    tmpHist2D_P1->Scale(scale_factor_P1);

                    if(tmpTotalMC2D_P1 == nullptr)
                    {
                        if(debuglevel >= 3)
                        {
                            std::cout << "Clone()" << std::endl;
                        }
                        tmpTotalMC2D_P1 = (TH2D*)tmpHist2D_P1->Clone(TString(hname_P1));
                    }
                    else
                    {
                        if(debuglevel >= 3)
                        {
                            std::cout << "Add()" << std::endl;
                        }
                        tmpTotalMC2D_P1->Add(tmpHist2D_P1);
                    }
                }

                if(paramEnabledP2 == true)
                {
                    if(debuglevel >= 2)
                    {
                        std::cout << "enabled P2" << std::endl;
                    }

                    tmpHist2D_P2 = (TH2D*)allMCSamples2D[channel]->FindObject(search_object_P2.c_str())->Clone();
                
                    if(tmpHist2D_P2 == nullptr)
                    {
                        std::cout << "ERROR: Could not find object " << search_object_P2 << std::endl;
                        throw "problem";
                    }

                    Double_t scale_factor_P2 = paramInitValue;
                    if(debuglevel >= 2)
                    {
                        std::cout << "scale_factor_P2=" << scale_factor_P2 << std::endl;
                    }
                    
                    tmpHist2D_P2->Scale(scale_factor_P2);

                    if(tmpTotalMC2D_P2 == nullptr)
                    {
                        if(debuglevel >= 3)
                        {
                            std::cout << "Clone()" << std::endl;
                        }
                        tmpTotalMC2D_P2 = (TH2D*)tmpHist2D_P2->Clone(TString(hname_P2));
                    }
                    else
                    {
                        if(debuglevel >= 3)
                        {
                            std::cout << "Add()" << std::endl;
                        }
                        tmpTotalMC2D_P2->Add(tmpHist2D_P2);
                    }
                }

                // TODO: what is parameter is NOT enabled

            } // mc sample name iterator

        } // file_param iterator

        allFakeDataSamples2D->Add(tmpTotalMC2D_P1);
        allFakeDataSamples2D->Add(tmpTotalMC2D_P2);

    } // channel
































#if 0



    std::cout << "build_fake_data()" << std::endl;


    std::cout << "rebuild: xi_31_baseline=" << xi_31_baseline << std::endl;
    //rebuild_150Nd_MC(/*xi_31_baseline*/ 0.296, xi_31_baseline);
    rebuild_150Nd_MC(xi_31_baseline, xi_31_baseline);


    TH1D *hAllMC1D[number1DHists];
    TH2D *hAllMC2D[number2DHists]; // TODO

    std::cout << "debug: number of data samples: " << allDataSamples1D->GetEntries() << std::endl;
    std::cout << "debug: number of MC samples: " << allMCSamples1D[0]->GetEntries() << std::endl;
        
    allFakeDataSamples1D = new TObjArray();
    allFakeDataSamples2D = new TObjArray();


    // each channel 1D hists
    // this is for(i = 0; i < 1; ++ i)
    // TODO: this isn't right. should this be iterating over the "channel" ?
    for(int i = 0; i < allDataSamples1D->GetEntries(); i++)
    {
//        std::cout << "i=" << i << std::endl;

        // because this isn't right TODO
        // uses at(i), but i should always be zero and there should be an
        // additional array index
        //data1D[i] = (TH1D*)allDataSamples1D->At(i)->Clone();
        

        TH1D *tmpHist;

        // j list MC samples for this channel i
        //std::cout << "debug: number of MC samples (i=" << i << "): " << allMCSamples1D[i]->GetEntries() << std::endl;


        // allMCSamples1D[0] contains objects such as: "zr96_rot_k40_2e_P2"

        // TODO i should be channel here?
        for(int j = 0; j < allMCSamples1D[i]->GetEntries(); j++)
        {

            //std::cout << "j=" << j << std::endl;

            TString j_str;
            j_str.Form("%i", j);

            tmpHist = (TH1D*)allMCSamples1D[i]->At(j)->Clone();
            TString tmpName = tmpHist->GetName();

            //std::cout << "(1) tmpName=" << tmpName << std::endl;

            //std::cout << "looking for " << tmpName << std::endl;
            int which_param = -1;
            bool found_param = false;

            // search through parameters to find right one
            // the histogram names are formatted like:
            // hTotalE_bi214_mylar_fit
            // histogram_name + "_" + mc_sample_name + "_fit"
            
            // used later
            double activity_scale_branching_ratio = 1.0;

            {
                std::string tmp_hist_name(tmpName);
                auto i_start = tmp_hist_name.find('_') + 1;
                auto i_end = tmp_hist_name.rfind('_');
                if(i_end - i_start > 0)
                {
                    std::string tmp_sample_name = tmp_hist_name.substr(i_start, i_end - i_start);

                    // TODO
                    // do not have to scale by 0.36 after reading from file
                    // as scaling is done when reading
                    /*
                    // set branching ratio fraction
                    if(tmp_sample_name == std::string("tl208_int_rot") ||
                       tmp_sample_name == std::string("tl208_feShield") ||
                       tmp_sample_name == std::string("tl208_pmt"))
                    {
                        activity_scale_branching_ratio = 0.36;
                    }
                    */

                    if(MCNameToParamNumberMap.count(tmp_sample_name) > 0)
                    {
                        int paramNumber = MCNameToParamNumberMap.at(tmp_sample_name);
                        // TODO: removed std::string, change tmpName type to be std::string from TString
                    
                        //which_param = paramNumber;
                        which_param = paramNumberToMinuitParamNumberMap.at(paramNumber);
                        found_param = true;
//                        std::cout << "j=" << j << ": paramNumber=" << paramNumber << " -> tmp_sample_name=" << tmp_sample_name << " ~> tmpName=" << tmpName << " which_param=" << which_param << std::endl;
                    }
                    else
                    {
                       std::cout << __func__ << " ERROR: could not find " << tmp_sample_name << " in MCNameToParamNumberMap" << std::endl;
                    }
                }
            }

//            std::cin.get();

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



                // no error thrown, which_param is presumably the correct index
                //Double_t activity_scale = AdjustActs[which_param] * activity_scale_branching_ratio;
                Double_t activity_scale = paramInitValueMap[which_param]; // * activity_scale_branching_ratio;
                if(which_param == 0)
                {
                    std::cout << "in build_fake_data(): activity_scale=" << activity_scale << std::endl;

                    // SSD
                    //activity_scale *= 1.2;
                    //tmpHist->Scale(1.5);
                }
                else if(which_param == 1)
                {
                    std::cout << "ERROR: which_param == 1 !" << std::endl;
                    std::cin.get();
                }
                //std::cout << "activity_scale=" << activity_scale << std::endl;
//                tmpHist->Scale(activity_scale);
// TODO: has already been scaled by this activity when read in

                if(tmpHist->Integral() > 0)
                {

                    TString hname = tmpHist->GetName();

                    if(j == 0)
                    {
                        //std::cout << "Clone() done" << "j=" << j << std::endl;
                        // TODO: bug here if Integral() for j == 0 is zero
                        
                        hAllMC1D[i] = (TH1D*)tmpHist->Clone("Total MC Fake Data");
                        //hAllMC1D[i] = (TH1D*)tmpHist_drawpointer->Clone("Total MC");
                        
                        /*
                        std::cout << "j=" << j << std::endl;
                        for(int k = 0; k < tmpHist->GetNbinsX(); ++ k)
                        {
                            std::cout << "k=" << k << " " << tmpHist->GetBinContent(k) << std::endl;
                        }
                        */
                    }
                    else
                    {
                        hAllMC1D[i]->Add((TH1D*)tmpHist);
                        //hAllMC1D[i]->Add((TH1D*)tmpHist_drawpointer);
                        
                        /*
                        std::cout << "j=" << j << std::endl;
                        for(int k = 0; k < tmpHist->GetNbinsX(); ++ k)
                        {
                            std::cout << "k=" << k << " " << tmpHist->GetBinContent(k) << std::endl;
                        }
                        */
                    }
	            }
                else
                {
                    //std::cout << "not adding to stack, Integral() <= 0: " << tmpHist->GetName() << std::endl;
                }
            }
            else
            {
                std::cout << __func__ << " error could not find histogram: tmpName=" << tmpName << std::endl;
            } 

        }

        /*
        std::cout << "integral for fakedata sample " << i << " is " << hAllMC1D[i]->Integral() << std::endl;
        for(int k = 0; k < hAllMC1D[i]->GetNbinsX(); ++ k)
        {
            std::cout << "k=" << k << " " << hAllMC1D[i]->GetBinContent(k) << std::endl;
        }
        */
        allFakeDataSamples1D->Add((TH1D*)hAllMC1D[i]);


        std::cout << "The integral for fake data is " << hAllMC1D[i]->Integral() << std::endl;
    }







    // each channel 2D hists
    for(int i = 0; i < allDataSamples2D->GetEntries(); i++)
    {

        TH2D *tmpHist;

        // allMCSamples2D[0] contains objects such as: "zr96_rot_k40_2e_P2"

        for(int channel = 0; channel < allMCSamples2D[i]->GetEntries(); ++ channel)
        {

            TString channel_str;
            channel_str.Form("%i", channel);

            tmpHist = (TH2D*)allMCSamples2D[i]->At(channel)->Clone();
            TString tmpName = tmpHist->GetName();

            //std::cout << "looking for " << tmpName << std::endl;
            int which_param = -1;
            bool found_param = false;

            // search through parameters to find right one
            // the histogram names are formatted like:
            // hTotalE_bi214_mylar_fit
            // histogram_name + "_" + mc_sample_name + "_fit"
            
            // used later
            double activity_scale_branching_ratio = 1.0;

            {
                std::string tmp_hist_name(tmpName);
                auto i_start = tmp_hist_name.find('_') + 1;
                auto i_end = tmp_hist_name.rfind('_');
                if(i_end - i_start > 0)
                {
                    std::string tmp_sample_name = tmp_hist_name.substr(i_start, i_end - i_start);

                    // TODO
                    // do not have to scale by 0.36 after reading from file
                    // as scaling is done when reading
                    /*
                    // set branching ratio fraction
                    if(tmp_sample_name == std::string("tl208_int_rot") ||
                       tmp_sample_name == std::string("tl208_feShield") ||
                       tmp_sample_name == std::string("tl208_pmt"))
                    {
                        activity_scale_branching_ratio = 0.36;
                    }
                    */

                    if(MCNameToParamNumberMap.count(tmp_sample_name) > 0)
                    {
                        int paramNumber = MCNameToParamNumberMap.at(tmp_sample_name);
                        // TODO: removed std::string, change tmpName type to be std::string from TString
                    
                        //which_param = paramNumber;
                        which_param = paramNumberToMinuitParamNumberMap.at(paramNumber);
                        found_param = true;
//                        std::cout << "j=" << j << ": paramNumber=" << paramNumber << " -> tmp_sample_name=" << tmp_sample_name << " ~> tmpName=" << tmpName << " which_param=" << which_param << std::endl;
                    }
                    else
                    {
                       std::cout << __func__ << " ERROR: could not find " << tmp_sample_name << " in MCNameToParamNumberMap" << std::endl;
                    }
                }
            }

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

                // no error thrown, which_param is presumably the correct index
                //Double_t activity_scale = AdjustActs[which_param] * activity_scale_branching_ratio;
                Double_t activity_scale = paramInitValueMap[which_param]; // * activity_scale_branching_ratio;
                if(which_param == 0)
                {
                    std::cout << "in build_fake_data(): activity_scale=" << activity_scale << std::endl;

                    // SSD
                    //activity_scale *= 1.2;
                    //tmpHist->Scale(1.5);
                }
                else if(which_param == 1)
                {
                    std::cout << "ERROR: which_param == 1 !" << std::endl;
                    std::cin.get();
                }
                //std::cout << "activity_scale=" << activity_scale << std::endl;
//                tmpHist->Scale(activity_scale);
// TODO: has already been scaled by this activity when read in

                if(tmpHist->Integral() > 0)
                {

                    TString hname = tmpHist->GetName();

                    if(channel == 0)
                    {
                        //std::cout << "Clone() done" << "j=" << j << std::endl;
                        // TODO: bug here if Integral() for j == 0 is zero
                        
                        hAllMC2D[i] = (TH2D*)tmpHist->Clone("Total MC Fake Data");
                    }
                    else
                    {
                        hAllMC2D[i]->Add((TH2D*)tmpHist);
                    }
	            }
                else
                {
                    //std::cout << "not adding to stack, Integral() <= 0: " << tmpHist->GetName() << std::endl;
                }
            }
            else
            {
                std::cout << __func__ << " error could not find histogram: tmpName=" << tmpName << std::endl;
            } 

        }

        allFakeDataSamples2D->Add((TH2D*)hAllMC2D[i]);


        std::cout << "The integral for fake data is " << hAllMC2D[i]->Integral() << std::endl;
    }

#endif




}


















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




void rebuild_150Nd_MC(const double xi_31, const double xi_31_baseline)
{

    int debuglevel = 1;

    ///////////////////////////////////////////////////////////////////////
    // reweight hMinMaxEnergy_
    // reweight all
    ///////////////////////////////////////////////////////////////////////


    // new code to reweight 150Nd by xi_{31} parameter

    // pointers of histograms to pass to reweight function
    TH1D *hTotalE_P1 = nullptr;
    TH1D *hSingleEnergy_P1 = nullptr;
    TH1D *hHighEnergy_P1 = nullptr;
    TH1D *hLowEnergy_P1 = nullptr;
    TH1D *hEnergySum_P1 = nullptr;
    TH1D *hEnergyDiff_P1 = nullptr;
    TH2D *hHighLowEnergy_P1 = nullptr;

    TH1D *hTotalE_P2 = nullptr;
    TH1D *hSingleEnergy_P2 = nullptr;
    TH1D *hHighEnergy_P2 = nullptr;
    TH1D *hLowEnergy_P2 = nullptr;
    TH1D *hEnergySum_P2 = nullptr;
    TH1D *hEnergyDiff_P2 = nullptr;
    TH2D *hHighLowEnergy_P2 = nullptr;

    /*
    TH1D *hSingleEnergyClone = nullptr;
    */

    // search through each channel for 150nd samples
    //for(int channel = 0; channel < allDataSamples1D->GetEntries(); ++ channel)
    for(int channel = 0; channel < number1DHists; ++ channel)
    {
        if(debuglevel >= 3)
        {
            std::cout << "channel=" << channel << std::endl;
        }
 
        std::string histname = std::string(channel_histname_1D[channel]);
        std::string search_object_P1 = histname + std::string(Nd150Files[0]) + "_P1_fit";
        std::string search_object_P2 = histname + std::string(Nd150Files[0]) + "_P2_fit";

        if(debuglevel >= 3)
        {
            std::cout << "search_object_P1=" << search_object_P1 << std::endl;
            std::cout << "search_object_P2=" << search_object_P2 << std::endl;
        }

        TObject *tmpobj_P1 = allMCSamples1D[channel]->FindObject(search_object_P1.c_str());
        TObject *tmpobj_P2 = allMCSamples1D[channel]->FindObject(search_object_P2.c_str());
        if(tmpobj_P1 == nullptr)
        {
            std::cout << "could not find object: " << search_object_P1 << std::endl;
            /*for(int i = 0; i < allMCSamples1D[channel]->GetEntries(); ++ i)
            {
                std::cout << allMCSamples1D[channel]->At(i)->GetName() << std::endl;
            }*/
            throw "problem";
        }
        if(tmpobj_P2 == nullptr)
        {
            std::cout << "could not find object: " << search_object_P2 << std::endl;
            throw "problem";
        }

        if(debuglevel >= 5)
        {
            std::cout << "before ->Remove():" << std::endl;
            for(int i = 0; i < allMCSamples1D[channel]->GetEntries(); ++ i)
            {
                std::cout << allMCSamples1D[channel]->At(i)->GetName() << std::endl;
            }
        }
        // NOTE: this just removes the object from the TObjArray
        // it does not delete the object from gROOT memory
        if(allMCSamples1D[channel]->Remove(tmpobj_P1) == nullptr)
        {
            std::cout << "FAILED TO REMOVE" << std::endl;
            throw "failed to remove";
        }
        if(allMCSamples1D[channel]->Remove(tmpobj_P2) == nullptr)
        {
            std::cout << "FAILED TO REMOVE" << std::endl;
            throw "failed to remove";
        }
        if(debuglevel >= 5)
        {
            std::cout << "after ->Remove():" << std::endl;
            for(int i = 0; i < allMCSamples1D[channel]->GetEntries(); ++ i)
            {
                std::cout << allMCSamples1D[channel]->At(i)->GetName() << std::endl;
            }
        }
    } // channel

    // search through each channel for 150nd samples
    //for(int channel = 0; channel < allDataSamples2D->GetEntries(); ++ channel)
    for(int channel = 0; channel < number2DHists; ++ channel)
    {
        std::cout << "channel=" << channel << std::endl;

        std::string histname = std::string(channel_histname_2D[channel]);
        std::string search_object_P1 = histname + std::string(Nd150Files[0]) + "_P1_fit";
        std::string search_object_P2 = histname + std::string(Nd150Files[0]) + "_P2_fit";

        if(debuglevel >= 3)
        {
            std::cout << "search_object_P1=" << search_object_P1 << std::endl;
            std::cout << "search_object_P2=" << search_object_P2 << std::endl;
        }

        TObject *tmpobj_P1 = allMCSamples2D[channel]->FindObject(search_object_P1.c_str());
        TObject *tmpobj_P2 = allMCSamples2D[channel]->FindObject(search_object_P2.c_str());
        if(tmpobj_P1 == nullptr)
        {
            std::cout << "could not find object: " << search_object_P1 << std::endl;
            throw "problem";
        }
        if(tmpobj_P2 == nullptr)
        {
            std::cout << "could not find object: " << search_object_P2 << std::endl;
            throw "problem";
        }

        if(debuglevel >= 5)
        {
            std::cout << "before ->Remove():" << std::endl;
            for(int i = 0; i < allMCSamples2D[channel]->GetEntries(); ++ i)
            {
                std::cout << allMCSamples2D[channel]->At(i)->GetName() << std::endl;
            }
        }
        // NOTE: this just removes the object from the TObjArray
        // it does not delete the object from gROOT memory
        if(allMCSamples2D[channel]->Remove(tmpobj_P1) == nullptr)
        {
            std::cout << "FAILED TO REMOVE" << std::endl;
            throw "failed to remove";
        }
        if(allMCSamples2D[channel]->Remove(tmpobj_P2) == nullptr)
        {
            std::cout << "FAILED TO REMOVE" << std::endl;
            throw "failed to remove";
        }
        if(debuglevel >= 5)
        {
            std::cout << "after ->Remove():" << std::endl;
            for(int i = 0; i < allMCSamples2D[channel]->GetEntries(); ++ i)
            {
                std::cout << allMCSamples2D[channel]->At(i)->GetName() << std::endl;
            }
        }
    } // channel

    if(debuglevel >= 2)
    {
        std::cout << "xi_31=" << xi_31 << " xi_31_baseline=" << xi_31_baseline << std::endl;
    }
    //const double xi_31_baseline{0.296};
    // NOTE: 2020-06-17 this was a bug, removed

    TH1D *hWeight_P1 = nullptr;
    TH1D *hWeight_P2 = nullptr; // TODO: move?
    if(debuglevel >= 5)
    {
        std::cout << "before reweight_apply()" << std::endl;
    }
    reweight_apply(
        hTotalE_P1,
        hSingleEnergy_P1,
        hHighEnergy_P1,
        hLowEnergy_P1,
        hEnergySum_P1,
        hEnergyDiff_P1,
        hHighLowEnergy_P1,
        hWeight_P1,
        hTotalE_P2,
        hSingleEnergy_P2,
        hHighEnergy_P2,
        hLowEnergy_P2,
        hEnergySum_P2,
        hEnergyDiff_P2,
        hHighLowEnergy_P2,
        hWeight_P2,
        xi_31,
        xi_31_baseline,
        h_nEqNull,
        h_nEqTwo,
        psiN0,
        psiN2,
        bb_Q);
    if(debuglevel >= 5)
    {
        std::cout << "after reweight_apply()" << std::endl;
    }

    /*
    hSingleEnergyClone->Divide(hSingleEnergy);
    for(Int_t ii = 1; ii <= hSingleEnergyClone->GetNbinsX(); ++ ii)
    {
        float content = hSingleEnergyClone->GetBinContent(ii);
        std::cout << "content: " << ii << " " << content << std::endl;
    }
    std::cin.get();
    */


    // TODO: just another example of manual code edits
    // make a file describing the channels to fit as well as the parameters
    if(debuglevel >= 4)
    {
        std::cout << "adding P1 histograms" << std::endl;
    }
    allMCSamples1D[0]->Add(hTotalE_P1);
    allMCSamples1D[1]->Add(hSingleEnergy_P1);
    allMCSamples1D[2]->Add(hHighEnergy_P1);
    allMCSamples1D[3]->Add(hLowEnergy_P1);
    allMCSamples1D[4]->Add(hEnergySum_P1);
    allMCSamples1D[5]->Add(hEnergyDiff_P1);

    allMCSamples2D[0]->Add(hHighLowEnergy_P1);
    
    if(debuglevel >= 4)
    {
        std::cout << "adding P2 histograms" << std::endl;
    }
    allMCSamples1D[0]->Add(hTotalE_P2);
    allMCSamples1D[1]->Add(hSingleEnergy_P2);
    allMCSamples1D[2]->Add(hHighEnergy_P2);
    allMCSamples1D[3]->Add(hLowEnergy_P2);
    allMCSamples1D[4]->Add(hEnergySum_P2);
    allMCSamples1D[5]->Add(hEnergyDiff_P2);


    allMCSamples2D[0]->Add(hHighLowEnergy_P2);
}











#endif // NEWLOGLIKFITTER_BUILDFAKEDATA_H
