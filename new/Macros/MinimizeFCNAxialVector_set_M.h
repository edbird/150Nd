#ifndef MINIMIZEFCNAXIALVECTOR_SET_M_H
#define MINIMIZEFCNAXIALVECTOR_SET_M_H



void 
MinimizeFCNAxialVector::set_M(const std::vector<double> &param) const
{

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

        // reset M
        //for(Int_t bin_x{1}; bin_x <= M_1D_P1[channel]->GetNbinsX(); ++ bin_x)
        for(Int_t bin_x{0}; bin_x < M_1D_P1_data[channel]->size(); ++ bin_x)
        {
            //M_1D_P1[channel]->SetBinContent(bin_x, 1, 0.0);
            #if VECTOR_RANGE_CHECK
            M_1D_P1_data[channel]->at(bin_x) = 0.0;
            #else
            M_1D_P1_data[channel]->operator[](bin_x) = 0.0;
            #endif
        }
        
        // reset M
        //for(Int_t bin_x{1}; bin_x <= M_1D_P2[channel]->GetNbinsX(); ++ bin_x)
        for(Int_t bin_x{0}; bin_x < M_1D_P2_data[channel]->size(); ++ bin_x)
        {
            //M_1D_P2[channel]->SetBinContent(bin_x, 1, 0.0);
            #if VECTOR_RANGE_CHECK
            M_1D_P2_data[channel]->at(bin_x) = 0.0;
            #else
            M_1D_P2_data[channel]->operator[](bin_x) = 0.0;
            #endif
        }

    // can M change?
    // yes - if xi_31 parameter changes, which is detectable using g_pg
    // as shown above
    // can it change another way?
    // it depends on all "amplitude" parameters, so need to rebuild if any
    // of these parameters changes
    // just assume it always changes for now


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
            //std::cout << "paramNumber=" << paramNumber << std::endl;
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
                    //std::cout << __func__ << " param " << paramNumber <<  " is disabled (overall) skip" << std::endl;
                    //std::cin.get();
                }
                continue;
            }
            if((paramEnabledP1 || paramEnabledP2) == false)
            {
                // if both are false
                if(debuglevel > 0)
                {
                    std::cout << __func__ << " param " << paramNumber << " is disabled for P1 and P2, skip" << std::endl;
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
                        //Double_t content_input = M_1D_P1[channel]->GetBinContent(bin_x + 1, 1);
                        #if VECTOR_RANGE_CHECK
                        Double_t content_input = M_1D_P1_data[channel]->at(bin_x);
                        #else
                        Double_t content_input = M_1D_P1_data[channel]->operator[](bin_x);
                        #endif
                        Double_t content_add = scale_factor_P1 * tmpHist1D_P1->GetBinContent(bin_x + 1);
                        Double_t content_output = content_input + content_add;
                        //M_1D_P1[channel]->SetBinContent(bin_x + 1, 1, content_output);
                        #if VECTOR_RANGE_CHECK
                        M_1D_P1_data[channel]->at(bin_x) = content_output;
                        #else
                        M_1D_P1_data[channel]->operator[](bin_x) = content_output;
                        #endif
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
                        //Int_t super_index = channel * 2 * 50 + 50 + bin_x;
                        //Double_t content_input = M_1D_P2[channel]->GetBinContent(bin_x + 1, 1);
                        #if VECTOR_RANGE_CHECK
                        Double_t content_input = M_1D_P2_data[channel]->at(bin_x);
                        #else
                        Double_t content_input = M_1D_P2_data[channel]->operator[](bin_x);
                        #endif
                        Double_t content_add = scale_factor_P2 * tmpHist1D_P2->GetBinContent(bin_x + 1);
                        Double_t content_output = content_input + content_add;
                        //M_1D_P2[channel]->SetBinContent(bin_x + 1, 1, content_output);
                        #if VECTOR_RANGE_CHECK
                        M_1D_P2_data[channel]->at(bin_x) = content_output;
                        #else
                        M_1D_P2_data[channel]->operator[](bin_x) = content_output;
                        #endif
                        //std::cout << "debug: " << "super_index=" << super_index << " content_input=" << content_input << " content_add=" << content_add << " content_output=" << content_output << " M:" << M->GetBinContent(super_index + 1, 1) << std::endl;

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
    }
}



#endif // MINIMIZEFCNAXIALVECTOR_SET_M_H
