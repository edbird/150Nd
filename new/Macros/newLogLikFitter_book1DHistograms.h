#ifndef NEWLOGLIKFITTER_BOOK1DHISTOGRAMS_H
#define NEWLOGLIKFITTER_BOOK1DHISTOGRAMS_H




///////////////////////////////////////////////////////////////////////////////
// book1DHistograms
///////////////////////////////////////////////////////////////////////////////

#if 0
void book1DHistograms_helper(TFile *myFile, Int_t channel_counter, TString theChannel, TString thePhase_arg, TString theHistogram,
    const int nBkgs, TString *BkgFiles)//, TH1D *tmpHist)
{
        
    TH1D *tmpHist = nullptr;

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

                //tmpHist = (TH1D*)gDirectory->Get(fullname.c_str())->Clone();
                tmpHist = (TH1D*)myFile->Get(fullname.c_str())->Clone(new_name.c_str());

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
                    //tmpHist = (TH1D*)gDirectory->Get(name.c_str())->Clone(hist_name.c_str());
                    //tmpHist = (TH1D*)gDirectory->Get(fullname.c_str())->Clone();

                    // scale by activity

                    // convert parameter number to minuit parameter number
                    //minuit_param_number = paramNumberToMinuitParamNumberMap.at(param_number);

                    // TODO: change such that samples are pre-scaled by activity input value
                    // get initial parameter values and error
                    Double_t param_init_value = 0.;
                    Double_t param_init_error = 0.; 
                    get_paramInitValueError(thePhase, param_number, param_init_value, param_init_error);

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
                        // TODO: check that this is not already applied in
                        // fit_2e
                        // NOTE: it isn't
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

                    // NOTE: Scale factor
                    // samples scaled by param_init_value (activity in Bq)
                    // after being read from file
                    // however objects in file are scaled by
                    // TotalTime / sampleNGenMC
                    // also scale by 0.36 for relevant samples to account for
                    // branching ratio here
                    
                    // in loglikelihood function the getNumberMC() functions
                    // multiply the bin content by AdjustActs parameter
                    // (minuit parameter)

                    if(param_number == 1)
                    {
                        std::cout << "ERROR" << std::endl;
                        throw "Error";
                    }

                    // NOTE: do NOT apply xi reweighting here
                    // this section just LOADS histograms from file and we want to LOAD
                    // the default (not reweighted) nd150 spectra
                    // TODO: this may no longer be true

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
#endif


void book1DHistograms_helper(
    TFile *aFile_P1,
    TFile *aFile_P2,
    Int_t channel_counter,
    TString theChannel,
    TString theHistogram,
    const int nBkgs,
    TString *BkgFiles)//, TH1D *tmpHist)
{

    int debuglevel = 2;
        
    TH1D *tmpHist_P1 = nullptr;
    TH1D *tmpHist_P2 = nullptr;

    for(int i = 0; i < nBkgs; i++)
    {
        // check if parameter is enabled
        // convert parameter string name to index
        //std::cout << "searching for string " << ExternalBkgOneParamFiles[i] << " in paramNameMap" << std::endl;
        //TString search_object = ExternalBkgOneParamFiles[i];
        std::string mc_name = std::string(BkgFiles[i]);

        // convert mc_name to scale factor
        Double_t scale_factor_P2 = 0.0;
        Double_t scale_factor_P1 = 0.0;
        int param_number_P1 = -1;
        int param_number_P2 = -1;
        //bool success = convert_MC_name_to_scale_factor(mc_name, param_number, scale_factor);
        bool success = g_pg.convert_MC_name_to_param_number(mc_name, param_number_P1, param_number_P2);//, scale_factor);

        //std::cout << "param_number_P1=" << param_number_P1 << " param_number_P2=" << param_number_P2 << std::endl;


        if(success == true)
        {
            // should always be able to obtain 1 of these parameter numbers
            // but sometimes one phase is disabled hence the param number
            // is -1
            if(param_number_P1 != -1)
            {
                scale_factor_P1 = g_pg.file_params.at(param_number_P1).paramInitValue;
            }
            else
            {
                if(mc_name.find("air") == std::string::npos)
                {
                    std::cout << "WARNING: param_number_P1=" << param_number_P1 << std::endl;
                    std::cout << "mc_name=" << mc_name << std::endl;
                }
                else
                {
                    // ignore anything which contains "air"
                    // these things are only for P1
                }
                // this is typically things such as bi214_air, which is ZERO for P2
            }
            if(param_number_P2 != -1)
            {
                scale_factor_P2 = g_pg.file_params.at(param_number_P2).paramInitValue;   
            }
            else
            {
                if(mc_name.find("air") == std::string::npos)
                {
                    std::cout << "WARNING: param_number_P2=" << param_number_P2 << std::endl;
                    std::cout << "mc_name=" << mc_name << std::endl;
                }
                else
                {
                    // ignore anything which contains "air"
                    // these things are only for P1
                }
                // this is typically things such as bi214_air, which is ZERO for P2
            }
            //std::cout << "mc_name=" << mc_name << " scale_factor_P1=" << scale_factor_P1 << " scale_factor_P2=" << scale_factor_P2 << std::endl;

            // account for 208 Tl branching ratio of 36 %
            if((mc_name.find("tl208_int_rot") != std::string::npos) ||
               (mc_name.find("tl208_feShield") != std::string::npos) || // TODO: this doesn't seem to work
               (mc_name.find("tl208_pmt") != std::string::npos))
               // TODO: do not apply to tl208_air ?
            {
                if(debuglevel >= 5)
                {
                    std::cout << "mc_name=" << mc_name << " applying additional scaling factor of 0.36" << std::endl;
                }
                //std::cin.get();
                scale_factor_P1 *= 0.36;
                scale_factor_P2 *= 0.36;
                // TODO: check that this is not already applied in
                // fit_2e
                // NOTE: it isn't
            }
            

            if(debuglevel >= 3)
            {
                std::cout << "mc_name=" << mc_name << " param_number_P1=" << param_number_P1 << " param_number_P2=" << param_number_P2 << " scale_factor_P1=" << scale_factor_P1 << " scale_factor_P2=" << scale_factor_P2 << std::endl;
            }

            //std::string directory("scaled/hTotalE_/");
            std::string directory("scaled/" + theHistogram + "/");
            std::string name(theHistogram + BkgFiles[i] + "_fit_scaled");
            std::string fullname = directory + name;
            //std::string new_name(theHistogram + BkgFiles[i] + "_fit"); // TODO: probably need a different new_name for P1 and P2
            std::string new_name_P1(theHistogram + BkgFiles[i] + "_P1_fit"); // TODO: probably need a different new_name for P1 and P2
            std::string new_name_P2(theHistogram + BkgFiles[i] + "_P2_fit"); // TODO: probably need a different new_name for P1 and P2
            if(debuglevel >= 3)
            {
                std::cout << "fullname=" << fullname << std::endl;
            }

            //gDirectory->GetListOfKeys();

            //tmpHist = (TH1D*)gDirectory->Get(fullname.c_str())->Clone();
            tmpHist_P1 = (TH1D*)aFile_P1->Get(fullname.c_str())->Clone(new_name_P1.c_str());
            tmpHist_P2 = (TH1D*)aFile_P2->Get(fullname.c_str())->Clone(new_name_P2.c_str());
            // TODO: should not clone but setname() here?

            if(debuglevel >= 3)
            {
                std::cout << "check the name of histograms: (P1): " << tmpHist_P1->GetName() << " (P2): " << tmpHist_P2->GetName() << std::endl;
            }
            // TODO: suspect I need to add P1 or P2 to the histograms
            // TODO: need to be careful with some histograms because they have
            // a P1 IN, P2 IN, etc... they are different depending on P1, P2
            // as in the name is different... in fit_2e.C

            if(tmpHist_P1 != nullptr)
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
                //tmpHist = (TH1D*)gDirectory->Get(name.c_str())->Clone(hist_name.c_str());
                //tmpHist = (TH1D*)gDirectory->Get(fullname.c_str())->Clone();

                // scale by activity

                // NOTE: TODO
                // possible flaw with this method: error is no longer
                // pre-set using values from input file
                // TODO: note this in input file documentation
                // however, this may be an improvement because it
                // guarantees minuit is responsible for error estimation
                tmpHist_P1->Scale(scale_factor_P1);
                //tmpHist_P2->Scale(scale_factor); // TODO; really not at all convinced this is going to work for the histograms where P1 and P2 are seperate params in the file
                // samples are now scaled by activity
                // changed input, and pre-scaling, now need to change output

                // NOTE: Scale factor
                // samples scaled by param_init_value (activity in Bq)
                // after being read from file
                // however objects in file are scaled by
                // TotalTime / sampleNGenMC
                // also scale by 0.36 for relevant samples to account for
                // branching ratio here
                
                // in loglikelihood function the getNumberMC() functions
                // multiply the bin content by AdjustActs parameter
                // (minuit parameter)

                if((param_number_P1 == 1) || (param_number_P2 == 1))
                {
                    std::cout << "ERROR" << std::endl;
                    throw "Error";
                }

                // NOTE: do NOT apply xi reweighting here
                // this section just LOADS histograms from file and we want to LOAD
                // the default (not reweighted) nd150 spectra
                // TODO: this may no longer be true

                allMCSamples1D[channel_counter]->Add((TH1D*)tmpHist_P1);
                // TODO: does this work as expected for secular equlibrium samples?

                //std::cout << tmpHist->GetName() << std::endl;

            }
            else
            {
                std::cout << __func__ << " could not find histogram in file: " << fullname << " - disabling parameter number " << param_number_P1 << " " << param_number_P2 << std::endl;
                // cannot find histogram input data, so disable parameter
                //std::remove(enabled_params.begin(), enabled_params.end(), param_number);
                // TODO: may not re-implement this
                // if I do then I simply change the flags in the g_pg
            }

            if(tmpHist_P2 != nullptr)
            {
                tmpHist_P2->Scale(scale_factor_P2); // TODO; really not at all convinced this is going to work for the histograms where P1 and P2 are seperate params in the file

                if((param_number_P1 == 1) || (param_number_P2 == 1))
                {
                    std::cout << "ERROR" << std::endl;
                    throw "Error";
                }

                allMCSamples1D[channel_counter]->Add((TH1D*)tmpHist_P2);
            }
            else
            {
                std::cout << __func__ << " could not find histogram in file: " << fullname << " - disabling parameter number " << param_number_P1  << " " << param_number_P2 << std::endl;
                // cannot find histogram input data, so disable parameter
                //std::remove(enabled_params.begin(), enabled_params.end(), param_number);
                // TODO: may not re-implement this
                // if I do then I simply change the flags in the g_pg
            }

        }
        else
        {
            std::cerr << "success=false" << std::endl;
        }


        #if 0
        std::string search_object = MCNameToParamNameMap.at(mc_name);
        if(paramNameToNumberMap.count(search_object) > 0)
        {
            // convert from mc sample name to param number

            int param_number = paramNameToNumberMap.at(search_object);
            //std::cout << "parameber number " << param_number << " is in the paramNameToNumberMap" << std::endl;
            if(std::find(enabled_params.begin(), enabled_params.end(), param_number) != enabled_params.end())
            {
                // check if param number is enabled

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
        #endif
    }
}

// channel_counter = 0
// theChannel = "2e_"
// thePhase = "P1"
// theHistogram = "hTotalE_"
//void book1DHistograms(Int_t channel_counter, TString theChannel, TString thePhase_arg, TString theHistogram)
void book1DHistograms(Int_t channel_counter, TString theChannel, TString theHistogram)
{

    int debuglevel = 1;

    //std::cout << "booking 1D hists for " << theChannel << " " << thePhase_arg << std::endl;
    if(debuglevel >= 2)
    {
        std::cout << "book1DHistograms(" << channel_counter << ", 2e_, " << channel_histname_1D[channel_counter] << ")" << std::endl;
        std::cout << "booking 1D hists for " << theChannel << " " << "P1 and P2" << std::endl;
    }
    allMCSamples1D[channel_counter] = new TObjArray();

    //TFile *aFile = TFile::Open("/home/ebirdsall/NEMO3/Nd150_analysis/MeasureStuff/new/Macros/Nd150_" + theChannel + thePhase_arg + ".root");
    //TFile *aFile = TFile::Open("Nd150_" + theChannel + thePhase_arg + ".root");
    TString thePhase_arg_P1 = "P1";
    TString thePhase_arg_P2 = "P2";
    //TFile *aFileP1 = TFile::Open("Nd150_" + theChannel + thePhase_arg_P1 + ".root");
    //TFile *aFileP2 = TFile::Open("Nd150_" + theChannel + thePhase_arg_P2 + ".root");
    TFile *aFileP1 = new TFile("Nd150_" + theChannel + thePhase_arg_P1 + ".root");
    TFile *aFileP2 = new TFile("Nd150_" + theChannel + thePhase_arg_P2 + ".root");
    //gDirectory->cd("singleHistos");
    //gDirectory->ls();


    //TH1D *tmpHist = nullptr; //new TH1D("tmpHist_" + theChannel + thePhase_arg, "" , 1, 0, 1);

    if(debuglevel >= 2)
    {
        std::cout << "External" << std::endl;
    }
    book1DHistograms_helper(aFileP1, aFileP2,
                            channel_counter, theChannel,
                            theHistogram,
                            nExternalBkgs,
                            ExternalBkgFiles);//,
                            //tmpHist);

    // TODO: does this work as expected for secular equlibrium samples?

    if(debuglevel >= 2)
    {
        std::cout << "Internal" << std::endl;
    }
    book1DHistograms_helper(aFileP1, aFileP2,
                            channel_counter, theChannel,
                            theHistogram,
                            nInternalBkgs,
                            InternalBkgFiles);//,
                            //tmpHist);

    if(debuglevel >= 2)
    {
        std::cout << "Rn 222" << std::endl;
    }
    book1DHistograms_helper(aFileP1, aFileP2,
                            channel_counter, theChannel,
                            theHistogram,
                            nRn222Bkgs,
                            //Rn222BkgFiles);//,
                            Rn222BkgFilesNew);//,
                            //tmpHist);

    if(debuglevel >= 2)
    {
        std::cout << "Rn 220" << std::endl;
    }
    book1DHistograms_helper(aFileP1, aFileP2,
                            channel_counter, theChannel,
                            theHistogram,
                            nRn220Bkgs,
                            Rn220BkgFiles);//,
                            //tmpHist);

    if(debuglevel >= 2)
    {
        std::cout << "Neighbour" << std::endl;
    }
    book1DHistograms_helper(aFileP1, aFileP2,
                            channel_counter, theChannel,
                            theHistogram,
                            nNeighbours,
                            NeighbourFiles);//,
                            //tmpHist);

    if(debuglevel >= 2)
    {
        std::cout << "Nd150" << std::endl;
    }
    book1DHistograms_helper(aFileP1, aFileP2,
                            channel_counter, theChannel,
                            theHistogram,
                            nNd150Samples,
                            Nd150Files);//,
                            //tmpHist);

    // TODO: need P1 and P2 data here


    /*
    if(g_mode_fake_data == false)
    {
    */
        if(debuglevel >= 2)
        {
            std::cout << "Data" << std::endl;
        }
        // TODO here
        // what is name in other section of code
        //std::string name(theHistogram + "data_2e");
        //std::string directory("processeddata/hTotalE_/");
        std::string directory("processeddata/" + theHistogram + "/");
        std::string name(theHistogram + DataFile);
        //std::string fake_data_name(theHistogram + "data_2e_fake");
        std::string fullname = directory + name;
        if(debuglevel >= 5)
        {
            std::cout << "fullname=" << fullname << std::endl;
        }
        //if(gDirectory->GetListOfKeys()->Contains(fullname.c_str()))
        //TH1D *tmpHist = (TH1D*)gDirectory->Get(fullname.c_str())->Clone();
        std::string new_name_P1(theHistogram + DataFile + "_P1"); // TODO: probably need a different new_name for P1 and P2
        std::string new_name_P2(theHistogram + DataFile + "_P2"); // TODO: probably need a different new_name for P1 and P2 check works
        TH1D *tmpHist_P1 = (TH1D*)aFileP1->Get(fullname.c_str())->Clone(new_name_P1.c_str());
        TH1D *tmpHist_P2 = (TH1D*)aFileP2->Get(fullname.c_str())->Clone(new_name_P2.c_str());
                // TODO: suspect I need to add P1 or P2 to the histograms
                // TODO: need to be careful with some histograms because they have
                // a P1 IN, P2 IN, etc... they are different depending on P1, P2
                // as in the name is different... in fit_2e.C
        // TODO: should not clone but setname() here?
        if(debuglevel >= 3)
        {
            std::cout << "check the name of histograms: (P1): " << tmpHist_P1->GetName() << " (P2): " << tmpHist_P2->GetName() << std::endl;
        }
        if(tmpHist_P1 != nullptr)
        {
            //TH1D *tmpHist = nullptr;
            // 2020-04-03: removed changing of histogram name
            //std::string hist_name("data_" + theChannel + thePhase_arg);
            //std::cout << "Get() : " << name << " from file, Clone() : " << hist_name << std::endl;
            //tmpHist = (TH1D*)gDirectory->Get(name.c_str())->Clone(hist_name.c_str());
            //tmpHist = (TH1D*)gDirectory->Get(fullname.c_str())->Clone();
            allDataSamples1D->Add((TH1D*)tmpHist_P1);
        }
        else
        {
            std::cout << "gDirectory->GetListOfKeys() does not contain " << fullname << std::endl;
            std::cout << __func__ << " could not find histogram in file: " << fullname << " - disabling parameter number TODO fix this error" << std::endl;
        }

        if(tmpHist_P2 != nullptr)
        {
            allDataSamples1D->Add((TH1D*)tmpHist_P2);
        }
        else
        {
            std::cout << "gDirectory->GetListOfKeys() does not contain " << fullname << std::endl;
            std::cout << __func__ << " could not find histogram in file: " << fullname << " - disabling parameter number TODO fix this error" << std::endl;
        }
        /*
        if(gDirectory->GetListOfKeys()->Contains(theHistogram + "Data"))
        {
            std::string name(theHistogram + "Data");
            std::cout << "Get() : " << name << " from file, Clone() : " << "Data_" + theChannel + thePhase_arg << std::endl;
            tmpHist = (TH1D*)gDirectory->Get(name.c_str())->Clone("Data_" + theChannel + thePhase_arg);
            allDataSamples1D->Add((TH1D*)tmpHist);
        }
        else
        {
            std::cout << "gDirectory->GetListOfKeys() does not contain " << theHistogram + "Data" << std::endl;
        }
        */

    /*
    }
    else if(g_mode_fake_data == true)
    {
        if(debuglevel >= 2)
        {
            std::cout << "FakeData" << std::endl;
        }

        //systematic_energy_offset = 0.0;
        //rebuild_fake_data_systematics();
    }
    */



    // std::cout << tmpHist->GetName() << std::endl;
    // tmpHist->Delete();
    //aFile->Close();
    // aFile->Delete();
    aFileP1->Close();
    delete aFileP1;
    aFileP2->Close();
    delete aFileP2;
}




#endif // NEWLOGLIKFITTER_BOOK1DHISTOGRAMS_H
