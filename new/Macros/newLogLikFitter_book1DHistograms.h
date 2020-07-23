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

        // convert mc_name to scale factor
        Double_t scale_factor = 0.0;
        int param_number = -1;
        bool success = convert_MC_name_to_scale_factor(mc_name, param_number, scale_factor);

        if(success == true)
        {
            //std::string directory("scaled/hTotalE_/");
            std::string directory("scaled/" + theHistogram + "/");
            std::string name(theHistogram + BkgFiles[i] + "_fit_scaled");
            std::string fullname = directory + name;
            std::string new_name(theHistogram + BkgFiles[i] + "_fit");
            std::cout << "fullname=" << fullname << std::endl;

            //gDirectory->GetListOfKeys();

            //tmpHist = (TH1D*)gDirectory->Get(fullname.c_str())->Clone();
            tmpHist = (TH1D*)myFile->Get(fullname.c_str())->Clone(new_name.c_str());
            // TODO: should not clone but setname() here?

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
                std::cout << "could not find histogram in file: " << fullname << " - disabling parameter number " << param_number << std::endl;
                // cannot find histogram input data, so disable parameter
                std::remove(enabled_params.begin(), enabled_params.end(), param_number);
            }

        }
        else
        {
            std::cerr << "errormsg" << std::endl;
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
void book1DHistograms(Int_t channel_counter, TString theChannel, TString thePhase_arg, TString theHistogram)
{

    std::cout << "booking 1D hists for " << theChannel << " " << thePhase_arg << std::endl;
    allMCSamples1D[channel_counter] = new TObjArray();

    //TFile *aFile = TFile::Open("/home/ebirdsall/NEMO3/Nd150_analysis/MeasureStuff/new/Macros/Nd150_" + theChannel + thePhase_arg + ".root");
    TFile *aFile = TFile::Open("Nd150_" + theChannel + thePhase_arg + ".root");
    //gDirectory->cd("singleHistos");
    //gDirectory->ls();


    //TH1D *tmpHist = nullptr; //new TH1D("tmpHist_" + theChannel + thePhase_arg, "" , 1, 0, 1);
    
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
    //std::string fake_data_name(theHistogram + "data_2e_fake");
    std::string fullname = directory + name;
    std::cout << "fullname=" << fullname << std::endl;
    //if(gDirectory->GetListOfKeys()->Contains(fullname.c_str()))
    //TH1D *tmpHist = (TH1D*)gDirectory->Get(fullname.c_str())->Clone();
    TH1D *tmpHist = (TH1D*)aFile->Get(fullname.c_str())->Clone();
    if(tmpHist != nullptr)
    {
        //TH1D *tmpHist = nullptr;
        // 2020-04-03: removed changing of histogram name
        //std::string hist_name("data_" + theChannel + thePhase_arg);
        //std::cout << "Get() : " << name << " from file, Clone() : " << hist_name << std::endl;
        //tmpHist = (TH1D*)gDirectory->Get(name.c_str())->Clone(hist_name.c_str());
        //tmpHist = (TH1D*)gDirectory->Get(fullname.c_str())->Clone();
        allDataSamples1D->Add((TH1D*)tmpHist);
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
        tmpHist = (TH1D*)gDirectory->Get(name.c_str())->Clone("Data_" + theChannel + thePhase_arg);
        allDataSamples1D->Add((TH1D*)tmpHist);
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




#endif // NEWLOGLIKFITTER_BOOK1DHISTOGRAMS_H
