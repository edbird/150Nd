#ifndef NEWLOGLIKFITTER_BOOK2DHISTOGRAMS_H
#define NEWLOGLIKFITTER_BOOK2DHISTOGRAMS_H


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
                    get_paramInitValueError(thePhase, param_number, param_init_value, param_init_error);
                    /*
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
                    */
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



#endif // NEWLOGLIKFITTER_BOOK2DHISTOGRAMS_H
