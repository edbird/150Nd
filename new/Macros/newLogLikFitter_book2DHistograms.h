#ifndef NEWLOGLIKFITTER_BOOK2DHISTOGRAMS_H
#define NEWLOGLIKFITTER_BOOK2DHISTOGRAMS_H


///////////////////////////////////////////////////////////////////////////////
// book2DHistograms
///////////////////////////////////////////////////////////////////////////////

void book2DHistograms_helper(
    TFile *aFile_P1,
    TFile *aFile_P2,
    Int_t channel_counter,
    TString theChannel,
    TString theHistogram,
    const int nBkgs,
    TString *BkgFiles)
    //, TH1D *tmpHist)
{
        
    TH2D *tmpHist_P1 = nullptr;
    TH2D *tmpHist_P2 = nullptr;

    for(int i = 0; i < nBkgs; i++)
    {
        // check if parameter is enabled
        // convert parameter string name to index

        // convert TString to std::string
        std::string mc_name = std::string(BkgFiles[i]);

        // convert mc_name to scale factor
        Double_t scale_factor = 0.0;
        int param_number = -1;
        bool success = g_pg.convert_MC_name_to_param_number(mc_name, param_number);


        if(success == true)
        {
            scale_factor = g_pg.file_params.at(param_number).paramInitValue;
            std::cout << "mc_name=" << mc_name << " param_number=" << param_number << " scale_factor=" << scale_factor << std::endl;

            // TODO:
            //std::string directory("scaled/hHighLowEnergy_/");
            std::string directory("scaled/" + theHistogram + "/");
            std::string name(theHistogram + BkgFiles[i] + "_fit_scaled");
            std::string fullname = directory + name;
            std::string new_name_P1(theHistogram + BkgFiles[i] + "_P1_fit");
            std::string new_name_P2(theHistogram + BkgFiles[i] + "_P2_fit");
            std::cout << "fullname=" << fullname << std::endl;

            // TODO: try catch block
            // load sample
            tmpHist_P1 = (TH2D*)aFile_P1->Get(fullname.c_str())->Clone(new_name_P1.c_str());
            tmpHist_P2 = (TH2D*)aFile_P2->Get(fullname.c_str())->Clone(new_name_P2.c_str());
            
            if(tmpHist_P1 != nullptr)
            {
                // scale by activity
                tmpHist_P1->Scale(scale_factor);

                // NOTE: do NOT apply xi reweighting here
                // this section just LOADS histograms from file and we want to LOAD
                // the default (not reweighted) nd150 spectra

                if(param_number == 1)
                {
                    std::cout << "ERROR" << std::endl;
                    throw std::runtime_error("param_number=1");
                }

                allMCSamples2D[channel_counter]->Add(tmpHist_P1);
            }
            else
            {
                std::cout << __func__ << " could not find histogram in file: " << fullname << " - disabling parameter number " << param_number << std::endl;
                // cannot find histogram input data, so disable parameter
                //std::remove(enabled_params.begin(), enabled_params.end(), param_number); // TODO
            }

            if(tmpHist_P2 != nullptr)
            {
                // scale by activity
                tmpHist_P2->Scale(scale_factor);

                // NOTE: do NOT apply xi reweighting here
                // this section just LOADS histograms from file and we want to LOAD
                // the default (not reweighted) nd150 spectra

                if(param_number == 1)
                {
                    std::cout << "ERROR" << std::endl;
                    throw std::runtime_error("param_number=1");
                }

                allMCSamples2D[channel_counter]->Add(tmpHist_P2);
            }
            else
            {
                std::cout << __func__ << " could not find histogram in file: " << fullname << " - disabling parameter number " << param_number << std::endl;
                // cannot find histogram input data, so disable parameter
                //std::remove(enabled_params.begin(), enabled_params.end(), param_number); // TODO
            }
        }
        else
        {
            std::cerr << "success=false" << std::endl;
        }
    }













        #if 0


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
                tmpHist = (TH2D*)myFile->Get(fullname.c_str())->Clone(new_name.c_str());

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
        #endif
        #if 0
        else
        {
            std::cout << "!!!!! ERROR: search_object=" << search_object << " not found in paramNameToNumberMap" << std::endl;
            std::cout << "mc_name=" << mc_name << std::endl;

            print_map(paramNameToNumberMap, "paramNameToNumberMap");
            print_map(MCNameToParamNameMap, "MCNameToParamNameMap");
        }
        #endif
    //}
}



// channel_counter = 0
// theChannel = "2e_"
// thePhase = "P1"
// theHistogram = "hHighLowEnergy_"
void book2DHistograms(Int_t channel_counter, TString theChannel, TString theHistogram) {

    std::cout << "booking 2D hists for " << theChannel << " " << "P1 and P2" << std::endl;
    allMCSamples2D[channel_counter] = new TObjArray();

    TString thePhase_arg_P1 = "P1";
    TString thePhase_arg_P2 = "P2";
    TFile *aFileP1 = TFile::Open("Nd150_" + theChannel + thePhase_arg_P1 + ".root");
    TFile *aFileP2 = TFile::Open("Nd150_" + theChannel + thePhase_arg_P2 + ".root");
    
    std::cout << "External" << std::endl;
    book2DHistograms_helper(aFileP1, aFileP2,
                            channel_counter, theChannel,
                            theHistogram,
                            nExternalBkgs,
                            ExternalBkgFiles);

    std::cout << "Internal" << std::endl;
    book2DHistograms_helper(aFileP1, aFileP2,
                            channel_counter, theChannel,
                            theHistogram,
                            nInternalBkgs,
                            InternalBkgFiles);

    std::cout << "Rn 222" << std::endl;
    book2DHistograms_helper(aFileP1, aFileP2,
                            channel_counter, theChannel,
                            theHistogram,
                            nRn222Bkgs,
                            //Rn222BkgFiles);//,
                            Rn222BkgFilesNew);

    std::cout << "Rn 220" << std::endl;
    book2DHistograms_helper(aFileP1, aFileP2,
                            channel_counter, theChannel,
                            theHistogram,
                            nRn220Bkgs,
                            Rn220BkgFiles);

    std::cout << "Neighbour" << std::endl;
    book2DHistograms_helper(aFileP1, aFileP2,
                            channel_counter, theChannel,
                            theHistogram,
                            nNeighbours,
                            NeighbourFiles);

    std::cout << "Nd150" << std::endl;
    book2DHistograms_helper(aFileP1, aFileP2,
                            channel_counter, theChannel,
                            theHistogram,
                            nNd150Samples,
                            Nd150Files);

    std::cout << "Data" << std::endl;
    // TODO here
    // what is name in other section of code
    //std::string name(theHistogram + "data_2e");
    std::string directory("processeddata/" + theHistogram + "/");
    std::string name(theHistogram + DataFile);
    //std::string fake_data_name(theHistogram + "data_2e_fake");
    std::string fullname = directory + name;
    std::cout << "fullname=" << fullname << std::endl;

    // load histogram from file
    // TODO: try catch block
    std::string new_name_P1(theHistogram + DataFile + "_P1"); // TODO: probably need a different new_name for P1 and P2
    std::string new_name_P2(theHistogram + DataFile + "_P2"); // TODO: probably need a different new_name for P1 and P2 check works
    TH1D *tmpHist_P1 = (TH1D*)aFileP1->Get(fullname.c_str())->Clone(new_name_P1.c_str());
    TH1D *tmpHist_P2 = (TH1D*)aFileP2->Get(fullname.c_str())->Clone(new_name_P2.c_str());
    //TH1D *tmpHist = (TH1D*)aFile->Get(fullname.c_str())->Clone();
    if(tmpHist_P1 != nullptr)
    {
        allDataSamples2D->Add((TH2D*)tmpHist_P1);
    }
    else
    {
        std::cout << "gDirectory->GetListOfKeys() does not contain " << fullname << std::endl;
    }

    if(tmpHist_P2 != nullptr)
    {
        allDataSamples2D->Add((TH2D*)tmpHist_P2);
    }
    else
    {
        std::cout << __func__ << " ERROR: could not find " << fullname << " in file" << std::endl;
    }

    // aFile->Close();
    // aFile->Delete();
}


#endif // NEWLOGLIKFITTER_BOOK2DHISTOGRAMS_H
