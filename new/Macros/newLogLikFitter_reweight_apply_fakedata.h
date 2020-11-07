#ifndef NEWLOGLIKFITTER_REWEIGHT_APPLY_FAKEDATA_H
#define NEWLOGLIKFITTER_REWEIGHT_APPLY_FAKEDATA_H




void reweight_apply_fakedata(
    // P1
    TH1D *&hTotalE_output_P1,
    TH1D *&hSingleEnergy_output_P1,
    TH1D *&hHighEnergy_output_P1,
    TH1D *&hLowEnergy_output_P1,
    TH1D *&hEnergySum_output_P1,
    TH1D *&hEnergyDiff_output_P1,
    TH2D *&hHighLowEnergy_output_P1,
    // P2
    TH1D *&hTotalE_output_P2,
    TH1D *&hSingleEnergy_output_P2,
    TH1D *&hHighEnergy_output_P2,
    TH1D *&hLowEnergy_output_P2,
    TH1D *&hEnergySum_output_P2,
    TH1D *&hEnergyDiff_output_P2,
    TH2D *&hHighLowEnergy_output_P2,
    // additional parameters
    const Double_t xi_31,
    const Double_t xi_31_baseline,
    const TH2D* const h_nEqNull,
    const TH2D* const h_nEqTwo,
    const Double_t psiN0,
    const Double_t psiN2,
    const Double_t bb_Q)
{

    ///////////////////////////////////////////////////////////////////////////
    // init map
    // convert mc name to path
    ///////////////////////////////////////////////////////////////////////////

    std::map<TString, TString> MCNameToPathMap;
    for(Int_t index = 0; index < nRn222Bkgs; ++ index)
    {
        TString mcname = Rn222BkgFilesNew[index];
        TString mcpath = "externals/";
        MCNameToPathMap.insert(std::make_pair(mcname, mcpath));
    }
    for(Int_t index = 0; index < nRn220Bkgs; ++ index)
    {
        TString mcname = Rn220BkgFiles[index];
        TString mcpath = "externals/";
        MCNameToPathMap.insert(std::make_pair(mcname, mcpath));
    }
    for(Int_t index = 0; index < nExternalBkgs; ++ index)
    {
        TString mcname = ExternalBkgFiles[index];
        TString mcpath = "externals/";
        MCNameToPathMap.insert(std::make_pair(mcname, mcpath));
    }
    for(Int_t index = 0; index < nInternalBkgs; ++ index)
    {
        TString mcname = InternalBkgFiles[index];
        TString mcpath = "internals/";
        MCNameToPathMap.insert(std::make_pair(mcname, mcpath));
    }
    for(Int_t index = 0; index < nNeighbours; ++ index)
    {
        TString mcname = NeighbourFiles[index];
        TString mcpath = "neighbours/";
        MCNameToPathMap.insert(std::make_pair(mcname, mcpath));
    }
    for(Int_t index = 0; index < nNd150Samples; ++ index)
    {
        TString mcname = Nd150Files[index];
        TString mcpath = "nd150/";
        MCNameToPathMap.insert(std::make_pair(mcname, mcpath));
    }


    ///////////////////////////////////////////////////////////////////////////
    // start
    ///////////////////////////////////////////////////////////////////////////

    int debuglevel = 1;

    if(debuglevel >= 3)
    {
        std::cout << __func__ << std::endl;
    }

    TString tinput_filename = "Nd150_2eNg_output_postprocessed_small.root";


    ///////////////////////////////////////////////////////////////////////////
    // delete "total" fake data histograms
    ///////////////////////////////////////////////////////////////////////////

    const TString name_append = "";
    const TString sampleName = "fakedata";

    if(debuglevel >= 6)
    {
        std::cout << "calling deleter_helper" << std::endl;
    }

    // P1
    deleter_helper("hTotalE_" + sampleName + name_append + "_P1");
    deleter_helper("hSingleEnergy_" + sampleName + name_append + "_P1");
    deleter_helper("hHighEnergy_" + sampleName + name_append + "_P1");
    deleter_helper("hLowEnergy_" + sampleName + name_append + "_P1");
    deleter_helper("hEnergySum_" + sampleName + name_append + "_P1");
    deleter_helper("hEnergyDiff_" + sampleName + name_append + "_P1");
    deleter_helper("hHighLowEnergy_" + sampleName + name_append + "_P1");
    // P2
    deleter_helper("hTotalE_" + sampleName + name_append + "_P2");
    deleter_helper("hSingleEnergy_" + sampleName + name_append + "_P2");
    deleter_helper("hHighEnergy_" + sampleName + name_append + "_P2");
    deleter_helper("hLowEnergy_" + sampleName + name_append + "_P2");
    deleter_helper("hEnergySum_" + sampleName + name_append + "_P2");
    deleter_helper("hEnergyDiff_" + sampleName + name_append + "_P2");
    deleter_helper("hHighLowEnergy_" + sampleName + name_append + "_P2");

    if(debuglevel >= 6)
    {
        std::cout << "done calling deleter_helper" << std::endl;
    }



    if(debuglevel >= 6)
    {
        std::cout << "allocating new histograms" << std::endl;
    }

    ///////////////////////////////////////////////////////////////////
    // phase 1 histograms (output)
    ///////////////////////////////////////////////////////////////////

    hTotalE_output_P1          = new TH1D("hTotalE_" + sampleName + name_append + "_P1",
                                TString("Phase ") + "P1" + " " + sampleName + name_append + " total energy; Total Energy #SigmaE_{e} (MeV)",
                                50, 0.0, 5.0);

    hSingleEnergy_output_P1    = new TH1D("hSingleEnergy_" + sampleName + name_append + "_P1",
                                TString("Phase ") + "P1" + " " + sampleName  + name_append + " Single Energy",
                                50, 0.0, 5.0);
    
    //std::cout << "High new" << std::endl;
    hHighEnergy_output_P1     = new TH1D("hHighEnergy_" + sampleName + name_append + "_P1",
                                TString("Phase ") + "P1" + " " + sampleName + name_append + " High Energy; Energy (MeV)",
                                50, 0.0, 5.0);

    //std::cout << "Low new" << std::endl;
    hLowEnergy_output_P1     = new TH1D("hLowEnergy_" + sampleName + name_append + "_P1",
                                TString("Phase ") + "P1" + " " + sampleName + name_append + " Low Energy",
                                50, 0.0, 5.0);

    //std::cout << "Sum new" << std::endl;
    hEnergySum_output_P1     = new TH1D("hEnergySum_" + sampleName + name_append + "_P1",
                                TString("Phase ") + "P1" + " " + sampleName + name_append + " Low Energy",
                                50, 0.0, 5.0);

    //std::cout << "Diff new" << std::endl;
    hEnergyDiff_output_P1     = new TH1D("hEnergyDiff_" + sampleName + name_append + "_P1",
                                TString("Phase ") + "P1" + " " + sampleName + name_append + " Low Energy",
                                50, 0.0, 5.0);

    //std::cout << "HighLow new" << std::endl;
    hHighLowEnergy_output_P1     = new TH2D("hHighLowEnergy_" + sampleName + name_append + "_P1",
                                TString("Phase ") + "P1" + " " + sampleName + name_append + ";Low Energy Electron Energy (MeV);High Energy Electron Energy (MeV)",
                                50, 0.0, 5.0, 50, 0.0, 5.0);


    ///////////////////////////////////////////////////////////////////
    // phase 2 histograms (output)
    ///////////////////////////////////////////////////////////////////

    hTotalE_output_P2           = new TH1D("hTotalE_" + sampleName + name_append + "_P2",
                                //"Phase " + Phase + " " + sampleName + name_append + " total energy; #SigmaE_{e} (MeV)",
                                TString("Phase ") + "P2" + " " + sampleName + name_append + " total energy; Total Energy #SigmaE_{e} (MeV)",
                                50, 0.0, 5.0);

    //std::cout << "Single new" << std::endl;
    hSingleEnergy_output_P2    = new TH1D("hSingleEnergy_" + sampleName + name_append + "_P2",
                                TString("Phase ") + "P2" + " " + sampleName  + name_append + " Single Energy",
                                50, 0.0, 5.0);
    
    //std::cout << "High new" << std::endl;
    hHighEnergy_output_P2     = new TH1D("hHighEnergy_" + sampleName + name_append + "_P2",
                                TString("Phase ") + "P2" + " " + sampleName + name_append + " High Energy; Energy (MeV)",
                                50, 0.0, 5.0);

    //std::cout << "Low new" << std::endl;
    hLowEnergy_output_P2     = new TH1D("hLowEnergy_" + sampleName + name_append + "_P2",
                                TString("Phase ") + "P2" + " " + sampleName + name_append + " Low Energy",
                                50, 0.0, 5.0);

    //std::cout << "Sum new" << std::endl;
    hEnergySum_output_P2     = new TH1D("hEnergySum_" + sampleName + name_append + "_P2",
                                TString("Phase ") + "P2" + " " + sampleName + name_append + " Low Energy",
                                50, 0.0, 5.0);

    //std::cout << "Diff new" << std::endl;
    hEnergyDiff_output_P2     = new TH1D("hEnergyDiff_" + sampleName + name_append + "_P2",
                                TString("Phase ") + "P2" + " " + sampleName + name_append + " Low Energy",
                                50, 0.0, 5.0);

    //std::cout << "HighLow new" << std::endl;
    hHighLowEnergy_output_P2     = new TH2D("hHighLowEnergy_" + sampleName + name_append + "_P2",
                                TString("Phase ") + "P2" + " " + sampleName + name_append + ";Low Energy Electron Energy (MeV);High Energy Electron Energy (MeV)",
                                50, 0.0, 5.0, 50, 0.0, 5.0);


    if(debuglevel >= 6)
    {
        std::cout << "done allocating new histograms" << std::endl;
    }



    hTotalE_output_P1->Sumw2();
    hSingleEnergy_output_P1->Sumw2();
    hHighEnergy_output_P1->Sumw2();
    hLowEnergy_output_P1->Sumw2();
    hHighLowEnergy_output_P1->Sumw2();
    hEnergySum_output_P1->Sumw2();
    hEnergyDiff_output_P1->Sumw2();


    hTotalE_output_P2->Sumw2();
    hSingleEnergy_output_P2->Sumw2();
    hHighEnergy_output_P2->Sumw2();
    hLowEnergy_output_P2->Sumw2();
    hHighLowEnergy_output_P2->Sumw2();
    hEnergySum_output_P2->Sumw2();
    hEnergyDiff_output_P2->Sumw2();

    
    ///////////////////////////////////////////////////////////////////////////
    // code to load data from file
    ///////////////////////////////////////////////////////////////////////////

    // loop over all the parameters
    std::map<int, file_parameter>::iterator it{g_pg.file_params.begin()};
    for(; it != g_pg.file_params.end(); ++ it)
    {

        
        int paramNumber = it->first;
        bool paramEnabled = it->second.paramEnabled;
        bool paramEnabledP1 = it->second.paramEnabledP1;
        bool paramEnabledP2 = it->second.paramEnabledP2;
        double paramInitValue = it->second.paramInitValue;

        if(debuglevel >= 4)
        {
            std::cout << "paramNumber=" << paramNumber << std::endl;
        }
        if(debuglevel >= 5)
        {
            std::cout << "paramEnabled=" << paramEnabled << std::endl;
            std::cout << "paramEnabledP1=" << paramEnabledP1 << std::endl;
            std::cout << "paramEnabledP2=" << paramEnabledP2 << std::endl;
        }

        bool ok = false;
        if(paramEnabled == true)
        {
            if(gEnablePhase1 == true)
            {
                if(paramEnabledP1 == true)
                {
                    if(debuglevel >= 9)
                    {
                        std::cout << "set: ok=true" << std::endl;
                    }
                    ok = true;
                }
            }

            if(gEnablePhase2 == true)
            {
                if(paramEnabledP2 == true)
                {
                    if(debuglevel >= 9)
                    {
                        std::cout << "set: ok=true" << std::endl;
                    }
                    ok = true;
                }
            }
        }
        else
        {
            ok = false;
            if(debuglevel >= 8)
            {
                std::cout << "ok=false, continue" << std::endl;
            }
            continue;
        }
        if(ok == false)
        {
            if(debuglevel > 0)
            {
                std::cout << __func__ << " ok == false" << std::endl;
                std::cin.get();
            }
            if(debuglevel >= 8)
            {
                std::cout << "ok=false, continue" << std::endl;
            }
            continue;
        }

        if(debuglevel >= 7)
        {
            std::cout << "ok=true" << std::endl;
        }

        // loop over all MC names
        std::vector<std::string>::iterator mc_name_it{it->second.MCNameList.begin()};
        for(; mc_name_it != it->second.MCNameList.end(); ++ mc_name_it)
        {
            const TString sampleName = *mc_name_it;
            //const TString name_append = ""; // TODO: this has to be something for the samples where we split the data by IN/OUT/L0/etc ?
            // unless I split in fit_2e code
            // 2020-08-30: done

            // NOTE: set name_append to "_fakedata_tmp" for memory management
            // purposes
            const TString name_append = "_fakedata_tmp";

            if(debuglevel >= 4)
            {
                std::cout << "sampleName=" << sampleName << std::endl;
            }
            if(debuglevel >= 5)
            {
                std::cout << "name_append=" << name_append << std::endl;
            }

            // load from files applying systematics

            // normal procedure is to fill histograms with weight variable
            // weight is usually 1, however it is not 1 for some MC samples
            // create histograms for hTotalE, hSingleEnergy, etc...
            // then scale these histograms
            // and ->Add() to total to create fake data

            // P1
            if(paramEnabledP1)
            {
                if(debuglevel >= 8)
                {
                    std::cout << "P1 deleter helper" << std::endl;
                }
                deleter_helper("hTotalE_" + sampleName + name_append + "_P1");
                deleter_helper("hSingleEnergy_" + sampleName + name_append + "_P1");
                deleter_helper("hHighEnergy_" + sampleName + name_append + "_P1");
                deleter_helper("hLowEnergy_" + sampleName + name_append + "_P1");
                deleter_helper("hEnergySum_" + sampleName + name_append + "_P1");
                deleter_helper("hEnergyDiff_" + sampleName + name_append + "_P1");
                deleter_helper("hHighLowEnergy_" + sampleName + name_append + "_P1");
            }
            // P2
            if(paramEnabledP2)
            {
                if(debuglevel >= 8)
                {
                    std::cout << "P2 deleter helper" << std::endl;
                }
                deleter_helper("hTotalE_" + sampleName + name_append + "_P2");
                deleter_helper("hSingleEnergy_" + sampleName + name_append + "_P2");
                deleter_helper("hHighEnergy_" + sampleName + name_append + "_P2");
                deleter_helper("hLowEnergy_" + sampleName + name_append + "_P2");
                deleter_helper("hEnergySum_" + sampleName + name_append + "_P2");
                deleter_helper("hEnergyDiff_" + sampleName + name_append + "_P2");
                deleter_helper("hHighLowEnergy_" + sampleName + name_append + "_P2");
            }

            
            TH1D* hTotalE_P1_tmp = nullptr;
            TH1D* hSingleEnergy_P1_tmp = nullptr;
            TH1D* hHighEnergy_P1_tmp = nullptr;
            TH1D* hLowEnergy_P1_tmp = nullptr;
            TH2D* hHighLowEnergy_P1_tmp = nullptr;
            TH1D* hEnergySum_P1_tmp = nullptr;
            TH1D* hEnergyDiff_P1_tmp = nullptr;

            TH1D* hTotalE_P2_tmp = nullptr;
            TH1D* hSingleEnergy_P2_tmp = nullptr;
            TH1D* hHighEnergy_P2_tmp = nullptr;
            TH1D* hLowEnergy_P2_tmp = nullptr;
            TH2D* hHighLowEnergy_P2_tmp = nullptr;
            TH1D* hEnergySum_P2_tmp = nullptr;
            TH1D* hEnergyDiff_P2_tmp = nullptr;
            
            ///////////////////////////////////////////////////////////////////
            // phase 1 histograms (tmp)
            ///////////////////////////////////////////////////////////////////

            if(paramEnabledP1)
            {
                if(debuglevel >= 8)
                {
                    std::cout << "P1 alloc" << std::endl;
                }
                hTotalE_P1_tmp          = new TH1D("hTotalE_" + sampleName + name_append + "_P1",
                                            TString("Phase ") + "P1" + " " + sampleName + name_append + " total energy; Total Energy #SigmaE_{e} (MeV)",
                                            50, 0.0, 5.0);

                hSingleEnergy_P1_tmp    = new TH1D("hSingleEnergy_" + sampleName + name_append + "_P1",
                                            TString("Phase ") + "P1" + " " + sampleName  + name_append + " Single Energy",
                                            50, 0.0, 5.0);
                
                //std::cout << "High new" << std::endl;
                hHighEnergy_P1_tmp     = new TH1D("hHighEnergy_" + sampleName + name_append + "_P1",
                                            TString("Phase ") + "P1" + " " + sampleName + name_append + " High Energy; Energy (MeV)",
                                            50, 0.0, 5.0);

                //std::cout << "Low new" << std::endl;
                hLowEnergy_P1_tmp     = new TH1D("hLowEnergy_" + sampleName + name_append + "_P1",
                                            TString("Phase ") + "P1" + " " + sampleName + name_append + " Low Energy",
                                            50, 0.0, 5.0);

                //std::cout << "Sum new" << std::endl;
                hEnergySum_P1_tmp     = new TH1D("hEnergySum_" + sampleName + name_append + "_P1",
                                            TString("Phase ") + "P1" + " " + sampleName + name_append + " Low Energy",
                                            50, 0.0, 5.0);

                //std::cout << "Diff new" << std::endl;
                hEnergyDiff_P1_tmp     = new TH1D("hEnergyDiff_" + sampleName + name_append + "_P1",
                                            TString("Phase ") + "P1" + " " + sampleName + name_append + " Low Energy",
                                            50, 0.0, 5.0);

                //std::cout << "HighLow new" << std::endl;
                hHighLowEnergy_P1_tmp     = new TH2D("hHighLowEnergy_" + sampleName + name_append + "_P1",
                                            TString("Phase ") + "P1" + " " + sampleName + name_append + ";Low Energy Electron Energy (MeV);High Energy Electron Energy (MeV)",
                                            50, 0.0, 5.0, 50, 0.0, 5.0);
            }


            ///////////////////////////////////////////////////////////////////
            // phase 2 histograms (tmp)
            ///////////////////////////////////////////////////////////////////

            if(paramEnabledP2)
            {
                if(debuglevel >= 8)
                {
                    std::cout << "P2 alloc" << std::endl;
                }
                hTotalE_P2_tmp = new TH1D("hTotalE_" + sampleName + name_append + "_P2",
                                   //"Phase " + Phase + " " + sampleName + name_append + " total energy; #SigmaE_{e} (MeV)",
                                   TString("Phase ") + "P2" + " " + sampleName + name_append + " total energy; Total Energy #SigmaE_{e} (MeV)",
                                   50, 0.0, 5.0);
                                   // TODO: changed from 4

                //std::cout << "Single new" << std::endl;
                hSingleEnergy_P2_tmp    = new TH1D("hSingleEnergy_" + sampleName + name_append + "_P2",
                                            TString("Phase ") + "P2" + " " + sampleName  + name_append + " Single Energy",
                                            50, 0.0, 5.0);
                
                //std::cout << "High new" << std::endl;
                hHighEnergy_P2_tmp     = new TH1D("hHighEnergy_" + sampleName + name_append + "_P2",
                                            TString("Phase ") + "P2" + " " + sampleName + name_append + " High Energy; Energy (MeV)",
                                            50, 0.0, 5.0);

                //std::cout << "Low new" << std::endl;
                hLowEnergy_P2_tmp     = new TH1D("hLowEnergy_" + sampleName + name_append + "_P2",
                                            TString("Phase ") + "P2" + " " + sampleName + name_append + " Low Energy",
                                            50, 0.0, 5.0);

                //std::cout << "Sum new" << std::endl;
                hEnergySum_P2_tmp     = new TH1D("hEnergySum_" + sampleName + name_append + "_P2",
                                            TString("Phase ") + "P2" + " " + sampleName + name_append + " Low Energy",
                                            50, 0.0, 5.0);

                //std::cout << "Diff new" << std::endl;
                hEnergyDiff_P2_tmp     = new TH1D("hEnergyDiff_" + sampleName + name_append + "_P2",
                                            TString("Phase ") + "P2" + " " + sampleName + name_append + " Low Energy",
                                            50, 0.0, 5.0);

                //std::cout << "HighLow new" << std::endl;
                hHighLowEnergy_P2_tmp     = new TH2D("hHighLowEnergy_" + sampleName + name_append + "_P2",
                                            TString("Phase ") + "P2" + " " + sampleName + name_append + ";Low Energy Electron Energy (MeV);High Energy Electron Energy (MeV)",
                                            50, 0.0, 5.0, 50, 0.0, 5.0);
            }

            
            if(paramEnabledP1)
            {
                hTotalE_P1_tmp->Sumw2();
                hSingleEnergy_P1_tmp->Sumw2();
                hHighEnergy_P1_tmp->Sumw2();
                hLowEnergy_P1_tmp->Sumw2();
                hHighLowEnergy_P1_tmp->Sumw2();
                hEnergySum_P1_tmp->Sumw2();
                hEnergyDiff_P1_tmp->Sumw2();
            }


            if(paramEnabledP2)
            {
                hTotalE_P2_tmp->Sumw2();
                hSingleEnergy_P2_tmp->Sumw2();
                hHighEnergy_P2_tmp->Sumw2();
                hLowEnergy_P2_tmp->Sumw2();
                hHighLowEnergy_P2_tmp->Sumw2();
                hEnergySum_P2_tmp->Sumw2();
                hEnergyDiff_P2_tmp->Sumw2();
            }


            ///////////////////////////////////////////////////////////////////
            // load data from file
            ///////////////////////////////////////////////////////////////////

            TString filePath = "/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
            TString thePath = MCNameToPathMap.at(sampleName);

            TFile *inputFile = new TFile(filePath + thePath + sampleName + "/Nd150_2eNg_output_postprocessed_small.root");
            if(!inputFile->IsOpen())
            {
                std::cout << "failed to open" << std::endl;
            }
            TTree *inputTree = (TTree*)inputFile->Get("Nd150_2eNg/Nd150_2eNg");
            if(debuglevel >= 8)
            {
                std::cout << "read tree, " << filePath + thePath + sampleName + "/Nd150_2eNg_output_postprocessed_small.root" << std::endl;
            }

            int run;
            double eventTime;
            double radonWeight;
            double bi210Weight;
            double electronEnergy[2];
            double trueElectronEnergy[2];

            inputTree->SetBranchAddress("Run", &run);
            inputTree->SetBranchAddress("eventTime", &eventTime);
            inputTree->SetBranchAddress("radonWeight", &radonWeight);
            inputTree->SetBranchAddress("bi210Weight", &bi210Weight);
            inputTree->SetBranchAddress("electronEnergy", electronEnergy);
            inputTree->SetBranchAddress("trueElectronEnergy", trueElectronEnergy);

            double el_energy_0;
            double el_energy_1;

            double weight_total_150Nd = 0.0;
            double event_pass_count_150Nd = 0.0;

            Long_t events = (Long_t)inputTree->GetEntries();
            for(Long_t event_i = 0; event_i < events; ++ event_i)
            {
            
                inputTree->GetEvent(event_i);

                
                ///////////////////////////////////////////////////////////////
                // see file: fit_2e.C for documentation
                ///////////////////////////////////////////////////////////////
                double weight = 1.0;
                if(sampleName.CompareTo("bi214_swire") == 0)
                {
                    weight = radonWeight;
                }
                else if(sampleName.CompareTo("pb214_swire") == 0)
                {
                    weight = radonWeight;
                }
                else if(sampleName.Contains("bi210"))
                {
                    // halflife
                    // https://periodictable.com/Isotopes/083.210/index.p.full.dm.html
                    const double T21 = 86400. * 5.011574074074;
                    // lambda
                    const double lambda = std::log(2.0) / T21;
                    double arg = -lambda * eventTime;
                    // I assume event time is in SI units
                    weight = bi210Weight * std::exp(arg);
                }
                else if(sampleName.Contains("co60"))
                {
                    // halflife
                    // https://periodictable.com/Isotopes/027.60/index.p.html
                    const double T21 = 31557600. * 5.274923896499;
                    // lambda
                    const double lambda = std::log(2.0) / T21;
                    double arg = -lambda * eventTime;
                    // I assume event time is in SI units
                    weight = std::exp(arg);
                }
                else if(sampleName.Contains("eu152"))
                {
                    // halflife
                    // https://periodictable.com/Isotopes/063.152/index.p.html
                    const double T21 = 31557600. * 13.537;
                    // lambda
                    const double lambda = std::log(2.0) / T21;
                    double arg = -lambda * eventTime;
                    // I assume event time is in SI units
                    weight = std::exp(arg);
                }
                else if(sampleName.Contains("eu154"))
                {
                    // halflife
                    // https://periodictable.com/Isotopes/063.154/index.html
                    const double T21 = 31557600. * 8.593;
                    // lambda
                    const double lambda = std::log(2.0) / T21;
                    double arg = -lambda * eventTime;
                    // I assume event time is in SI units
                    weight = std::exp(arg);
                }



                int highE_index = -1;
                int lowE_index = -1;
                if(electronEnergy[0] > electronEnergy[1])
                {
                    highE_index = 0;
                    lowE_index = 1;
                }
                else
                {
                    highE_index = 1;
                    lowE_index = 0;
                }



                if(sampleName.Contains("nd150"))
                {
                    const Double_t trueT0{trueElectronEnergy[highE_index]};
                    const Double_t trueT1{trueElectronEnergy[lowE_index]};

                    Double_t _weight{ReWeight3(trueT0, trueT1, xi_31_baseline, xi_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true")}; // TODO remove true?

                    /*
                    if(_weight != 1.0)
                    {
                        std::cout << "_weight=" << _weight << std::endl;
                    }
                    */

                    weight *= _weight;

                    //weight_total_150Nd += _weight;
                }


                // set phase
                // 1869 - 3395 is P1
                // 3396 - 9186 is P2
                bool cut = true;
                Int_t thePhase = -1;
                if((1869 <= run) && (run <= 3395))
                {
                    cut = false;
                    thePhase = 0;
                }
                else if((3396 <= run) && (run <= 9186))
                {
                    cut = false;
                    thePhase = 1;
                }
                if(cut == true) continue;


                ///////////////////////////////////////////////////////////////
                // apply systematics
                ///////////////////////////////////////////////////////////////

                Double_t el_energy_0{electronEnergy[highE_index]};
                Double_t el_energy_1{electronEnergy[lowE_index]};

                // energy scale factor (multiply)
                el_energy_0 = el_energy_0 * (1.0 + gSystematics.systematic_energy_scale);
                el_energy_1 = el_energy_1 * (1.0 + gSystematics.systematic_energy_scale);

                //std::cout << "gSystematics.systematic_energy_offset=" << gSystematics.systematic_energy_offset << std::endl;
                //std::cout << "el_energy_0=" << el_energy_0 << std::endl;
                //std::cout << "el_energy_1=" << el_energy_1 << std::endl;

                // linear energy offset systematic
                el_energy_0 = el_energy_0 + gSystematics.systematic_energy_offset;
                el_energy_1 = el_energy_1 + gSystematics.systematic_energy_offset;

                // linear energy offsetsmall systematic
                el_energy_0 = el_energy_0 + gSystematics.systematic_energy_offsetsmall;
                el_energy_1 = el_energy_1 + gSystematics.systematic_energy_offsetsmall;

                //std::cout << "-> el_energy_0=" << el_energy_0 << std::endl;
                //std::cout << "-> el_energy_1=" << el_energy_1 << std::endl;

                // efficiency and enrichment
                // if sample is nd150, then apply both
                if(sampleName.Contains("nd150"))
                {
                //std::cout << "sampleName=" << sampleName << std::endl;
                //std::cin.get();
                //if(gSystematics.systematic_enrichment != 0.0)
                //{
                //std::cout << gSystematics.systematic_enrichment << std::endl;
                //std::cout << "weight=" << weight << std::endl;
                //}
                    weight = weight *
                        (1.0 + gSystematics.systematic_efficiency
                             + gSystematics.systematic_enrichment);
                             
                //if(gSystematics.systematic_enrichment != 0.0)
                //{
                //std::cout << "weight=" << weight << std::endl;
                //std::cin.get();
                //}
                }
                // if some other MC sample, then only apply efficiency effect
                else
                {
                    weight = weight *
                        (1.0 + gSystematics.systematic_efficiency);
                }


                if(gSystematics.systematic_foil_thickness_virtual != 0.0)
                {
                    if(gSystematics.systematic_foil_thickness_virtual > 0.0)
                    {
                        Int_t bin_x{gSystematics.h_systematic_foil_thickness_virtual_h->GetXaxis()->FindBin(el_energy_0)};
                        Int_t bin_y{gSystematics.h_systematic_foil_thickness_virtual_h->GetYaxis()->FindBin(el_energy_1)};
                        Double_t content{gSystematics.h_systematic_foil_thickness_virtual_h->GetBinContent(bin_x, bin_y)};
                        weight = weight * (gSystematics.systematic_foil_thickness_virtual * content);
                    }
                    else if(gSystematics.systematic_foil_thickness_virtual < 0.0)
                    {
                        Int_t bin_x{gSystematics.h_systematic_foil_thickness_virtual_l->GetXaxis()->FindBin(el_energy_0)};
                        Int_t bin_y{gSystematics.h_systematic_foil_thickness_virtual_l->GetYaxis()->FindBin(el_energy_1)};
                        Double_t content{gSystematics.h_systematic_foil_thickness_virtual_l->GetBinContent(bin_x, bin_y)};
                        weight = weight * (std::abs(gSystematics.systematic_foil_thickness_virtual) * content);
                    }
                }

                if(gSystematics.systematic_dEdX_virtual != 0.0)
                {
                    if(gSystematics.systematic_dEdX_virtual > 0.0)
                    {
                        Int_t bin_x{gSystematics.h_systematic_dEdX_virtual_h->GetXaxis()->FindBin(el_energy_0)};
                        Int_t bin_y{gSystematics.h_systematic_dEdX_virtual_h->GetYaxis()->FindBin(el_energy_1)};
                        Double_t content{gSystematics.h_systematic_dEdX_virtual_h->GetBinContent(bin_x, bin_y)};
                        weight = weight * (gSystematics.systematic_dEdX_virtual * content);
                    }
                    else if(gSystematics.systematic_dEdX_virtual < 0.0)
                    {
                        Int_t bin_x{gSystematics.h_systematic_dEdX_virtual_l->GetXaxis()->FindBin(el_energy_0)};
                        Int_t bin_y{gSystematics.h_systematic_dEdX_virtual_l->GetYaxis()->FindBin(el_energy_1)};
                        Double_t content{gSystematics.h_systematic_dEdX_virtual_l->GetBinContent(bin_x, bin_y)};
                        weight = weight * (std::abs(gSystematics.systematic_dEdX_virtual) * content);
                    }
                }

                if(gSystematics.systematic_brem_virtual != 0.0)
                {
                    if(gSystematics.systematic_brem_virtual > 0.0)
                    {
                        Int_t bin_x{gSystematics.h_systematic_brem_virtual_h->GetXaxis()->FindBin(el_energy_0)};
                        Int_t bin_y{gSystematics.h_systematic_brem_virtual_h->GetYaxis()->FindBin(el_energy_1)};
                        Double_t content{gSystematics.h_systematic_brem_virtual_h->GetBinContent(bin_x, bin_y)};
                        weight = weight * (gSystematics.systematic_brem_virtual * content);
                    }
                    else if(gSystematics.systematic_brem_virtual < 0.0)
                    {
                        Int_t bin_x{gSystematics.h_systematic_brem_virtual_l->GetXaxis()->FindBin(el_energy_0)};
                        Int_t bin_y{gSystematics.h_systematic_brem_virtual_l->GetYaxis()->FindBin(el_energy_1)};
                        Double_t content{gSystematics.h_systematic_brem_virtual_l->GetBinContent(bin_x, bin_y)};
                        weight = weight * (std::abs(gSystematics.systematic_brem_virtual) * content);
                    }
                }

                if(gSystematics.systematic_foil_thickness_nominal != 0.0)
                {
                    if(gSystematics.systematic_foil_thickness_nominal > 0.0)
                    {
                        Int_t bin_x{gSystematics.h_systematic_foil_thickness_nominal_h->GetXaxis()->FindBin(el_energy_0)};
                        Int_t bin_y{gSystematics.h_systematic_foil_thickness_nominal_h->GetYaxis()->FindBin(el_energy_1)};
                        Double_t content{gSystematics.h_systematic_foil_thickness_nominal_h->GetBinContent(bin_x, bin_y)};
                        weight = weight * (gSystematics.systematic_foil_thickness_nominal * content);
                    }
                    else if(gSystematics.systematic_foil_thickness_nominal < 0.0)
                    {
                        Int_t bin_x{gSystematics.h_systematic_foil_thickness_nominal_l->GetXaxis()->FindBin(el_energy_0)};
                        Int_t bin_y{gSystematics.h_systematic_foil_thickness_nominal_l->GetYaxis()->FindBin(el_energy_1)};
                        Double_t content{gSystematics.h_systematic_foil_thickness_nominal_l->GetBinContent(bin_x, bin_y)};
                        weight = weight * (std::abs(gSystematics.systematic_foil_thickness_nominal) * content);
                    }
                }

                if(gSystematics.systematic_dEdX_nominal != 0.0)
                {
                    if(gSystematics.systematic_dEdX_nominal > 0.0)
                    {
                        Int_t bin_x{gSystematics.h_systematic_dEdX_nominal_h->GetXaxis()->FindBin(el_energy_0)};
                        Int_t bin_y{gSystematics.h_systematic_dEdX_nominal_h->GetYaxis()->FindBin(el_energy_1)};
                        Double_t content{gSystematics.h_systematic_dEdX_nominal_h->GetBinContent(bin_x, bin_y)};
                        weight = weight * (gSystematics.systematic_dEdX_nominal * content);
                    }
                    else if(gSystematics.systematic_dEdX_nominal < 0.0)
                    {
                        Int_t bin_x{gSystematics.h_systematic_dEdX_nominal_l->GetXaxis()->FindBin(el_energy_0)};
                        Int_t bin_y{gSystematics.h_systematic_dEdX_nominal_l->GetYaxis()->FindBin(el_energy_1)};
                        Double_t content{gSystematics.h_systematic_dEdX_nominal_l->GetBinContent(bin_x, bin_y)};
                        weight = weight * (std::abs(gSystematics.systematic_dEdX_nominal) * content);
                    }
                }

                if(gSystematics.systematic_brem_nominal != 0.0)
                {
                    if(gSystematics.systematic_brem_nominal > 0.0)
                    {
                        Int_t bin_x{gSystematics.h_systematic_brem_nominal_h->GetXaxis()->FindBin(el_energy_0)};
                        Int_t bin_y{gSystematics.h_systematic_brem_nominal_h->GetYaxis()->FindBin(el_energy_1)};
                        Double_t content{gSystematics.h_systematic_brem_nominal_h->GetBinContent(bin_x, bin_y)};
                        weight = weight * (gSystematics.systematic_brem_nominal * content);
                    }
                    else if(gSystematics.systematic_brem_nominal < 0.0)
                    {
                        Int_t bin_x{gSystematics.h_systematic_brem_nominal_l->GetXaxis()->FindBin(el_energy_0)};
                        Int_t bin_y{gSystematics.h_systematic_brem_nominal_l->GetYaxis()->FindBin(el_energy_1)};
                        Double_t content{gSystematics.h_systematic_brem_nominal_l->GetBinContent(bin_x, bin_y)};
                        weight = weight * (std::abs(gSystematics.systematic_brem_nominal) * content);
                    }
                }

                ///////////////////////////////////////////////////////////////
                // standard cuts not implemented in preprocessing
                ///////////////////////////////////////////////////////////////

                if((el_energy_0 < 0.3) || (el_energy_1 < 0.3))
                {
                    continue;
                }

                if(sampleName.Contains("nd150"))
                {
                    weight_total_150Nd += 1.0 * weight;
                    event_pass_count_150Nd += 1.0;
                }

                ///////////////////////////////////////////////////////////////
                // fill histograms
                ///////////////////////////////////////////////////////////////

                if(thePhase == 0)
                {
                    if(paramEnabledP1)
                    {
                        hTotalE_P1_tmp->Fill(el_energy_0 + el_energy_1, 1.0 * weight);
                        hSingleEnergy_P1_tmp->Fill(el_energy_1, 1.0 * weight);
                        hSingleEnergy_P1_tmp->Fill(el_energy_0, 1.0 * weight);
                        hHighEnergy_P1_tmp->Fill(el_energy_0, 1.0 * weight);
                        hLowEnergy_P1_tmp->Fill(el_energy_1, 1.0 * weight);
                        hHighLowEnergy_P1_tmp->Fill(el_energy_1, el_energy_0, 1.0 * weight);
                        hEnergySum_P1_tmp->Fill(el_energy_0 + el_energy_1, 1.0 * weight);
                        hEnergyDiff_P1_tmp->Fill(el_energy_0 - el_energy_1, 1.0 * weight);
                    }
                }
                else if(thePhase == 1)
                {
                    if(paramEnabledP2)
                    {
                        hTotalE_P2_tmp->Fill(el_energy_0 + el_energy_1, 1.0 * weight);
                        hSingleEnergy_P2_tmp->Fill(el_energy_0, 1.0 * weight);
                        hSingleEnergy_P2_tmp->Fill(el_energy_1, 1.0 * weight);
                        hHighEnergy_P2_tmp->Fill(el_energy_0, 1.0 * weight);
                        hLowEnergy_P2_tmp->Fill(el_energy_1, 1.0 * weight);
                        hHighLowEnergy_P2_tmp->Fill(el_energy_1, el_energy_0, 1.0 * weight);
                        hEnergySum_P2_tmp->Fill(el_energy_0 + el_energy_1, 1.0 * weight);
                        hEnergyDiff_P2_tmp->Fill(el_energy_0 - el_energy_1, 1.0 * weight);
                    }
                }
                else
                {
                    std::cout << "Error: thePhase=" << thePhase << " invalid value" << std::endl;
                }


            }

            /*
            if(sampleName.Contains("nd150"))
            {
                for(Int_t i = 1; i <= hSingleEnergy_P2_tmp->GetNbinsX(); ++ i)
                {
                    std::cout << "bin " << i << " -> " << hSingleEnergy_P2_tmp->GetBinContent(i) << std::endl;
                }
            }
            */

            inputFile->Close();
            if(debuglevel >= 8)
            {
                std::cout << "read tree done" << std::endl;
            }
            //std::cout << hSingleEnergy_P2_tmp->GetName() << " -> " << hSingleEnergy_P2_tmp->GetBinContent(10) << std::endl;

            /* 
            if(REWEIGHT_ADD_RANDOM_VARIANCE_150ND == 1)
            {
                if(sampleName.Contains("nd150"))
                {
                    std::cout << "ADD VARIANCE TO 150ND MC" << std::endl;

                    TRandom3 rng;
                    for(Int_t bin_ix = 1; bin_ix <= hTotalE_P1_tmp->GetNbinsX(); ++ bin_ix)
                    {
                        Double_t content = hTotalE_P1_tmp->GetBinContent(bin_ix);
                        if(content != 0.0)
                        {
                            Double_t error = std::sqrt(content);
                            Double_t rgauss = rng.Gaus(0.0, 1.0);
                            Double_t content_out = content + 5.0 * error * rgauss;
                            hTotalE_P1_tmp->SetBinContent(bin_ix, content_out);
                        }
                    }

                    for(Int_t bin_ix = 1; bin_ix <= hTotalE_P1_tmp->GetNbinsX(); ++ bin_ix)
                    {
                        Double_t content = hTotalE_P2_tmp->GetBinContent(bin_ix);
                        if(content != 0.0)
                        {
                            Double_t error = std::sqrt(content);
                            Double_t rgauss = rng.Gaus(0.0, 1.0);
                            Double_t content_out = content + 5.0 * error * rgauss;
                            hTotalE_P2_tmp->SetBinContent(bin_ix, content_out);
                        }
                    }

                    for(Int_t bin_ix = 1; bin_ix <= hSingleEnergy_P1_tmp->GetNbinsX(); ++ bin_ix)
                    {
                        Double_t content = hSingleEnergy_P1_tmp->GetBinContent(bin_ix);
                        if(content != 0.0)
                        {
                            Double_t error = std::sqrt(content);
                            Double_t rgauss = rng.Gaus(0.0, 1.0);
                            Double_t content_out = content + 5.0 * error * rgauss;
                            hSingleEnergy_P1_tmp->SetBinContent(bin_ix, content_out);
                        }
                    }

                    for(Int_t bin_ix = 1; bin_ix <= hSingleEnergy_P1_tmp->GetNbinsX(); ++ bin_ix)
                    {
                        Double_t content = hSingleEnergy_P2_tmp->GetBinContent(bin_ix);
                        if(content != 0.0)
                        {
                            Double_t error = std::sqrt(content);
                            Double_t rgauss = rng.Gaus(0.0, 1.0);
                            Double_t content_out = content + 5.0 * error * rgauss;
                            hSingleEnergy_P2_tmp->SetBinContent(bin_ix, content_out);
                        }
                    }
                }
            }
            else
            {
                if(sampleName.Contains("nd150"))
                {
                    std::cout << "DO NOT ADD VARIANCE TO 150ND MC" << std::endl;
                }
            }
            */



            // tests
            if(sampleName.Contains("nd150"))
            {
                const double scale_factor = ND150_FAKEDATA_SCALE_FACTOR;

                // P1
                if(paramEnabledP1)
                {
                    hTotalE_P1_tmp->Scale(scale_factor);
                    hSingleEnergy_P1_tmp->Scale(scale_factor);
                    hHighEnergy_P1_tmp->Scale(scale_factor);
                    hLowEnergy_P1_tmp->Scale(scale_factor);
                    hHighLowEnergy_P1_tmp->Scale(scale_factor);
                    hEnergySum_P1_tmp->Scale(scale_factor);
                    hEnergyDiff_P1_tmp->Scale(scale_factor);
                }
                // P2
                if(paramEnabledP2)
                {
                    hTotalE_P2_tmp->Scale(scale_factor);
                    hSingleEnergy_P2_tmp->Scale(scale_factor);
                    hHighEnergy_P2_tmp->Scale(scale_factor);
                    hLowEnergy_P2_tmp->Scale(scale_factor);
                    hHighLowEnergy_P2_tmp->Scale(scale_factor);
                    hEnergySum_P2_tmp->Scale(scale_factor);
                    hEnergyDiff_P2_tmp->Scale(scale_factor);
                }
            }





            // scale
            std::ifstream inFile;
            const TString typedir = thePath;
            const TString sample_file_name = sampleName;
            inFile.open(filePath + typedir + sample_file_name + "/JobSummary.txt");
            //std::cout << "file: " << filePath + typedir + sample_file_name + "/JobSummary.txt" << std::endl;
            if(!inFile.is_open())
            {
                std::cout << "file not open" << std::endl;
                std::cin.get();
            }
            std::string dummy;
            double sampleNGenMC;
            inFile >> dummy >> sampleNGenMC;
            inFile.close();

            const double TotalTime = 167629292.;
            
            //std::cout << "sampleNGenMC=" << sampleNGenMC << std::endl;

            //std::cout << hSingleEnergy_P2_tmp->GetName() << " -> scale=" << TotalTime / sampleNGenMC << std::endl;

            // note: was below ->Scale
            if(debuglevel >= 15) // 5
            {
                if(paramEnabledP2)
                {
                    std::cout << hTotalE_P2_tmp->GetName() << " : bin(7) -> " << hTotalE_P2_tmp->GetBinContent(7) << std::endl;
                    std::cout << hTotalE_P2_tmp->GetName() << " : bin(8) -> " << hTotalE_P2_tmp->GetBinContent(8) << std::endl;
                    std::cout << hTotalE_P2_tmp->GetName() << " : bin(9) -> " << hTotalE_P2_tmp->GetBinContent(9) << std::endl;
                    std::cout << hTotalE_P2_tmp->GetName() << " : bin(10) -> " << hTotalE_P2_tmp->GetBinContent(10) << std::endl;
                    std::cout << hTotalE_P2_tmp->GetName() << " : bin(11) -> " << hTotalE_P2_tmp->GetBinContent(11) << std::endl;
                    std::cout << hTotalE_P2_tmp->GetName() << " : bin(12) -> " << hTotalE_P2_tmp->GetBinContent(12) << std::endl;
                    std::cout << hTotalE_P2_tmp->GetName() << " : bin(13) -> " << hTotalE_P2_tmp->GetBinContent(13) << std::endl;
                
                    /*
                    std::cout << hSingleEnergy_P2_tmp->GetName() << " : bin(7) -> " << hSingleEnergy_P2_tmp->GetBinContent(7) << std::endl;
                    std::cout << hSingleEnergy_P2_tmp->GetName() << " : bin(8) -> " << hSingleEnergy_P2_tmp->GetBinContent(8) << std::endl;
                    std::cout << hSingleEnergy_P2_tmp->GetName() << " : bin(9) -> " << hSingleEnergy_P2_tmp->GetBinContent(9) << std::endl;
                    std::cout << hSingleEnergy_P2_tmp->GetName() << " : bin(10) -> " << hSingleEnergy_P2_tmp->GetBinContent(10) << std::endl;
                    std::cout << hSingleEnergy_P2_tmp->GetName() << " : bin(11) -> " << hSingleEnergy_P2_tmp->GetBinContent(11) << std::endl;
                    std::cout << hSingleEnergy_P2_tmp->GetName() << " : bin(12) -> " << hSingleEnergy_P2_tmp->GetBinContent(12) << std::endl;
                    std::cout << hSingleEnergy_P2_tmp->GetName() << " : bin(13) -> " << hSingleEnergy_P2_tmp->GetBinContent(13) << std::endl;
                    */
                }
            }

            // P1
            if(paramEnabledP1)
            {
                hTotalE_P1_tmp->Scale(TotalTime / sampleNGenMC);
                hSingleEnergy_P1_tmp->Scale(TotalTime / sampleNGenMC);
                hHighEnergy_P1_tmp->Scale(TotalTime / sampleNGenMC);
                hLowEnergy_P1_tmp->Scale(TotalTime / sampleNGenMC);
                hHighLowEnergy_P1_tmp->Scale(TotalTime / sampleNGenMC);
                hEnergySum_P1_tmp->Scale(TotalTime / sampleNGenMC);
                hEnergyDiff_P1_tmp->Scale(TotalTime / sampleNGenMC);
            }
            // P2
            if(paramEnabledP2)
            {
                hTotalE_P2_tmp->Scale(TotalTime / sampleNGenMC);
                hSingleEnergy_P2_tmp->Scale(TotalTime / sampleNGenMC);
                hHighEnergy_P2_tmp->Scale(TotalTime / sampleNGenMC);
                hLowEnergy_P2_tmp->Scale(TotalTime / sampleNGenMC);
                hHighLowEnergy_P2_tmp->Scale(TotalTime / sampleNGenMC);
                hEnergySum_P2_tmp->Scale(TotalTime / sampleNGenMC);
                hEnergyDiff_P2_tmp->Scale(TotalTime / sampleNGenMC);
            }


            //std::cout << "paramInitValue=" << paramInitValue << std::endl;

            if(sampleName.Contains("tl208"))
            {
                const double scale_factor = 0.36;

                // P1
            if(paramEnabledP1)
            {
                hTotalE_P1_tmp->Scale(scale_factor);
                hSingleEnergy_P1_tmp->Scale(scale_factor);
                hHighEnergy_P1_tmp->Scale(scale_factor);
                hLowEnergy_P1_tmp->Scale(scale_factor);
                hHighLowEnergy_P1_tmp->Scale(scale_factor);
                hEnergySum_P1_tmp->Scale(scale_factor);
                hEnergyDiff_P1_tmp->Scale(scale_factor);
            }
                // P2
            if(paramEnabledP2)
            {
                hTotalE_P2_tmp->Scale(scale_factor);
                hSingleEnergy_P2_tmp->Scale(scale_factor);
                hHighEnergy_P2_tmp->Scale(scale_factor);
                hLowEnergy_P2_tmp->Scale(scale_factor);
                hHighLowEnergy_P2_tmp->Scale(scale_factor);
                hEnergySum_P2_tmp->Scale(scale_factor);
                hEnergyDiff_P2_tmp->Scale(scale_factor);
            }
            }

            // realistic scaling for amplitude
            // debatable whether this is a good thing to include
            // makes the MPS have the same y axis range
            /*
            if(sampleName.Contains("nd150"))
            {
                const double scale_factor = event_pass_count_150Nd / weight_total_150Nd;

                // P1
                hTotalE_P1_tmp->Scale(scale_factor);
                hSingleEnergy_P1_tmp->Scale(scale_factor);
                hHighEnergy_P1_tmp->Scale(scale_factor);
                hLowEnergy_P1_tmp->Scale(scale_factor);
                hHighLowEnergy_P1_tmp->Scale(scale_factor);
                hEnergySum_P1_tmp->Scale(scale_factor);
                hEnergyDiff_P1_tmp->Scale(scale_factor);
                // P2
                hTotalE_P2_tmp->Scale(scale_factor);
                hSingleEnergy_P2_tmp->Scale(scale_factor);
                hHighEnergy_P2_tmp->Scale(scale_factor);
                hLowEnergy_P2_tmp->Scale(scale_factor);
                hHighLowEnergy_P2_tmp->Scale(scale_factor);
                hEnergySum_P2_tmp->Scale(scale_factor);
                hEnergyDiff_P2_tmp->Scale(scale_factor);
            }
            */

            // P1
            if(paramEnabledP1)
            {
            hTotalE_P1_tmp->Scale(paramInitValue);
            hSingleEnergy_P1_tmp->Scale(paramInitValue);
            hHighEnergy_P1_tmp->Scale(paramInitValue);
            hLowEnergy_P1_tmp->Scale(paramInitValue);
            hHighLowEnergy_P1_tmp->Scale(paramInitValue);
            hEnergySum_P1_tmp->Scale(paramInitValue);
            hEnergyDiff_P1_tmp->Scale(paramInitValue);
            }
            // P2
            if(paramEnabledP2)
            {
            hTotalE_P2_tmp->Scale(paramInitValue);
            hSingleEnergy_P2_tmp->Scale(paramInitValue);
            hHighEnergy_P2_tmp->Scale(paramInitValue);
            hLowEnergy_P2_tmp->Scale(paramInitValue);
            hHighLowEnergy_P2_tmp->Scale(paramInitValue);
            hEnergySum_P2_tmp->Scale(paramInitValue);
            hEnergyDiff_P2_tmp->Scale(paramInitValue);
            }

            /*
            if(paramEnabledP2)
            {
                std::cout << "hSingleEnergy_P2->GetName()=" << hSingleEnergy_P2_tmp->GetName() << " scale_factor_P2=" << paramInitValue << std::endl;
                for(Int_t i = 1; i < hSingleEnergy_P2_tmp->GetNbinsX(); ++ i)
                {
                    if(hSingleEnergy_P2_tmp->GetBinContent(i) != 0.0)
                    {
                        std::cout << "hSingleEnergy_P2->GetBinContents(" << i << ")=" << hSingleEnergy_P2_tmp->GetBinContent(i) << std::endl;
                    }
                }
            }
            */

            
            //if(TString(hSingleEnergy_P2_tmp->GetName()).Contains("tl208"))
            //{
            //    std::cout << __func__ << " GetName()->" << hSingleEnergy_P2_tmp->GetName() << ", bin 13: " << hSingleEnergy_P2_tmp->GetBinContent(13) << std::endl;
            //}

            /*
            std::cout << "name: " << hSingleEnergy_P2_tmp->GetName() << std::endl;
            std::cout << "integral: " << hSingleEnergy_P2_tmp->Integral() << std::endl;
            if(hSingleEnergy_P2_tmp->Integral() > 0.0)
            {
                TCanvas *c = new TCanvas("ctmpdeleteme", "ctmpdeleteme");
                hSingleEnergy_P2_tmp->Draw("hist");
                hSingleEnergy_P1_tmp->Draw("histsame");
                c->SaveAs("ctmpdeleteme.png");
            }
            std::cin.get();
            */


            ///////////////////////////////////////////////////////////////////
            // stack
            ///////////////////////////////////////////////////////////////////

            // P1
            if(paramEnabledP1 == true)
            {
                hTotalE_output_P1->Add(hTotalE_P1_tmp);
                hSingleEnergy_output_P1->Add(hSingleEnergy_P1_tmp);
                hHighEnergy_output_P1->Add(hHighEnergy_P1_tmp);
                hLowEnergy_output_P1->Add(hLowEnergy_P1_tmp);
                hHighLowEnergy_output_P1->Add(hHighLowEnergy_P1_tmp);
                hEnergySum_output_P1->Add(hEnergySum_P1_tmp);
                hEnergyDiff_output_P1->Add(hEnergyDiff_P1_tmp);
            }
            // P2
            if(paramEnabledP2 == true)
            {
                hTotalE_output_P2->Add(hTotalE_P2_tmp);
                hSingleEnergy_output_P2->Add(hSingleEnergy_P2_tmp);
                hHighEnergy_output_P2->Add(hHighEnergy_P2_tmp);
                hLowEnergy_output_P2->Add(hLowEnergy_P2_tmp);
                hHighLowEnergy_output_P2->Add(hHighLowEnergy_P2_tmp);
                hEnergySum_output_P2->Add(hEnergySum_P2_tmp);
                hEnergyDiff_output_P2->Add(hEnergyDiff_P2_tmp);

                //std::cout << __func__ << " GetName()->" << hSingleEnergy_P2_tmp->GetName() << ", bin 13: " << hSingleEnergy_P2_tmp->GetBinContent(13) << " scale_factor_P2=" << paramInitValue <<std::endl;
            }
            // NOTE: 2020-09-10:
            // 
            // This section was missing the check on paramEnabledPx
            //
            // There may be other sections of code missing this? (Check
            // draw for example)
            // Checked draw: drawchannel function is ok

            // what to do with gEnabledPhasex?
        }

    }
    

}

#if 0
    // loop over all the parameters
    std::map<int, file_parameter>::iterator it{g_pg.file_params.begin()};
    for(; it != g_pg.file_params.end(); ++ it)
    {
        int paramNumberInt = -1;

        int paramNumber = it->second.paramNumber;
        double paramInitValue = it->second.paramInitValue;
        double paramInitError = it->second.paramInitError;

        // loop over all MC names
        std::vector<std::string>::iterator mc_name_it{it->second.MCNameList.begin()};
        for(; mc_name_it != it->second.MCNameList.end(); ++ mc_name_it)
        {
            std::string mc_name = *mc_name_it;
            std::string histname = std::string(channel_histname_1D[channel]);
            std::string search_object_P1 = histname + mc_name + "_P1_fit";
            std::string search_object_P2 = histname + mc_name + "_P2_fit";
            TH1D *tmpHist1D_P1 = nullptr;
            TH1D *tmpHist1D_P2 = nullptr;

            if(debuglevel >= 5)
            {
                std::cout << "search_object_P1=" << search_object_P1
                          << " search_object_P2=" << search_object_P2 << std::endl;
            }

            paramNumberInt = g_pg.ExtToIntParamNumberMap.at(paramNumber);
            if(debuglevel >= 5)
            {
                std::cout << "paramNumber=" << paramNumber << " -> " << paramNumberInt << std::endl;
            }
#endif

#endif // NEWLOGLIKFITTER_REWEIGHT_APPLY_FAKEDATA_H
