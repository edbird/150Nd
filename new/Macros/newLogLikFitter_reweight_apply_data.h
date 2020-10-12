#ifndef NEWLOGLIKEFITTER_REWEIGHT_APPLY_DATA_H
#define NEWLOGLIKEFITTER_REWEIGHT_APPLY_DATA_H



void reweight_apply_data(
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
    TH2D *&hHighLowEnergy_output_P2)
{

    ///////////////////////////////////////////////////////////////////////////
    // this function reloads the data histograms from data in root trees
    // and applies systematics (currently disabled)
    ///////////////////////////////////////////////////////////////////////////


    std::cout << __func__ << std::endl;
    std::cout << "this function should not be called" << std::endl;
    std::cout << "systematics should be applied to fake data only!" << std::endl;
    std::cout << "systematics disabled in this function" << std::endl;
    std::cin.get();


    const int debuglevel = 10;

    TString tinput_filename = "Nd150_2eNg_output_postprocessed_small.root";

    if(debuglevel >= 5)
    {
    /*
        if(gSystematics.systematic_energy_offset != gSystematics.systematic_energy_offset_last)
        {
            std::cout << "gSystematics.systematic_energy_offset=" << gSystematics.systematic_energy_offset << std::endl;
            std::cin.get();
            gSystematics.systematic_energy_offset_last = gSystematics.systematic_energy_offset;
        }
        if(gSystematics.systematic_energy_scale != gSystematics.systematic_energy_scale_last)
        {
            std::cout << "gSystematics.systematic_energy_scale=" << gSystematics.systematic_energy_scale << std::endl;
            std::cin.get();
            gSystematics.systematic_energy_scale_last = gSystematics.systematic_energy_scale;
        }
    */
    }


    const TString sampleName = DataFile;
    const TString name_append = "";


    // TODO: should these have "reweight" at the end? is this not a (minor)
    // problem in that these histograms do not exist on the first call
    // to this function?
    // is the P1/P2 in the right place?

    if(debuglevel >= 3)
    {
        std::cout << "search and remove P1 histograms" << std::endl;
    }
 
    // this is required to remove the objects from gROOT memory
    deleter_helper("hTotalE_" + sampleName + name_append + "_P1");
    deleter_helper("hSingleEnergy_" + sampleName + name_append + "_P1");
    deleter_helper("hHighEnergy_" + sampleName + name_append + "_P1"); 
    deleter_helper("hLowEnergy_" + sampleName + name_append + "_P1"); 
    deleter_helper("hEnergySum_" + sampleName + name_append + "_P1"); 
    deleter_helper("hEnergyDiff_" + sampleName + name_append + "_P1"); 
    deleter_helper("hHighLowEnergy_" + sampleName + name_append + "_P1");
    deleter_helper("hWeight_" + sampleName + name_append + "_P1"); 


    if(debuglevel >= 3)
    {
        std::cout << "search and remove P2 histograms" << std::endl;
    }

    // this is required to remove the objects from gROOT memory
    deleter_helper("hTotalE_" + sampleName + name_append + "_P2");
    deleter_helper("hSingleEnergy_" + sampleName + name_append + "_P2");
    deleter_helper("hHighEnergy_" + sampleName + name_append + "_P2"); 
    deleter_helper("hLowEnergy_" + sampleName + name_append + "_P2"); 
    deleter_helper("hEnergySum_" + sampleName + name_append + "_P2"); 
    deleter_helper("hEnergyDiff_" + sampleName + name_append + "_P2"); 
    deleter_helper("hHighLowEnergy_" + sampleName + name_append + "_P2");
    deleter_helper("hWeight_" + sampleName + name_append + "_P2"); 




    if(debuglevel >= 3)
    {
        std::cout << "create P1 histograms" << std::endl;
    }
    //std::cin.get();


    // TODO: this does not work, need to re-Fill histogram using file
    //output = input->Clone(input->GetName() + "_reweight");
    //std::cout << "Total new" << std::endl;
    hTotalE_output_P1 = new TH1D("hTotalE_" + sampleName + name_append + "_P1",
                       //"Phase " + Phase + " " + sampleName + name_append + " total energy; #SigmaE_{e} (MeV)",
                       TString("Phase ") + "P1" + " " + sampleName + name_append + " total energy; Total Energy #SigmaE_{e} (MeV)",
                       50, 0.0, 5.0);
                       // TODO: changed from 4

    //std::cout << "Single new" << std::endl;
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


    if(debuglevel >= 3)
    {
        std::cout << "create P2 histograms" << std::endl;
    }

    hTotalE_output_P2 = new TH1D("hTotalE_" + sampleName + name_append + "_P2",
                       //"Phase " + Phase + " " + sampleName + name_append + " total energy; #SigmaE_{e} (MeV)",
                       TString("Phase ") + "P2" + " " + sampleName + name_append + " total energy; Total Energy #SigmaE_{e} (MeV)",
                       50, 0.0, 5.0);
                       // TODO: changed from 4

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


    //std::cout << "end of new" << std::endl;

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



    if(debuglevel >= 5)
    {
        std::cout << "filling... (data)" << std::endl;
    }

    TString filePath = "/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
    TString thePath = "betabeta/";
    //TString sampleName = DataName;
    TFile *finput = new TFile(filePath + thePath + sampleName + "/" + tinput_filename);
    if(!finput->IsOpen())
    {
        std::cout << "failed to open" << std::endl;
    }
    TTree *tinput = (TTree*)finput->Get("Nd150_2eNg/Nd150_2eNg");

    Int_t run;
    double electronEnergy[2];
    double trueElectronEnergy[2];

    tinput->SetBranchAddress("Run"                  , &run);  
    tinput->SetBranchAddress("electronEnergy"       , electronEnergy);
    tinput->SetBranchAddress("trueElectronEnergy"   , trueElectronEnergy);


    for(Long64_t ix{0}; ix < tinput->GetEntries(); ++ ix)
    {
        std::cout << "ix=" << ix << std::endl;

        tinput->GetEntry(ix);

        // Phase cut moved to be final cut
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

        Double_t el_energy_0{electronEnergy[highE_index]};
        Double_t el_energy_1{electronEnergy[lowE_index]};

        // linear energy offset systematic
        //el_energy_0 = el_energy_0 + gSystematics.systematic_energy_offset;
        //el_energy_1 = el_energy_1 + gSystematics.systematic_energy_offset;

        // now file contains energy down to about ~ 0.2 MeV
        //if((electronEnergy[0] < 0.3) || (electronEnergy[1] < 0.3)) continue;

        if(el_energy_0 < 0.3)
        {
            continue;
        }
        if(el_energy_1 < 0.3)
        {
            continue;
        }

        const double weight = 1.0;

        if(thePhase == 0)
        {
            hTotalE_output_P1->Fill(el_energy_0 + el_energy_1, 1.0 * weight);
            hSingleEnergy_output_P1->Fill(el_energy_1, 1.0 * weight);
            hSingleEnergy_output_P1->Fill(el_energy_0, 1.0 * weight);
            hHighEnergy_output_P1->Fill(el_energy_0, 1.0 * weight);
            hLowEnergy_output_P1->Fill(el_energy_1, 1.0 * weight);
            hHighLowEnergy_output_P1->Fill(el_energy_1, el_energy_0, 1.0 * weight);
            hEnergySum_output_P1->Fill(el_energy_0 + el_energy_1, 1.0 * weight);
            hEnergyDiff_output_P1->Fill(el_energy_0 - el_energy_1, 1.0 * weight);
        }
        else if(thePhase == 1)
        {
            hTotalE_output_P2->Fill(el_energy_0 + el_energy_1, 1.0 * weight);
            hSingleEnergy_output_P2->Fill(el_energy_1, 1.0 * weight);
            hSingleEnergy_output_P2->Fill(el_energy_0, 1.0 * weight);
            hHighEnergy_output_P2->Fill(el_energy_0, 1.0 * weight);
            hLowEnergy_output_P2->Fill(el_energy_1, 1.0 * weight);
            hHighLowEnergy_output_P2->Fill(el_energy_1, el_energy_0, 1.0 * weight);
            hEnergySum_output_P2->Fill(el_energy_0 + el_energy_1, 1.0 * weight);
            hEnergyDiff_output_P2->Fill(el_energy_0 - el_energy_1, 1.0 * weight);
        }
        else
        {
            std::cout << "Error: thePhase=" << thePhase << " invalid value" << std::endl;
        }

    }

    finput->Close();

    std::cout << "finished building the data" << std::endl;
    std::cout << "P2 channel1 integral: " << hSingleEnergy_output_P2->Integral() << std::endl;
    std::cout << "seo=" << gSystematics.systematic_energy_offset << std::endl;
    std::cout << "sem=" << gSystematics.systematic_energy_scale << std::endl;
    std::cout << "sef=" << gSystematics.systematic_efficiency << std::endl;
    std::cout << "sen=" << gSystematics.systematic_enrichment << std::endl;
    std::cout << "seos=" << gSystematics.systematic_energy_offsetsmall << std::endl;
    std::cin.get();


    // TODO: may need to apply some visual modifcations here
    // for example setting line width etc
    // done differently for data to MC, check fit_2e.C code





}
    



#endif // NEWLOGLIKEFITTER_REWEIGHT_APPLY_DATA_H
