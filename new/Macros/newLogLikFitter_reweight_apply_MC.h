#ifndef NEWLOGLIKFITTER_REWEIGHT_APPLY_H
#define NEWLOGLIKFITTER_REWEIGHT_APPLY_H



void deleter_helper(TString hname)
{
    int debuglevel = 1;

    TObject *hist = gROOT->FindObject(hname);
    if(hist == nullptr)
    {
        if(debuglevel >= 3)
        {
            std::cout << hname << " nullptr!" << std::endl;
        }
    }
    else
    {
        if(debuglevel >= 3)
        {
            std::cout << "delete " << hist->GetName() << std::endl;
        }
        delete hist;
    }
}






// xi31 reweighting function for 150Nd spectra
//void reweight_apply(TH2D *output, const TH1D *const input,
void reweight_apply_MC(
    TH1D *&hTotalE_output_P1,
    TH1D *&hSingleEnergy_output_P1,
    TH1D *&hHighEnergy_output_P1,
    TH1D *&hLowEnergy_output_P1,
    TH1D *&hEnergySum_output_P1,
    TH1D *&hEnergyDiff_output_P1,
    TH2D *&hHighLowEnergy_output_P1,
    TH1D *&hWeight_output_P1,
    TH1D *&hTotalE_output_P2,
    TH1D *&hSingleEnergy_output_P2,
    TH1D *&hHighEnergy_output_P2,
    TH1D *&hLowEnergy_output_P2,
    TH1D *&hEnergySum_output_P2,
    TH1D *&hEnergyDiff_output_P2,
    TH2D *&hHighLowEnergy_output_P2,
    TH1D *&hWeight_output_P2,
//    const std::string& tinput_filename,
    const Double_t xi_31,
    const Double_t xi_31_baseline,
    const TH2D* const h_nEqNull,
    const TH2D* const h_nEqTwo,
    const Double_t psiN0,
    const Double_t psiN2,
    const Double_t bb_Q)
{

    //std::cout << "__func__ : " << "(xi_31=" << xi_31 << ")" << std::endl;
    //std::cout << "are you happy about this?" << std::endl;
    //std::cin.get();


    ///////////////////////////////////////////////////////////////////////////
    // this function reloads 150Nd MC from file, reweight applied depending on
    // xi_31 parameter
    ///////////////////////////////////////////////////////////////////////////

    //TString tinput_filename = "Nd150_2eNg_output_truth_postprocessed_small.root";
    TString tinput_filename = "Nd150_2eNg_output_postprocessed_small.root";

    int debuglevel = 1;

    if(debuglevel >= 2)
    {
        std::cout << "reweight_apply: xi_31=" << xi_31 << " xi_31_baseline=" << xi_31_baseline << std::endl;
    }

    /*
    if(debuglevel >= 5)
    {
        if(systematic_energy_offset != systematic_energy_offset_last)
        {
            std::cout << "systematic_energy_offset=" << systematic_energy_offset << std::endl;
            std::cin.get();
            systematic_energy_offset_last = systematic_energy_offset;
        }
    }
    */

    //std::cout << "reweight_apply" << std::endl;

    //const TString sampleName = "nd150_rot_2n2b_m4";
    const TString sampleName = Nd150Files[0];
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
    deleter_helper("hTotalE_" + sampleName + name_append + "_P1_fit");
    deleter_helper("hSingleEnergy_" + sampleName + name_append + "_P1_fit");
    deleter_helper("hHighEnergy_" + sampleName + name_append + "_P1_fit"); 
    deleter_helper("hLowEnergy_" + sampleName + name_append + "_P1_fit"); 
    deleter_helper("hEnergySum_" + sampleName + name_append + "_P1_fit"); 
    deleter_helper("hEnergyDiff_" + sampleName + name_append + "_P1_fit"); 
    deleter_helper("hHighLowEnergy_" + sampleName + name_append + "_P1_fit");
    deleter_helper("hWeight_" + sampleName + name_append + "_P1_fit"); 


    if(debuglevel >= 3)
    {
        std::cout << "search and remove P2 histograms" << std::endl;
    }

    // this is required to remove the objects from gROOT memory
    deleter_helper("hTotalE_" + sampleName + name_append + "_P2_fit");
    deleter_helper("hSingleEnergy_" + sampleName + name_append + "_P2_fit");
    deleter_helper("hHighEnergy_" + sampleName + name_append + "_P2_fit"); 
    deleter_helper("hLowEnergy_" + sampleName + name_append + "_P2_fit"); 
    deleter_helper("hEnergySum_" + sampleName + name_append + "_P2_fit"); 
    deleter_helper("hEnergyDiff_" + sampleName + name_append + "_P2_fit"); 
    deleter_helper("hHighLowEnergy_" + sampleName + name_append + "_P2_fit");
    deleter_helper("hWeight_" + sampleName + name_append + "_P2_fit"); 

    if(debuglevel >= 4)
    {
        std::cout << "gROOT->FindObject test:" << std::endl;
        std::cout << "gROOT->FindObject(\"hTotalE_nd150_rot_2n2b_m4_P2_fit\") -> " << gROOT->FindObject("hTotalE_nd150_rot_2n2b_m4_P2_fit") << std::endl;
        std::cout << "gROOT->FindObject(\"hHighLowEnergy_nd150_rot_2n2b_m4_P2_fit\") -> " << gROOT->FindObject("hHighLowEnergy_nd150_rot_2n2b_m4_P2_fit") << std::endl;
        delete (TH1D*)gROOT->FindObject("hTotalE_nd150_rot_2n2b_m4_P2_fit");
        delete (TH2D*)gROOT->FindObject("hHighLowEnergy_nd150_rot_2n2b_m4_P2_fit");
        std::cout << "gROOT->FindObject(\"hTotalE_nd150_rot_2n2b_m4_P2_fit\") -> " << gROOT->FindObject("hTotalE_nd150_rot_2n2b_m4_P2_fit") << std::endl;
        std::cout << "gROOT->FindObject(\"hHighLowEnergy_nd150_rot_2n2b_m4_P2_fit\") -> " << gROOT->FindObject("hHighLowEnergy_nd150_rot_2n2b_m4_P2_fit") << std::endl;
        std::cout << "end of test" << std::endl;
        //std::cin.get();
    }

    //std::cout << "gROOT" << std::endl;
    //std::cout << gROOT->GetListOfClasses() << std::endl;
    //gROOT->GetListOfClasses()->Print();
    //gDirectory->ls();

    if(debuglevel >= 3)
    {
        std::cout << "create P1 histograms" << std::endl;
    }
    //std::cin.get();


    // TODO: this does not work, need to re-Fill histogram using file
    //output = input->Clone(input->GetName() + "_reweight");
    //std::cout << "Total new" << std::endl;
    hTotalE_output_P1 = new TH1D("hTotalE_" + sampleName + name_append + "_P1_fit",
                       //"Phase " + Phase + " " + sampleName + name_append + " total energy; #SigmaE_{e} (MeV)",
                       TString("Phase ") + "P1" + " " + sampleName + name_append + " total energy; Total Energy #SigmaE_{e} (MeV)",
                       50, 0.0, 5.0);
                       // TODO: changed from 4

    //std::cout << "Single new" << std::endl;
    hSingleEnergy_output_P1    = new TH1D("hSingleEnergy_" + sampleName + name_append + "_P1_fit",
                                TString("Phase ") + "P1" + " " + sampleName  + name_append + " Single Energy",
                                50, 0.0, 5.0);
    
    //std::cout << "High new" << std::endl;
    hHighEnergy_output_P1     = new TH1D("hHighEnergy_" + sampleName + name_append + "_P1_fit",
                                TString("Phase ") + "P1" + " " + sampleName + name_append + " High Energy; Energy (MeV)",
                                50, 0.0, 5.0);

    //std::cout << "Low new" << std::endl;
    hLowEnergy_output_P1     = new TH1D("hLowEnergy_" + sampleName + name_append + "_P1_fit",
                                TString("Phase ") + "P1" + " " + sampleName + name_append + " Low Energy",
                                50, 0.0, 5.0);

    //std::cout << "Sum new" << std::endl;
    hEnergySum_output_P1     = new TH1D("hEnergySum_" + sampleName + name_append + "_P1_fit",
                                TString("Phase ") + "P1" + " " + sampleName + name_append + " Low Energy",
                                50, 0.0, 5.0);

    //std::cout << "Diff new" << std::endl;
    hEnergyDiff_output_P1     = new TH1D("hEnergyDiff_" + sampleName + name_append + "_P1_fit",
                                TString("Phase ") + "P1" + " " + sampleName + name_append + " Low Energy",
                                50, 0.0, 5.0);

    //std::cout << "HighLow new" << std::endl;
    hHighLowEnergy_output_P1     = new TH2D("hHighLowEnergy_" + sampleName + name_append + "_P1_fit",
                                TString("Phase ") + "P1" + " " + sampleName + name_append + ";Low Energy Electron Energy (MeV);High Energy Electron Energy (MeV)",
                                50, 0.0, 5.0, 50, 0.0, 5.0);

    //std::cout << "Weight new" << std::endl;
    hWeight_output_P1       = new TH1D("hWeight_" + sampleName + name_append + "_P1_fit",
                                TString("Phase ") + "P1" + " " + sampleName + name_append + ";Weight",
                                //50, -2.0, 4.0);
                                50, 0.0, 0.0);

    if(debuglevel >= 3)
    {
        std::cout << "create P2 histograms" << std::endl;
    }

    hTotalE_output_P2 = new TH1D("hTotalE_" + sampleName + name_append + "_P2_fit",
                       //"Phase " + Phase + " " + sampleName + name_append + " total energy; #SigmaE_{e} (MeV)",
                       TString("Phase ") + "P2" + " " + sampleName + name_append + " total energy; Total Energy #SigmaE_{e} (MeV)",
                       50, 0.0, 5.0);
                       // TODO: changed from 4

    //std::cout << "Single new" << std::endl;
    hSingleEnergy_output_P2    = new TH1D("hSingleEnergy_" + sampleName + name_append + "_P2_fit",
                                TString("Phase ") + "P2" + " " + sampleName  + name_append + " Single Energy",
                                50, 0.0, 5.0);
    
    //std::cout << "High new" << std::endl;
    hHighEnergy_output_P2     = new TH1D("hHighEnergy_" + sampleName + name_append + "_P2_fit",
                                TString("Phase ") + "P2" + " " + sampleName + name_append + " High Energy; Energy (MeV)",
                                50, 0.0, 5.0);

    //std::cout << "Low new" << std::endl;
    hLowEnergy_output_P2     = new TH1D("hLowEnergy_" + sampleName + name_append + "_P2_fit",
                                TString("Phase ") + "P2" + " " + sampleName + name_append + " Low Energy",
                                50, 0.0, 5.0);

    //std::cout << "Sum new" << std::endl;
    hEnergySum_output_P2     = new TH1D("hEnergySum_" + sampleName + name_append + "_P2_fit",
                                TString("Phase ") + "P2" + " " + sampleName + name_append + " Low Energy",
                                50, 0.0, 5.0);

    //std::cout << "Diff new" << std::endl;
    hEnergyDiff_output_P2     = new TH1D("hEnergyDiff_" + sampleName + name_append + "_P2_fit",
                                TString("Phase ") + "P2" + " " + sampleName + name_append + " Low Energy",
                                50, 0.0, 5.0);

    //std::cout << "HighLow new" << std::endl;
    hHighLowEnergy_output_P2     = new TH2D("hHighLowEnergy_" + sampleName + name_append + "_P2_fit",
                                TString("Phase ") + "P2" + " " + sampleName + name_append + ";Low Energy Electron Energy (MeV);High Energy Electron Energy (MeV)",
                                50, 0.0, 5.0, 50, 0.0, 5.0);

    //std::cout << "Weight new" << std::endl;
    hWeight_output_P2       = new TH1D("hWeight_" + sampleName + name_append + "_P2_fit",
                                TString("Phase ") + "P2" + " " + sampleName + name_append + ";Weight",
                                //50, -2.0, 4.0);
                                50, 0.0, 0.0);

    //std::cout << "end of new" << std::endl;

    hTotalE_output_P1->Sumw2();
    hSingleEnergy_output_P1->Sumw2();
    hHighEnergy_output_P1->Sumw2();
    hLowEnergy_output_P1->Sumw2();
    hHighLowEnergy_output_P1->Sumw2();
    hEnergySum_output_P1->Sumw2();
    hEnergyDiff_output_P1->Sumw2();

    hWeight_output_P1->Sumw2();

    hTotalE_output_P2->Sumw2();
    hSingleEnergy_output_P2->Sumw2();
    hHighEnergy_output_P2->Sumw2();
    hLowEnergy_output_P2->Sumw2();
    hHighLowEnergy_output_P2->Sumw2();
    hEnergySum_output_P2->Sumw2();
    hEnergyDiff_output_P2->Sumw2();

    hWeight_output_P2->Sumw2();

    //std::cout << "houtput->GetName() -> " << houtput->GetName() << std::endl;

    //for(Int_t i{1}; i <= output->GetNbinsX(); ++ i)
    //{
    //    Double_t content = output->GetBinContent(i);
    //    Double_t reweight_factor = ReWeight3(T1, T2, xi_31, xi_31_baseline,
    //                                         hNEqNull, hNEqTwo, psiN0, psiN2, "false");
    //    content *= reweight_factor;
    //    output->SetBinContent(i, content);
    //}


//    std::cout << "xi_31=" << xi_31 << ", xi_31_baseline=" << xi_31_baseline << std::endl;

    //const double &epsilon_31_baseline{xi_31_baseline};
    //const double &epsilon_31{xi_31};
    //const Double_t &bb_Q{staticsgroup.Get_bb_Q()};
    //std::cout << "bb_Q=" << bb_Q << std::endl;

    //const TH2D* const h_nEqNull{staticsgroup.Get_h_nEqNull()};
    //const TH2D* const h_nEqTwo{staticsgroup.Get_h_nEqTwo()};
    //const double& psiN0{staticsgroup.Get_psiN0()};
    //const double& psiN2{staticsgroup.Get_psiN2()};

    //TTree *t{staticsgroup.Get_tree()};

    //const Int_t &nElectrons{staticsgroup.Get_nElectrons()};
    //const Double_t &trueT1{staticsgroup.Get_trueT1()};
    //const Double_t &trueT2{staticsgroup.Get_trueT2()};
    //const Double_t* const el_energy_{staticsgroup.Get_el_energy_()};
    //const Double_t &gen_weight{staticsgroup.Get_gen_weight()};
    
    if(debuglevel >= 5)
    {
        std::cout << "filling... (MC)" << std::endl;
    }

    TString filePath = "/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
    TString thePath = "nd150/";
    //TString sampleName = Nd150Names[0];
    TFile *finput = new TFile(filePath + thePath + sampleName + "/" + tinput_filename);
    TTree *tinput = (TTree*)finput->Get("Nd150_2eNg/Nd150_2eNg");

    Int_t run;
    //Int_t nElectrons;
    double electronEnergy[2];
    double trueElectronEnergy[2];

    tinput->SetBranchAddress("Run"                  , &run);  
    //tinput->SetBranchAddress("nElectrons"           , &nElectrons);  
    tinput->SetBranchAddress("electronEnergy"       , electronEnergy);
    tinput->SetBranchAddress("trueElectronEnergy"   , trueElectronEnergy);

    bool doprint = false;

    //const Double_t &trueT1{trueElectronEnergy[0]};
    //const Double_t &trueT2{trueElectronEnergy[1]};

    // TODO: there appears to be a very slight difference in the number of data/MC events
    // compared to what is produced by fit_2e.C
    // this is due to background MC scaling differences (in activities.txt / parameter_names.lst)
    // or is it?
    for(Long64_t ix{0}; ix < tinput->GetEntries(); ++ ix)
    {
//        std::cout << "ix=" << ix << std::endl;

//        if(ix == 10) break;

        tinput->GetEntry(ix);

        // analysis only valid for 2 electron events
        //if(nElectrons != 2) continue;


        // Phase cut moved to be final cut
        // 1869 - 3395 is P1
        // 3396 - 9186 is P2
        //if(mode_flag == 0)
        //{
        bool cut = true;
        Int_t thePhase = -1;
        //if(thePhase == 0)
        //{
        if((1869 <= run) && (run <= 3395))
        {
            cut = false;
            thePhase = 0;
        }
        //}
        //else if(thePhase == 1)
        //{
        if((3396 <= run) && (run <= 9186))
        {
            cut = false;
            thePhase = 1;
        }
        //}
        //else if(nocutonphase == 1)
        //{
        //    cut = false;
        //}
        //else
        //{
        //    cut = true;
        //}
        // accept P1 and P2 for testing purposes
        //cut = false;
        if(cut == true) continue;
        //}

        // now file contains energy down to about ~ 0.2 MeV
        //if((electronEnergy[0] < 0.3) || (electronEnergy[1] < 0.3)) continue;

        //Double_t el_energy_0{electronEnergy[0]}; //{el_energy_[0]};
        //Double_t el_energy_1{electronEnergy[1]}; //{el_energy_[1]}; // * systematic_energy_mult

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

        if(doprint)
        {
        std::cout << "electronEnergy[0]=" << electronEnergy[0]
                  << " electronEnergy[1]=" << electronEnergy[1]
                  << " highE_index=" << highE_index
                  << " lowE_index=" << lowE_index
                  << std::endl;
        }

        //if(el_energy_0 > el_energy_1)
        //{
        //    std::cout << "electron energy in wrong order" << std::endl;
        //    std::cin.get();
        //}
        
        // true energy
        // TODO: presumably this does not exist for data so need to search for
        // all instances of the trueT1, trueT2 variable and remove/replace
        // note name change
        //Double_t T0{trueT1 / bb_Q};
        //Double_t T1{trueT2 / bb_Q};

        const Double_t trueT0{trueElectronEnergy[highE_index]};
        const Double_t trueT1{trueElectronEnergy[lowE_index]};
        //const Double_t T0{trueElectronEnergy[highE_index] / bb_Q};
        //const Double_t T1{trueElectronEnergy[lowE_index] / bb_Q};
        Double_t el_energy_0{electronEnergy[highE_index]};
        Double_t el_energy_1{electronEnergy[lowE_index]};

        if(doprint)
        {
        std::cout << "trueT0=" << trueT0
                  << " trueT1=" << trueT1
                  << " el_energy_0=" << el_energy_0
                  << " el_energy_1=" << el_energy_1
                  << std::endl;
        }

        //if(trueT1 > trueT2)
        //{
        //    std::cout << "true electron energy in wrong order" << std::endl;
        //    std::cin.get();
        //}

        /*
        // if MC apply energy degradation (correction)
        // only apply to MC
        if(staticsgroup.mc_flag)
        {
            if(staticsgroup.optical_correction_enable)
            {
                // apply energy degradation (optical correction)
                Double_t visible_true_ratio_0{1.0};
                Double_t visible_true_ratio_1{1.0};

                // if energy correction systematic is enabled, choose
                // energy correction value depending on which subanalysis
                // class this is
                // 0 is default
                // 1 is high systematic
                // -1 is low systematic
                if(staticsgroup.systematic_enable_optical_correction_statistical)
                {
                    if(staticsgroup.systematic_optical_correction_statistical_direction == 0)
                    {
                        visible_true_ratio_0 = staticsgroup.g_optical_correction->Eval(1000.0 * T0);
                        visible_true_ratio_1 = staticsgroup.g_optical_correction->Eval(1000.0 * T1);
                    }
                    else if(staticsgroup.systematic_optical_correction_statistical_direction == 1)
                    {
                        visible_true_ratio_0 = staticsgroup.g_optical_correction_systematic_high->Eval(1000.0 * T0);
                        visible_true_ratio_1 = staticsgroup.g_optical_correction_systematic_high->Eval(1000.0 * T1);
                    }
                    else if(staticsgroup.systematic_optical_correction_statistical_direction == -1)
                    {
                        visible_true_ratio_0 = staticsgroup.g_optical_correction_systematic_low->Eval(1000.0 * T0);
                        visible_true_ratio_1 = staticsgroup.g_optical_correction_systematic_low->Eval(1000.0 * T1);
                    }
                    else
                    {
                        // NOOP
                        //std::cout << "invalid value of systematic_optical_correction_statistical" << std::endl;
                        //throw "invalid value of systematic_optical_correction_statistical";
                    }
                }
                else
                {
                    //std::cout << "energy correction systematic is disabled" << std::endl;
                    visible_true_ratio_0 = staticsgroup.g_optical_correction->Eval(1000.0 * T0);
                    visible_true_ratio_1 = staticsgroup.g_optical_correction->Eval(1000.0 * T1);
                }

                //std::cout << "visible_true_ratio = " << visible_true_ratio_0 << ", " << visible_true_ratio_1 << std::endl;

                // TODO this goes elsewhere
                // apply energy correction with systematics if enabled
                el_energy_0 = el_energy_0 * visible_true_ratio_0;
                el_energy_1 = el_energy_1 * visible_true_ratio_1;
            }
            else
            {
                // optical correction is diabled
                // NOOP
            }

            // TODO: other types of optical correction systematic
        
        }
        */



        /*** SYSTEMATICS **********************************************************/

        // new implementation of a linear systematic
        if(true)
        {
            //const double systematic_energy_offset = +0.001;
            //el_energy_0 = el_energy_0 + systematic_energy_offset;
            //el_energy_1 = el_energy_1 + systematic_energy_offset;
            // apply systematics to data instead
        }

    /*
        if(staticsgroup.systematic_enable_energy_multiply)
        {
            el_energy_0 *= staticsgroup.systematic_energy_multiply;
            el_energy_1 *= staticsgroup.systematic_energy_multiply;
        }

        if(staticsgroup.systematic_enable_energy_add)
        {
            el_energy_0 += staticsgroup.systematic_energy_add;
            el_energy_1 += staticsgroup.systematic_energy_add;
        }

        if(staticsgroup.systematic_enable_weight_multiply); // does nothing


        // check electron energy threshold
        if(staticsgroup.threshold_low_energy_enable)
        {
            if(el_energy_0 < staticsgroup.threshold_low_energy) continue;
            if(el_energy_1 < staticsgroup.threshold_low_energy) continue;
        }
        if(staticsgroup.threshold_low_energy_sum_enable)
        {
            if(el_energy_0 + el_energy_1 < staticsgroup.threshold_low_energy_sum) continue;
        }
    */

        // note: this was disabled in code I copied it from
        /*
        // this if statement sorts out the logical problem of having different
        // high/low sysematic energy multipliers for the purpose of using them
        // as labels to address the SubAnalysis entries in the map inside Analysis,
        // and simultaniously allowing the systematic energy mult systematic to be
        // turned off while another systematic is on
        if(systematic_energy_mult_enable == true)
        {
            el_energy_0 = el_energy_0 * systematic_energy_mult;
            el_energy_1 = el_energy_1 * systematic_energy_mult;
        }
            
        // linear energy offset systematic
        el_energy_0 = el_energy_0 + systematic_energy_offset;
        el_energy_1 = el_energy_1 + systematic_energy_offset;
        */

   /*
        // efficiency systematic
        // TODO: can remove, as weight_efficiency = systematic_efficiency
        Double_t weight_efficiency = 1.0;
        //weight_efficiency = weight_efficiency * systematic_efficiency;
        weight_efficiency = weight_efficiency * 1.0;

        // generator weight (MC weight) multiplied by weight efficiency
        Double_t aux_weight{gen_weight};
        aux_weight = aux_weight * weight_efficiency;
        
        // TODO; what happens if the energy shift / systematics move the energy
        // out of a valid range
        // answer: nothing, reweight function depends only on T1 T2
        // TODO: should T1 and T2 be shifted by systematic?
    */

        // note: already disabled
        // TODO: energy degratation systematic
        // cut both electrons > 300 keV
        /*
        if(_energy_cut_enable_)
        {
            if(el_energy_0 < 0.3) return;
            if(el_energy_1 < 0.3) return;
        }
        */


        // NOTE: more logical to set variables
        // weight_1, weight_2
        // for baseline and reweighted (now "baseline" and "test" / "universe")
        // then fill each histogram with each different weight
        // ?
        // NOTE: why would we reweight at all, why not use the decay rates from the
        // theorists directly?

        // ReWeight = baseline 0.0, ReWeight2 = baseline = 0.382
        //Double_t weight{ReWeight3(trueT1, trueT2, epsilon_31_baseline, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true")}; // TODO remove true?
        
        /*
        if(xi_31_baseline != 0.0)
        {
            std::cout << "xi_31_baseline=" << xi_31_baseline << " != 0.0" << std::endl;
            std::cin.get();
        }
        
        if((xi_31 < -0.4) || (xi_31 > 1.0))
        {
            std::cout << "xi_31=" << xi_31 << " considered to be an unusual value" << std::endl;
            std::cin.get();
        }
        */


        if(el_energy_0 < 0.3)
        {
            continue;
        }
        if(el_energy_1 < 0.3)
        {
            continue;
        }


        Double_t weight{ReWeight3(trueT0, trueT1, xi_31_baseline, xi_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true")}; // TODO remove true?

        /*
        if(weight < 1.0 - 1.0e-3)
        {
            std::cout << "weight=" << weight << std::endl;
            std::cin.get();
        }
        else if(weight > 1.0 + 1.0e-3)
        {
            std::cout << "weight=" << weight << std::endl;
            std::cin.get();
        }
        */

        if(doprint)
        {
        std::cout << "weight=" << weight << std::endl;
        }


        //Double_t weight{ReWeight3(T0, T1, epsilon_31_baseline, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true")};
        //Double_t weight{ReWeight2(T1, T2, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true")};
        //Double_t weight{ReWeight(T1, T2, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true")};

        // reweight
        //h_el_energy_reweight->Fill(el_energy_0, weight * aux_weight);
        //h_el_energy_reweight->Fill(el_energy_1, weight * aux_weight);

        //h_el_energy_sum_reweight->Fill(el_energy_0 + el_energy_1, weight * aux_weight);
        //houtput->Fill(el_energy_0 + el_energy_1, weight * aux_weight);
        //hTotalE_output->Fill(el_energy_0 + el_energy_1, 1.0 * weight);
        //hSingleEnergy_output->Fill(el_energy_0, 1.0 * weight);
        //hSingleEnergy_output->Fill(el_energy_1, 1.0 * weight);
        //hLowEnergy_output->Fill(el_energy_1, 1.0 * weight);
        //hHighEnergy_output->Fill(el_energy_0, 1.0 * weight);
        //hHighLowEnergy_output->Fill(el_energy_1, el_energy_0, 1.0 * weight);
        
        if(thePhase == 0)
        {
        /*
            hTotalE_output_P1->Fill(electronEnergy[0] + electronEnergy[1], 1.0 * weight);
            hSingleEnergy_output_P1->Fill(electronEnergy[lowE_index], 1.0 * weight);
            hSingleEnergy_output_P1->Fill(electronEnergy[highE_index], 1.0 * weight);
            hHighEnergy_output_P1->Fill(electronEnergy[highE_index], 1.0 * weight);
            hLowEnergy_output_P1->Fill(electronEnergy[lowE_index], 1.0 * weight);
            hHighLowEnergy_output_P1->Fill(electronEnergy[lowE_index], electronEnergy[highE_index], 1.0 * weight);
            hEnergySum_output_P1->Fill(electronEnergy[highE_index] + electronEnergy[lowE_index], 1.0 * weight);
            hEnergyDiff_output_P1->Fill(electronEnergy[highE_index] - electronEnergy[lowE_index], 1.0 * weight);
        */
            hTotalE_output_P1->Fill(el_energy_0 + el_energy_1, 1.0 * weight);
            hSingleEnergy_output_P1->Fill(el_energy_1, 1.0 * weight);
            hSingleEnergy_output_P1->Fill(el_energy_0, 1.0 * weight);
            hHighEnergy_output_P1->Fill(el_energy_0, 1.0 * weight);
            hLowEnergy_output_P1->Fill(el_energy_1, 1.0 * weight);
            hHighLowEnergy_output_P1->Fill(el_energy_1, el_energy_0, 1.0 * weight);
            hEnergySum_output_P1->Fill(el_energy_0 + el_energy_1, 1.0 * weight);
            hEnergyDiff_output_P1->Fill(el_energy_0 - el_energy_1, 1.0 * weight);

            hWeight_output_P1->Fill(1.0 * weight);

           // std::cout << "thePhase=0" << std::endl;
            //std::cin.get();
        }
        else if(thePhase == 1)
        {
        /*
            hTotalE_output_P2->Fill(electronEnergy[0] + electronEnergy[1], 1.0 * weight);
            hSingleEnergy_output_P2->Fill(electronEnergy[lowE_index], 1.0 * weight);
            hSingleEnergy_output_P2->Fill(electronEnergy[highE_index], 1.0 * weight);
            hHighEnergy_output_P2->Fill(electronEnergy[highE_index], 1.0 * weight);
            hLowEnergy_output_P2->Fill(electronEnergy[lowE_index], 1.0 * weight);
            hHighLowEnergy_output_P2->Fill(electronEnergy[lowE_index], electronEnergy[highE_index], 1.0 * weight);
            hEnergySum_output_P2->Fill(electronEnergy[highE_index] + electronEnergy[lowE_index], 1.0 * weight);
            hEnergyDiff_output_P2->Fill(electronEnergy[highE_index] - electronEnergy[lowE_index], 1.0 * weight);
        */
            hTotalE_output_P2->Fill(el_energy_0 + el_energy_1, 1.0 * weight);
            hSingleEnergy_output_P2->Fill(el_energy_1, 1.0 * weight);
            hSingleEnergy_output_P2->Fill(el_energy_0, 1.0 * weight);
            hHighEnergy_output_P2->Fill(el_energy_0, 1.0 * weight);
            hLowEnergy_output_P2->Fill(el_energy_1, 1.0 * weight);
            hHighLowEnergy_output_P2->Fill(el_energy_1, el_energy_0, 1.0 * weight);
            hEnergySum_output_P2->Fill(el_energy_0 + el_energy_1, 1.0 * weight);
            hEnergyDiff_output_P2->Fill(el_energy_0 - el_energy_1, 1.0 * weight);

            hWeight_output_P2->Fill(1.0 * weight);

//            std::cout << "thePhase=1" << std::endl;
        }
        else
        {
            std::cout << "Error: thePhase=" << thePhase << " invalid value" << std::endl;
        }

        // note: already removed
        /*
        if(el_energy_0 <= el_energy_1)
        {
            h_el_energy_2d_reweight->Fill(el_energy_0, el_energy_1, weight * aux_weight);
        }
        else
        {
            h_el_energy_2d_reweight->Fill(el_energy_1, el_energy_0, weight * aux_weight);
        }
        */

        // note: no systematic energy shift in this class
        //Double_t el_energy_0{el_energy_[0] * systematic_energy_mult};
        //Double_t el_energy_1{el_energy_[1] * systematic_energy_mult};
        //Double_t el_energy_0{el_energy_[0]};
        //Double_t el_energy_1{el_energy_[1]};
        // NOTE: removed, defined above, need to check


        
        // note: all below already removed

        // true energy
        // TODO: presumably this does not exist for data so need to search for
        // all instances of the trueT1, trueT2 variable and remove/replace
        //Double_t T1{trueT1 / bb_Q};
        //Double_t T2{trueT2 / bb_Q};
        // NOTE: removed, defined above, need to check

        // if MC apply energy degradation (correction)
        /*
        if(m_mode == MODE_FLAG::MODE_MC)
        {
            Double_t visible_true_ratio_1{g_energy_correction->Eval(T1)};
            Double_t visible_true_ratio_2{g_energy_correction->Eval(T2)};

            el_energy_0 = el_energy_0 * visible_true_ratio_1;
            el_energy_1 = el_energy_0 * visible_true_ratio_2;
        }
        */

        // TODO: not sure if this should be removed
        // cut both electrons > 300 keV
        /*
        if(_energy_cut_enable_)
        {
            if(el_energy_0 < 0.3) continue;
            if(el_energy_1 < 0.3) continue;
        }
        */

        // end of note, all


    }


    //print_histo_text(std::cout, "h_el_energy_reweight", h_el_energy_reweight);

    /*
    TFile *fout = new TFile("fout.root", "recreate");
    h_el_energy_reweight->Write();
    h_el_energy_sum_reweight->Write();
    fout->Close();
    std::cout << "fout" << std::endl;
    std::cin.get();

    TCanvas *c = new TCanvas("c_srw", "c_srw", 800, 600);
    h_el_energy_sum_reweight->Draw("e");
    c->SaveAs("c_srw.png");
    */

    finput->Close();

    // before scale or stack check histogram
    /*
    TString ctmp2_name;
    ctmp2_name.Form("ctmp2_%i", ctmp_counter);
    TCanvas *ctmp2 = new TCanvas(ctmp2_name);
    houtput->DrawClone();
    */



    const double TotalTime = 167629292.; // P1+P2
    const double AcceptedTime0 = 33859178.;
    const double AcceptedTime1 = 133770114.;
    const double AcceptedTime2 = TotalTime;

    // scale

    // stack

    hTotalE_output_P1->SetLineWidth(2);
    hSingleEnergy_output_P1->SetLineWidth(2);
    hHighEnergy_output_P1->SetLineWidth(2);
    hLowEnergy_output_P1->SetLineWidth(2);
    hHighLowEnergy_output_P1->SetLineWidth(2);
    hEnergySum_output_P1->SetLineWidth(2);
    hEnergyDiff_output_P1->SetLineWidth(2);
    //houtput->SetMarkerStyle(20);
    //TODO: other stuff here
    hTotalE_output_P2->SetLineWidth(2);
    hSingleEnergy_output_P2->SetLineWidth(2);
    hHighEnergy_output_P2->SetLineWidth(2);
    hLowEnergy_output_P2->SetLineWidth(2);
    hHighLowEnergy_output_P2->SetLineWidth(2);
    hEnergySum_output_P2->SetLineWidth(2);
    hEnergyDiff_output_P2->SetLineWidth(2);

    //hTotalE_output->Sumw2();
    hTotalE_output_P1->SetFillColor(Nd150Colors[0]);
    hTotalE_output_P1->SetLineColor(Nd150Colors[0]);
    hTotalE_output_P1->SetTitle(Nd150Names[0]);
    hTotalE_output_P2->SetFillColor(Nd150Colors[0]);
    hTotalE_output_P2->SetLineColor(Nd150Colors[0]);
    hTotalE_output_P2->SetTitle(Nd150Names[0]);

    //hSingleEnergy_output->Sumw2();
    hSingleEnergy_output_P1->SetFillColor(Nd150Colors[0]);
    hSingleEnergy_output_P1->SetLineColor(Nd150Colors[0]);
    hSingleEnergy_output_P1->SetTitle(Nd150Names[0]);
    hSingleEnergy_output_P2->SetFillColor(Nd150Colors[0]);
    hSingleEnergy_output_P2->SetLineColor(Nd150Colors[0]);
    hSingleEnergy_output_P2->SetTitle(Nd150Names[0]);

    //hHighEnergy_output->Sumw2();
    hHighEnergy_output_P1->SetFillColor(Nd150Colors[0]);
    hHighEnergy_output_P1->SetLineColor(Nd150Colors[0]);
    hHighEnergy_output_P1->SetTitle(Nd150Names[0]);
    hHighEnergy_output_P2->SetFillColor(Nd150Colors[0]);
    hHighEnergy_output_P2->SetLineColor(Nd150Colors[0]);
    hHighEnergy_output_P2->SetTitle(Nd150Names[0]);

    //hLowEnergy_output->Sumw2();
    hLowEnergy_output_P1->SetFillColor(Nd150Colors[0]);
    hLowEnergy_output_P1->SetLineColor(Nd150Colors[0]);
    hLowEnergy_output_P1->SetTitle(Nd150Names[0]);
    hLowEnergy_output_P2->SetFillColor(Nd150Colors[0]);
    hLowEnergy_output_P2->SetLineColor(Nd150Colors[0]);
    hLowEnergy_output_P2->SetTitle(Nd150Names[0]);

    //hHighLowEnergy_output->Sumw2();
/*
    hHighLowEnergy_output->SetFillColor(Nd150Colors[0]);
    hHighLowEnergy_output->SetLineColor(Nd150Colors[0]);
    hHighLowEnergy_output->SetTitle(Nd150Names[0]);
*/

    //hEnergySum_output->Sumw2();
    hEnergySum_output_P1->SetFillColor(Nd150Colors[0]);
    hEnergySum_output_P1->SetLineColor(Nd150Colors[0]);
    hEnergySum_output_P1->SetTitle(Nd150Names[0]);
    hEnergySum_output_P2->SetFillColor(Nd150Colors[0]);
    hEnergySum_output_P2->SetLineColor(Nd150Colors[0]);
    hEnergySum_output_P2->SetTitle(Nd150Names[0]);

    //hEnergyDiff_output->Sumw2();
    hEnergyDiff_output_P1->SetFillColor(Nd150Colors[0]);
    hEnergyDiff_output_P1->SetLineColor(Nd150Colors[0]);
    hEnergyDiff_output_P1->SetTitle(Nd150Names[0]);
    hEnergyDiff_output_P2->SetFillColor(Nd150Colors[0]);
    hEnergyDiff_output_P2->SetLineColor(Nd150Colors[0]);
    hEnergyDiff_output_P2->SetTitle(Nd150Names[0]);


    std::ifstream inFile;
    //std::string filePath = ; //in header
    //std::string filePath = "/media/ecb/Maxtor/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
    //std::string filePath = "/media/ecb/backup/users/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
    // change to using md0 raid over network
    //std::string filePath = "/mnt/md0/users/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
    
    // disabled this when putting variable further up in function
    ////std::string filePath = "/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
    
    // TODO: this will break when changing to 2eNg
    // TODO: MOVE INTO INCLUDABLE HEADER FOR BOTH fit_2e.h and newLog...h
    TString sampleFiles_150Nd = Nd150Files[0];
    //std::cout << "sampleFiles_150Nd=" << sampleFiles_150Nd << std::endl;
    // THIS WILL GO WRONG IF NOT READING 2e !
    std::string typedir = "nd150/";
    std::string fullpath = std::string(filePath.Data()) + typedir + std::string(sampleFiles_150Nd.Data()) + "/JobSummary.txt";
    //std::cout << "filepath=" << filePath << std::endl;
    //std::cout << "fullpath=" << fullpath << std::endl;
    inFile.open(fullpath);
    if(!inFile.is_open())
    {
        std::cout << "could not open: " << fullpath << std::endl;
        throw "file not open for reading number of generated MC";
    }
    std::string dummy;
    double sampleNGenMC_150Nd; // TODO: wrong type in header?
    inFile >> dummy >> sampleNGenMC_150Nd;
    inFile.close();
    //std::cout << "sampleNGemMC_150Nd=" << sampleNGenMC_150Nd << std::endl;

    hTotalE_output_P1->Scale(TotalTime / sampleNGenMC_150Nd);
    hSingleEnergy_output_P1->Scale(TotalTime / sampleNGenMC_150Nd);
    hHighEnergy_output_P1->Scale(TotalTime / sampleNGenMC_150Nd);
    hLowEnergy_output_P1->Scale(TotalTime / sampleNGenMC_150Nd);
    hHighLowEnergy_output_P1->Scale(TotalTime / sampleNGenMC_150Nd);
    hEnergySum_output_P1->Scale(TotalTime / sampleNGenMC_150Nd);
    hEnergyDiff_output_P1->Scale(TotalTime / sampleNGenMC_150Nd);
    //std::cout << "after Scale(), Integral is: " << houtput->Integral() << " scaled by " << TotalTime / sampleNGenMC_150Nd << std::endl;
    //std::cin.get();
    hTotalE_output_P2->Scale(TotalTime / sampleNGenMC_150Nd);
    hSingleEnergy_output_P2->Scale(TotalTime / sampleNGenMC_150Nd);
    hHighEnergy_output_P2->Scale(TotalTime / sampleNGenMC_150Nd);
    hLowEnergy_output_P2->Scale(TotalTime / sampleNGenMC_150Nd);
    hHighLowEnergy_output_P2->Scale(TotalTime / sampleNGenMC_150Nd);
    hEnergySum_output_P2->Scale(TotalTime / sampleNGenMC_150Nd);
    hEnergyDiff_output_P2->Scale(TotalTime / sampleNGenMC_150Nd);
    
    std::string mc_name = std::string(Nd150Files[0]);
    Double_t scale_factor_P1 = 0.0;
    Double_t scale_factor_P2 = 0.0;
    int param_number_P1 = -1;
    int param_number_P2 = -1;
    bool success = g_pg.convert_MC_name_to_param_number(mc_name, param_number_P1, param_number_P2);
    //std::cout << "mc_name=" << mc_name << " param_number=" << param_number << std::endl;
    //std::cin.get();
    if(success == true)
    //std::string search_object = MCNameToParamNameMap.at(mc_name);
    //if(paramNameToNumberMap.count(search_object) > 0)
    {
        // convert from mc sample name to param number

        //int param_number = paramNameToNumberMap.at(search_object);
        //std::cout << "parameber number " << param_number << " is in the paramNameToNumberMap" << std::endl;
        //std::cin.get();

        // TODO: change such that samples are pre-scaled by activity input value
        // get initial parameter values and error
        //Double_t param_init_value = 0.;
        //Double_t param_init_error = 0.; 
        //get_paramInitValueError(thePhase, param_number, param_init_value, param_init_error);

        //Double_t scale_factor = param_init_value;
        scale_factor_P1 = g_pg.file_params.at(param_number_P1).paramInitValue;
        scale_factor_P2 = g_pg.file_params.at(param_number_P2).paramInitValue;

        // NOTE: TODO
        // possible flaw with this method: error is no longer
        // pre-set using values from input file
        // TODO: note this in input file documentation
        // however, this may be an improvement because it
        // guarantees minuit is responsible for error estimation
        //std::cout << "scaling by: " << scale_factor << std::endl;
        //std::cin.get();
        hTotalE_output_P1->Scale(scale_factor_P1);
        hSingleEnergy_output_P1->Scale(scale_factor_P1);
        hHighEnergy_output_P1->Scale(scale_factor_P1);
        hLowEnergy_output_P1->Scale(scale_factor_P1);
        hHighLowEnergy_output_P1->Scale(scale_factor_P1);
        hEnergySum_output_P1->Scale(scale_factor_P1);
        hEnergyDiff_output_P1->Scale(scale_factor_P1);

        hTotalE_output_P2->Scale(scale_factor_P2);
        hSingleEnergy_output_P2->Scale(scale_factor_P2);
        hHighEnergy_output_P2->Scale(scale_factor_P2);
        hLowEnergy_output_P2->Scale(scale_factor_P2);
        hHighLowEnergy_output_P2->Scale(scale_factor_P2);
        hEnergySum_output_P2->Scale(scale_factor_P2);
        hEnergyDiff_output_P2->Scale(scale_factor_P2);

        // NOTE: Scale factor
        // samples are scaled by
        // TotalTime / sampleNGenMC
        // such that scale factor at this point matches that of
        // the scale factor applied in fit_2e.C to other histograms
        // and then
        // activity (Bq)
        // such that overall scale factor matches that of histograms
        // in book1DHistograms function

        /*
        std::cout << "Integrals after scaling: " << hTotalE_output->Integral()
                  << " " << hSingleEnergy_output->Integral()
                  << " " << hLowEnergy_output->Integral()
                  << " " << hHighEnergy_output->Integral()
                  << " " << hHighLowEnergy_output->Integral()
                  << std::endl;
        */
    }
    else
    {
        throw "problem";
    }
 

    /*
    TString ctmp_name;
    ctmp_name.Form("ctmp_%i", ctmp_counter);
    TCanvas *ctmp = new TCanvas(ctmp_name);
    houtput->DrawClone();
    //ctmp->Draw();
    ++ ctmp_counter;
    std::cout << "draw new hTotalE" << std::endl;
    //std:cin.get();
    */


}



#endif // NEWLOGLIKFITTER_REWEIGHT_APPLY_H
