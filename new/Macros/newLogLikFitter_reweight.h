#ifndef NEWLOGLIKFITTER_REWEIGHT_H
#define NEWLOGLIKFITTER_REWEIGHT_H





void LD_REWEIGHT_DATA(TTree *&outputG0, TTree *&outputG2, Long64_t &countG0, Long64_t &countG2)
{
    //TH1F* out = nullptr;

    TTree *t_G0 = new TTree();
    Long64_t count_G0 = t_G0->ReadFile("../../data/150Nd-data/dG150Nd/G0/dG0.dat", "electronEnergy1:electronEnergy2:electronWeight");
    std::cout << "count_G0=" << count_G0 << std::endl;

    TTree *t_G2 = new TTree();
    Long64_t count_G2 = t_G2->ReadFile("../../data/150Nd-data/dG150Nd/G2/dG2.dat", "electronEnergy1:electronEnergy2:electronWeight");
    std::cout << "count_G2=" << count_G2 << std::endl;

    //output = out;
    outputG0 = t_G0;
    outputG2 = t_G2;
    countG0 = count_G0;
    countG2 = count_G2;

    return;
}

// function dec
Double_t ReWeight3(
    const Double_t T1,
    const Double_t T2,
    const Double_t epsilon_baseline,
    const Double_t epsilon,
    const TH2D* const h_nEqNull,
    const TH2D* const h_nEqTwo,
    const Double_t psiN0,
    const Double_t psiN2,
    const std::string& debug);


// used to draw hTotalE_ after each call to reweight
// for debugging purposes
static int ctmp_counter = 0;



// xi31 reweighting function for 150Nd spectra
//void reweight_apply(TH2F *output, const TH1F *const input,
void reweight_apply(
    TH1F *&hTotalE_output,
    TH1F *&hSingleEnergy_output,
    TH1F *&hHighEnergy_output,
    TH1F *&hLowEnergy_output,
    TH1F *&hEnergySum_output,
    TH1F *&hEnergyDiff_output,
    TH2F *&hHighLowEnergy_output,
    TH1F *&hWeight_output,
    const std::string& tinput_filename,
    const Double_t xi_31,
    const Double_t xi_31_baseline,
    const TH2D* const h_nEqNull,
    const TH2D* const h_nEqTwo,
    const Double_t psiN0,
    const Double_t psiN2,
    const Double_t bb_Q)
{

    //const TString sampleName = "nd150_rot_2n2b_m4";
    const TString sampleName = Nd150Files[0];
    const TString name_append = "";


    delete gROOT->FindObject("hTotalE_" + sampleName + name_append + "_reweight");
    delete gROOT->FindObject("hSingleEnergy_" + sampleName + name_append + "_reweight");
    delete gROOT->FindObject("hHighEnergy_" + sampleName + name_append + "_reweight");
    delete gROOT->FindObject("hLowEnergy_" + sampleName + name_append + "_reweight");
    delete gROOT->FindObject("hEnergySum_" + sampleName + name_append + "_reweight");
    delete gROOT->FindObject("hEnergyDiff_" + sampleName + name_append + "_reweight");
    delete gROOT->FindObject("hHighLowEnergy_" + sampleName + name_append + "_reweight");
    delete gROOT->FindObject("hWeight_" + sampleName + name_append + "_reweight");


    // TODO: this does not work, need to re-Fill histogram using file
    //output = input->Clone(input->GetName() + "_reweight");
    hTotalE_output = new TH1F("hTotalE_" + sampleName + name_append + "_reweight",
                       //"Phase " + Phase + " " + sampleName + name_append + " total energy; #SigmaE_{e} (MeV)",
                       "Phase " + Phase + " " + sampleName + name_append + " total energy; Total Energy #SigmaE_{e} (MeV)",
                       50, 0.0, 5.0);
                       // TODO: changed from 4

    hSingleEnergy_output    = new TH1F("hSingleEnergy_" + sampleName + name_append + "_reweight",
                                "Phase " + Phase + " " + sampleName  + name_append + " Single Energy",
                                50, 0.0, 5.0);

    hHighEnergy_output     = new TH1F("hHighEnergy_" + sampleName + name_append + "_reweight",
                                "Phase " + Phase + " " + sampleName + name_append + " High Energy; Energy (MeV)",
                                50, 0.0, 5.0);

    hLowEnergy_output     = new TH1F("hLowEnergy_" + sampleName + name_append + "_reweight",
                                "Phase " + Phase + " " + sampleName + name_append + " Low Energy",
                                50, 0.0, 5.0);

    hEnergySum_output     = new TH1F("hEnergySum_" + sampleName + name_append + "_reweight",
                                "Phase " + Phase + " " + sampleName + name_append + " Low Energy",
                                50, 0.0, 5.0);

    hEnergyDiff_output     = new TH1F("hEnergyDiff_" + sampleName + name_append + "_reweight",
                                "Phase " + Phase + " " + sampleName + name_append + " Low Energy",
                                50, 0.0, 5.0);

    hHighLowEnergy_output     = new TH2F("hHighLowEnergy_" + sampleName + name_append + "_reweight",
                                "Phase " + Phase + " " + sampleName + name_append + ";Low Energy Electron Energy (MeV);High Energy Electron Energy (MeV)",
                                50, 0.0, 5.0, 50, 0.0, 5.0);

    hWeight_output       = new TH1F("hWeight_" + sampleName + name_append + "_reweight",
                                "Phase " + Phase + " " + sampleName + name_append + ";Weight",
                                //50, -2.0, 4.0);
                                50, 0.0, 0.0);

    hTotalE_output->Sumw2();
    hSingleEnergy_output->Sumw2();
    hHighEnergy_output->Sumw2();
    hLowEnergy_output->Sumw2();
    hHighLowEnergy_output->Sumw2();
    hEnergySum_output->Sumw2();
    hEnergyDiff_output->Sumw2();

    hWeight_output->Sumw2();

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
    
    TFile *finput = new TFile(tinput_filename.c_str());
    TTree *tinput = (TTree*)finput->Get("Nd150_2eNg/Nd150_2eNg");

    Int_t nElectrons;
    double electronEnergy[2];
    double trueElectronEnergy[2];

    tinput->SetBranchAddress("nElectrons"          , &nElectrons);  
    tinput->SetBranchAddress("electronEnergy"      , electronEnergy);
    tinput->SetBranchAddress("trueElectronEnergy"  , trueElectronEnergy);

    bool doprint = false;

    //const Double_t &trueT1{trueElectronEnergy[0]};
    //const Double_t &trueT2{trueElectronEnergy[1]};

    // TODO: there appears to be a very slight difference in the number of data/MC events
    // compared to what is produced by fit_2e.C
    // this is due to background MC scaling differences (in activities.txt / parameter_names.lst)
    // or is it?
    for(Long64_t ix{0}; ix < tinput->GetEntries(); ++ ix)
    {

//        if(ix == 10) break;

        tinput->GetEntry(ix);

        // analysis only valid for 2 electron events
        //if(nElectrons != 2) continue;

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
        const Double_t el_energy_0{electronEnergy[highE_index]};
        const Double_t el_energy_1{electronEnergy[lowE_index]};

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
        
        hTotalE_output->Fill(electronEnergy[0] + electronEnergy[1], 1.0 * weight);
        hSingleEnergy_output->Fill(electronEnergy[lowE_index], 1.0 * weight);
        hSingleEnergy_output->Fill(electronEnergy[highE_index], 1.0 * weight);
        hHighEnergy_output->Fill(electronEnergy[highE_index], 1.0 * weight);
        hLowEnergy_output->Fill(electronEnergy[lowE_index], 1.0 * weight);
        hHighLowEnergy_output->Fill(electronEnergy[lowE_index], electronEnergy[highE_index], 1.0 * weight);
        hEnergySum_output->Fill(electronEnergy[highE_index] + electronEnergy[lowE_index], 1.0 * weight);
        hEnergyDiff_output->Fill(electronEnergy[highE_index] - electronEnergy[lowE_index], 1.0 * weight);

        hWeight_output->Fill(1.0 * weight);

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

    hTotalE_output->SetLineWidth(2);
    hSingleEnergy_output->SetLineWidth(2);
    hHighEnergy_output->SetLineWidth(2);
    hLowEnergy_output->SetLineWidth(2);
    hHighLowEnergy_output->SetLineWidth(2);
    hEnergySum_output->SetLineWidth(2);
    hEnergyDiff_output->SetLineWidth(2);
    //houtput->SetMarkerStyle(20);
    //TODO: other stuff here

    //hTotalE_output->Sumw2();
    hTotalE_output->SetFillColor(Nd150Colors[0]);
    hTotalE_output->SetLineColor(Nd150Colors[0]);
    hTotalE_output->SetTitle(Nd150Names[0]);

    //hSingleEnergy_output->Sumw2();
    hSingleEnergy_output->SetFillColor(Nd150Colors[0]);
    hSingleEnergy_output->SetLineColor(Nd150Colors[0]);
    hSingleEnergy_output->SetTitle(Nd150Names[0]);

    //hHighEnergy_output->Sumw2();
    hHighEnergy_output->SetFillColor(Nd150Colors[0]);
    hHighEnergy_output->SetLineColor(Nd150Colors[0]);
    hHighEnergy_output->SetTitle(Nd150Names[0]);

    //hLowEnergy_output->Sumw2();
    hLowEnergy_output->SetFillColor(Nd150Colors[0]);
    hLowEnergy_output->SetLineColor(Nd150Colors[0]);
    hLowEnergy_output->SetTitle(Nd150Names[0]);

    //hHighLowEnergy_output->Sumw2();
/*
    hHighLowEnergy_output->SetFillColor(Nd150Colors[0]);
    hHighLowEnergy_output->SetLineColor(Nd150Colors[0]);
    hHighLowEnergy_output->SetTitle(Nd150Names[0]);
*/

    //hEnergySum_output->Sumw2();
    hEnergySum_output->SetFillColor(Nd150Colors[0]);
    hEnergySum_output->SetLineColor(Nd150Colors[0]);
    hEnergySum_output->SetTitle(Nd150Names[0]);

    //hEnergyDiff_output->Sumw2();
    hEnergyDiff_output->SetFillColor(Nd150Colors[0]);
    hEnergyDiff_output->SetLineColor(Nd150Colors[0]);
    hEnergyDiff_output->SetTitle(Nd150Names[0]);


    std::ifstream inFile;
    //std::string filePath = ; //in header
    //std::string filePath = "/media/ecb/Maxtor/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
    //std::string filePath = "/media/ecb/backup/users/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
    // change to using md0 raid over network
    //std::string filePath = "/mnt/md0/users/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
    std::string filePath = "/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
    // TODO: this will break when changing to 2eNg
    // TODO: MOVE INTO INCLUDABLE HEADER FOR BOTH fit_2e.h and newLog...h
    TString sampleFiles_150Nd = Nd150Files[0];
    //std::cout << "sampleFiles_150Nd=" << sampleFiles_150Nd << std::endl;
    // THIS WILL GO WRONG IF NOT READING 2e !
    std::string typedir = "nd150/";
    std::string fullpath = filePath + typedir + std::string(sampleFiles_150Nd.Data()) + "/JobSummary.txt";
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

    hTotalE_output->Scale(TotalTime / sampleNGenMC_150Nd);
    hSingleEnergy_output->Scale(TotalTime / sampleNGenMC_150Nd);
    hHighEnergy_output->Scale(TotalTime / sampleNGenMC_150Nd);
    hLowEnergy_output->Scale(TotalTime / sampleNGenMC_150Nd);
    hHighLowEnergy_output->Scale(TotalTime / sampleNGenMC_150Nd);
    hEnergySum_output->Scale(TotalTime / sampleNGenMC_150Nd);
    hEnergyDiff_output->Scale(TotalTime / sampleNGenMC_150Nd);
    //std::cout << "after Scale(), Integral is: " << houtput->Integral() << " scaled by " << TotalTime / sampleNGenMC_150Nd << std::endl;
    //std::cin.get();
    
    std::string mc_name = std::string(Nd150Files[0]);
    std::string search_object = MCNameToParamNameMap.at(mc_name);
    if(paramNameToNumberMap.count(search_object) > 0)
    {
        // convert from mc sample name to param number

        int param_number = paramNameToNumberMap.at(search_object);
        //std::cout << "parameber number " << param_number << " is in the paramNameToNumberMap" << std::endl;
        //std::cin.get();

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

        // NOTE: TODO
        // possible flaw with this method: error is no longer
        // pre-set using values from input file
        // TODO: note this in input file documentation
        // however, this may be an improvement because it
        // guarantees minuit is responsible for error estimation
        //std::cout << "scaling by: " << scale_factor << std::endl;
        //std::cin.get();
        hTotalE_output->Scale(scale_factor);
        hSingleEnergy_output->Scale(scale_factor);
        hLowEnergy_output->Scale(scale_factor);
        hHighEnergy_output->Scale(scale_factor);
        hHighLowEnergy_output->Scale(scale_factor);
        hEnergySum_output->Scale(scale_factor);
        hEnergyDiff_output->Scale(scale_factor);

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




// electron energies: T1, T2 (in units of Q value)
// theory parameter: epsilon_baseline
// theory parameter: epsilon
// data: G0, G2 (units of Q value)
//Double_t ReWeight(const Double_t T1, const Double_t T2, const Double_t epsilon,
//                  const std::vector<std::vector<double>>& data_nEqNull,
//                  const std::vector<std::vector<double>>& data_nEqTwo)
Double_t ReWeight3(const Double_t T1, const Double_t T2,
                  const Double_t epsilon_baseline, const Double_t epsilon,
                  const TH2D* const h_nEqNull,
                  const TH2D* const h_nEqTwo,
                  const Double_t psiN0, const Double_t psiN2,
                  const std::string& debug = "")
{

    // TODO: NO INTERPOLATION DONE YET

    // find bin corresponding to energies T1, T2
    Int_t bin_x{h_nEqNull->GetXaxis()->FindBin(T1)};
    Int_t bin_y{h_nEqNull->GetYaxis()->FindBin(T2)};

    //std::cout << "T1=" << T1 << " T2=" << T2 << std::endl;
    //std::cout << "the bin is " << bin_x << " " << bin_y << std::endl;
    //std::cin.get();

    // bin content corresponding to energies T1, T2
    Double_t h_nEqNull_c{h_nEqNull->GetBinContent(bin_x, bin_y)};
    Double_t h_nEqTwo_c{h_nEqTwo->GetBinContent(bin_x, bin_y)};
  
    // get the weight for this T1, T2
    // the input data is for epsilon = 0.0
    //Double_t phase_1{1.0 / psiN0};
    //Double_t weight_1{phase_1 * h_nEqNull->GetBinContent(bin_x, bin_y)};
    Double_t phase_1{1.0 / (psiN0 + epsilon_baseline * psiN2)};
    Double_t weight_1{phase_1 * (h_nEqNull_c + epsilon_baseline * h_nEqTwo_c)};

    // get the weight for this T1, T2
    // the input data is for epsilon = some arbitary value
    Double_t phase_2{1.0 / (psiN0 + epsilon * psiN2)};
    Double_t weight_2{phase_2 * (h_nEqNull_c + epsilon * h_nEqTwo_c)};

    if(debug == "true")
    {
        if(std::isnan(weight_2 / weight_1))
        {
            std::cout << "T1=" << T1 << " T2=" << T2 << std::endl;
            std::cout << "bin_x=" << bin_x << " bin_y=" << bin_y << std::endl;
            std::cout << "h_nEqNull->GetBinContent(bin_x, bin_y)=" << h_nEqNull->GetBinContent(bin_x, bin_y) << std::endl;
            std::cout << "h_nEqTwo->GetBinContent(bin_x, bin_y)=" << h_nEqTwo->GetBinContent(bin_x, bin_y) << std::endl;
            std::cout << "weight_1=" << weight_1 << " weight_2=" << weight_2 << std::endl;
            std::cout << "phase_1=" << phase_1 << " phase_2=" << phase_2 << std::endl;
            std::cout << psiN0 << " " << psiN2 << " " << epsilon << " " << epsilon_baseline << std::endl;
            std::cin.get();
        }
    }


    Double_t weight = weight_2 / weight_1;

/*
    if(std::abs(weight - 1.0) > 1.0e-5)
    {
        std::cout << std::scientific;
        std::cout << "weight=" << weight << " weight_1=" << weight_1 << " weight_2=" << weight_2 << std::endl;
        std::cout << h_nEqNull << " " << h_nEqTwo << std::endl;
        std::cout << "c: " << h_nEqNull_c << ", " << h_nEqTwo_c << std::endl;
        std::cout << "eps: " << epsilon << ", " << epsilon_baseline << std::endl;
        std::cout << "psiN0=" << psiN0 << " psiN2=" << psiN2 << std::endl;
        std::cout << "phase_1=" << phase_1 << " phase_2=" << phase_2 << std::endl;
        std::cout << "T1=" << T1 << " T2=" << T2 << " bin_x=" << bin_x << " bin_y=" << bin_y << std::endl;
        std::cin.get();
    }
*/

    // return re-weighting factor
    return weight;


}

// TODO: move elsewhere
//void reweight_apply(TH1F *output, const TH1F *const input, ...)
//{
//
//    for(Int_t bin_ix{1}; bin_ix <= input->GetNBinsX(); ++ bin_ix)
//    {
//        Double_t bin_energy = ???
//        // >>> cannot reweight in this manner
//        // >>> do not know T1, T2, need to rebuild histogram by reading from TTree
//        // or by storing as 2d histogram
//    }
//
//    ReWeight3(T1, T2, epsilon_baseline, h_nEqNull, h_nEqTwo, psiN0, psiN2, "false");
//    
//
//}

#endif // NEWLOGLIKFITTER_REWEIGHT_H
