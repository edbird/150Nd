#ifndef NEWLOGLIKFITTER_DRAW_CHANNEL_PHASE_WITH_PULL_H
#define NEWLOGLIKFITTER_DRAW_CHANNEL_PHASE_WITH_PULL_H



const int AXISMODE = 2;


// thePhase_flag:
// 0 = P1
// 1 = P2
// 2 = P1+P2
void draw_channel_phase_with_pull(
    const int channel,
    const std::string& phase_arg_str,
    const draw_input_data &drawinputdata,
    const std::vector<double> &params,
    const std::vector<double> &param_errs
    )
{

/*
    std::cout << __func__ << std::endl;
    std::cout << "params: ";
    for(int i = 0; i < params.size(); ++ i)
    {
        std::cout << params.at(i) << " ";
    }
    std::cout << std::endl;
*/

    int debuglevel = 1;

    std::vector<std::string> phase_arg_strs;
    if(phase_arg_str == "P1")
    {
        phase_arg_strs.push_back("P1");
    }
    else if(phase_arg_str == "P2")
    {
        phase_arg_strs.push_back("P2");
    }
    else if(phase_arg_str == "PBOTH")
    {
        phase_arg_strs.push_back("P1");
        phase_arg_strs.push_back("P2");
    }



    {


        // WHAT OBJECTS DO I WANT TO DRAW
        //
        // TOTAL MC as a BLACK LINE (THIS IS A STACK)
        // TOTAL DATA/FAKEDATA as BLACK ERROR BAR
        // STACK UP TO 2nbb
        // STACK UP TO 208Tl
        // STACK UP TO 214Bi
        // STACK UP TO 207Bi
        // STACK UP TO INTERNAL
        // STACK UP TO NEIGHBOUR
        // STACK UP TO RADON
        // STACK UP TO EXTERNAL
        //
        // RATIO OF DATA/FAKE DATA DIV BY TOTAL MC
        //
        // PULL HISTOGRAM: P = (DATA/FAKEDATA - TOTAL MC) / sqrt(TOTAL MC)


        ///////////////////////////////////////////////////////////////////////////
        // Draw phase as specified by phase_arg_str
        // which can be one of:
        // P1, P2, PBOTH
        ///////////////////////////////////////////////////////////////////////////

        //THStack *stacks1D_Px;
        /*
        TH1D *h_2nubb_Px = nullptr;
        TH1D *h_tl208_int_Px = nullptr;
        TH1D *h_bi214_int_Px = nullptr;
        TH1D *h_bi207_int_Px = nullptr;
        TH1D *h_internal_Px = nullptr;
        TH1D *h_external_Px = nullptr;
        TH1D *h_radon_Px = nullptr; 
        TH1D *h_neighbour_Px = nullptr;
        TH1D *h_other_Px = nullptr;
        */

        // custom stack
        TH1D *h_stack_external = nullptr;
        TH1D *h_stack_radon = nullptr;
        TH1D *h_stack_neighbour = nullptr;
        TH1D *h_stack_internal = nullptr;
        TH1D *h_stack_bi207_int = nullptr;
        TH1D *h_stack_bi214_int = nullptr;
        TH1D *h_stack_tl208_int = nullptr;
        TH1D *h_stack_2nubb = nullptr;
        TH1D *h_stack_total_MC = nullptr;

        TH1D *h_data_fakedata = nullptr;

        TH1D *h_ratio = nullptr;
        TH1D *h_pull = nullptr;

        TCanvas *canvas_Px = nullptr;
        TPad *p0_Px = nullptr;
        TPad *p1_Px = nullptr;
        TPad *p2_Px = nullptr;
        /*
        TH1D *hRatio_Px = nullptr;
        TH1D *hSig_Px = nullptr;
//        TH1D *hAllMC1D_Px = nullptr;
        TH1D *data1D_Px = nullptr;
        TH1D *fakeData1D_Px = nullptr;
        */



        ///////////////////////////////////////////////////////////////////////
        // ALLOC OBJECTS
        ///////////////////////////////////////////////////////////////////////
    
        TString the_histname = channel_histname_1D[channel];

        TString name_stack_external = the_histname + "_stack_external_" + phase_arg_str;
        h_stack_external = new TH1D(name_stack_external, name_stack_external, 50, 0.0, 5.0);
    
        TString name_stack_radon = the_histname + "_stack_radon_" + phase_arg_str;
        h_stack_radon = new TH1D(name_stack_radon, name_stack_radon, 50, 0.0, 5.0);
    
        TString name_stack_neighbour = the_histname + "_stack_neighbour_" + phase_arg_str;
        h_stack_neighbour = new TH1D(name_stack_neighbour, name_stack_neighbour, 50, 0.0, 5.0);
    
        TString name_stack_internal = the_histname + "_stack_internal_" + phase_arg_str;
        h_stack_internal = new TH1D(name_stack_internal, name_stack_internal, 50, 0.0, 5.0);
    
        TString name_stack_bi207_int = the_histname + "_stack_bi207_int_" + phase_arg_str;
        h_stack_bi207_int = new TH1D(name_stack_bi207_int, name_stack_bi207_int, 50, 0.0, 5.0);
    
        TString name_stack_bi214_int = the_histname + "_stack_bi214_int_" + phase_arg_str;
        h_stack_bi214_int = new TH1D(name_stack_bi214_int, name_stack_bi214_int, 50, 0.0, 5.0);
    
        TString name_stack_tl208_int = the_histname + "_stack_tl208_int_" + phase_arg_str;
        h_stack_tl208_int = new TH1D(name_stack_tl208_int, name_stack_tl208_int, 50, 0.0, 5.0);
    
        TString name_stack_2nubb = the_histname + "_stack_2nubb_" + phase_arg_str;
        h_stack_2nubb = new TH1D(name_stack_2nubb, name_stack_2nubb, 50, 0.0, 5.0);
    
        TString name_stack_total_MC = the_histname + "_stack_total_MC_" + phase_arg_str;
        h_stack_total_MC = new TH1D(name_stack_total_MC, name_stack_total_MC, 50, 0.0, 5.0);
    
        TString name_data_fakedata = the_histname + "_data_fakedata_" + phase_arg_str;
        h_data_fakedata = new TH1D(name_data_fakedata, name_data_fakedata, 50, 0.0, 5.0);
    
        TString name_ratio = the_histname + "_ratio_" + phase_arg_str;
        h_ratio = new TH1D(name_ratio, name_ratio, 50, 0.0, 5.0);
    
        TString name_pull = the_histname + "_pull_" + phase_arg_str;
        h_pull = new TH1D(name_pull, name_pull, 50, 0.0, 5.0);


        ///////////////////////////////////////////////////////////////////////
        // SET COMMON PROPERTIES
        ///////////////////////////////////////////////////////////////////////

        h_stack_external->SetLineWidth(0);
        h_stack_external->SetFillStyle(1001);
        h_stack_external->SetFillColor(ExternalBkgColor);
        h_stack_external->SetTitle(0);
        h_stack_external->SetStats(0);

        h_stack_radon->SetLineWidth(0);
        h_stack_radon->SetFillStyle(1001);
        h_stack_radon->SetFillColor(RadonBkgColor);
        h_stack_radon->SetTitle(0);
        h_stack_radon->SetStats(0);

        h_stack_neighbour->SetLineWidth(0);
        h_stack_neighbour->SetFillStyle(1001);
        h_stack_neighbour->SetFillColor(NeighbourColor);
        h_stack_neighbour->SetTitle(0);
        h_stack_neighbour->SetStats(0);

        h_stack_internal->SetLineWidth(0);
        h_stack_internal->SetFillStyle(1001);
        h_stack_internal->SetFillColor(InternalBkgColor);
        h_stack_internal->SetTitle(0);
        h_stack_internal->SetStats(0);

        h_stack_bi207_int->SetLineWidth(0);
        h_stack_bi207_int->SetFillStyle(1001);
        h_stack_bi207_int->SetFillColor(bi207InternalBkgColor);
        h_stack_bi207_int->SetTitle(0);
        h_stack_bi207_int->SetStats(0);

        h_stack_bi214_int->SetLineWidth(0);
        h_stack_bi214_int->SetFillStyle(1001);
        h_stack_bi214_int->SetFillColor(bi214InternalBkgColor);
        h_stack_bi214_int->SetTitle(0);
        h_stack_bi214_int->SetStats(0);

        h_stack_tl208_int->SetLineWidth(0);
        h_stack_tl208_int->SetFillStyle(1001);
        h_stack_tl208_int->SetFillColor(tl208InternalBkgColor);
        h_stack_tl208_int->SetTitle(0);
        h_stack_tl208_int->SetStats(0);

        h_stack_2nubb->SetLineWidth(0);
        h_stack_2nubb->SetFillStyle(1001);
        h_stack_2nubb->SetFillColor(Nd150Color);
        h_stack_2nubb->SetTitle(0);
        h_stack_2nubb->SetStats(0);

        h_stack_total_MC->SetLineWidth(2);
        h_stack_total_MC->SetLineColor(kBlack);
        h_stack_total_MC->SetTitle(0);
        h_stack_total_MC->SetStats(0);

        h_data_fakedata->SetLineWidth(2);
        h_data_fakedata->SetLineColor(kBlack);
        h_data_fakedata->SetTitle(0);
        h_data_fakedata->SetStats(0);

        h_ratio->SetLineWidth(2);
        h_ratio->SetLineColor(kBlack);
        h_ratio->SetTitle(0);
        h_ratio->SetStats(0);

        h_pull->SetLineWidth(2);
        h_pull->SetLineColor(kBlack);
        h_pull->SetTitle(0);
        h_pull->SetStats(0);





        //std::cout << "debug: number of data samples: " << allDataSamples1D->GetEntries() << std::endl;
        //std::cout << "debug: number of MC samples: " << allMCSamples1D[0]->GetEntries() << std::endl;


        //std::cout << "g_mode_fake_data=" << g_mode_fake_data << std::endl;
        //std::cout << allDataSamples1D->GetEntries() << std::endl;
        //std::cin.get();


        // additional array index
        // TODO: fix this to get Px data
        std::string histname = std::string(channel_histname_1D[channel]);
        for(std::vector<std::string>::iterator it{phase_arg_strs.begin()};
            it != phase_arg_strs.end(); ++ it)
        {
            std::string search_object_Px;
            TH1D *tmpHist1D = nullptr;
            if(g_mode_fake_data == false)
            {
                search_object_Px = histname + std::string(DataFile) + "_" + (*it);
                tmpHist1D = (TH1D*)allDataSamples1D->FindObject(search_object_Px.c_str());
            }
            else if(g_mode_fake_data == true)
            {
                search_object_Px = histname + std::string("fakedata") + "_" + (*it);
                tmpHist1D = (TH1D*)allFakeDataSamples1D->FindObject(search_object_Px.c_str());
            }
            //std::cout << "search_object_Px=" << search_object_Px << std::endl;
            //std::string search_object_P2 = histname + std::string(DataFile) + "_P2";
            if(tmpHist1D != nullptr)
            {
                h_data_fakedata->Add((TH1D*)tmpHist1D, 1.0);
            }
            else
            {
                std::cout << "ERROR: Could not find histogram: " << search_object_Px << std::endl;
            }

            ///////////////////////////////////////////////////////////////////////
            // manually create a THStack since the programmers of ROOT were
            // incapable of programming such a thing properly
            ///////////////////////////////////////////////////////////////////////
        }

        // init errors for data/fakedata stack
        for(Int_t i = 0; i < h_data_fakedata->GetNbinsX(); ++ i)
        {
            Double_t content = h_data_fakedata->GetBinContent(i);
            h_data_fakedata->SetBinError(i, std::sqrt(content));
        }

        TString channel_str;
        channel_str.Form("ch%i", channel);

        canvas_Px = new TCanvas("c_" + channel_str + "_" + phase_arg_str, "c_" + channel_str + "_" + phase_arg_str);
        canvas_Px->SetFillColor(kWhite);

        
        // these are the minor stacks?
//        stacks1D_Px = new THStack("stacks1D_" + channel_str + "_" + phase_arg_str, "stacks1D_" + channel_str + "_" + phase_arg_str);
//        stacks1D_Px = nullptr;


        // loop over all the parameters
        std::map<int, file_parameter>::iterator it{g_pg.file_params.begin()};
        for(; it != g_pg.file_params.end(); ++ it)
        {
            int paramNumberInt = -1;

            int paramNumber = it->second.paramNumber;
            bool paramEnabled = it->second.paramEnabled;
            bool paramEnabledP1 = false;
            if(phase_arg_str == "P1" || phase_arg_str == "PBOTH")
            {
                if(it->second.paramEnabledP1 == true)
                {
                    paramEnabledP1 = true;
                }
            }
            bool paramEnabledP2 = false;
            if(phase_arg_str == "P2" || phase_arg_str == "PBOTH")
            {
                if(it->second.paramEnabledP2 == true)
                {
                    paramEnabledP2 = true;
                }
            }

            if(debuglevel >= 2)
            {
                std::cout << "paramNumber=" << paramNumber << " P1: " << paramEnabledP1 << " P2: " << paramEnabledP2 << std::endl;
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
                continue;
            }

            std::vector<std::string>::iterator mc_name_it{it->second.MCNameList.begin()};
            for(; mc_name_it != it->second.MCNameList.end(); ++ mc_name_it)
            {
                std::string mc_name = *mc_name_it;
                std::string histname = std::string(channel_histname_1D[channel]);
                std::string search_object_P1 = histname + mc_name + "_P1_fit";
                std::string search_object_P2 = histname + mc_name + "_P2_fit";
                TH1D *tmpHist1D_P1 = nullptr;
                TH1D *tmpHist1D_P2 = nullptr;

                paramNumberInt = g_pg.ExtToIntParamNumberMap.at(paramNumber);
                if(debuglevel >= 2)
                {
                    std::cout << "paramNumber=" << paramNumber << " -> " << paramNumberInt << std::endl;
                }


                if(paramEnabledP1 == true)
                {
                    //TString name = TString(allMCSamples1D[channel]->FindObject(search_object_P1.c_str())->GetName()) + "_tmp";
                    //tmpHist1D_P1 = (TH1D*)allMCSamples1D[channel]->FindObject(search_object_P1.c_str())->Clone(name); // need to clone to scale
                    tmpHist1D_P1 = (TH1D*)allMCSamples1D[channel]->FindObject(search_object_P1.c_str());
                    
                    if(tmpHist1D_P1 == nullptr)
                    {
                        std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
                        throw "problem";
                    }

                    // TODO: there was a check here for Integral() > 0
                    // however I removed it

                    Double_t scale_factor_P1 = params.at(paramNumberInt);
                    //tmpHist1D_P1->Scale(scale_factor_P1);

                    new_stacker_helper(
                        tmpHist1D_P1,
                        h_stack_external,
                        h_stack_radon,
                        h_stack_neighbour,
                        h_stack_internal,
                        h_stack_bi207_int,
                        h_stack_bi214_int,
                        h_stack_tl208_int,
                        h_stack_2nubb,
                        h_stack_total_MC,
                        scale_factor_P1
                    );
                }

                
                if(paramEnabledP2 == true)
                {
                    //TString name = TString(allMCSamples1D[channel]->FindObject(search_object_P2.c_str())->GetName()) + "_tmp";
                    //tmpHist1D_P2 = (TH1D*)allMCSamples1D[channel]->FindObject(search_object_P2.c_str())->Clone(name); // need clone to scale
                    tmpHist1D_P2 = (TH1D*)allMCSamples1D[channel]->FindObject(search_object_P2.c_str());
                
                    if(tmpHist1D_P2 == nullptr)
                    {
                        std::cout << "ERROR: Could not find object " << search_object_P2 << std::endl;
                        throw "problem";
                    }

                    // TODO: there was a check here for Integral() > 0
                    // however I removed it

                    Double_t scale_factor_P2 = params.at(paramNumberInt);
                    
                    //if(tmpHist1D_P2->Integral() > 0)
                    //{
                    new_stacker_helper(
                        tmpHist1D_P2,
                        h_stack_external,
                        h_stack_radon,
                        h_stack_neighbour,
                        h_stack_internal,
                        h_stack_bi207_int,
                        h_stack_bi214_int,
                        h_stack_tl208_int,
                        h_stack_2nubb,
                        h_stack_total_MC,
                        scale_factor_P2
                    );
                    //}
                    //else
                    //{
                        //std::cout << "not adding to stack, Integral() <= 0: " << tmpHist_draw1D->GetName() << std::endl;
                    //}
                }

            }
        }



        ///////////////////////////////////////////////////////////////////////
        // END OF STACKING CODE
        // START OF CANVAS DRAWING CODE
        ///////////////////////////////////////////////////////////////////////


        double PAD_U_Y_MIN_P1_P2 = 0.0;
        double PAD_U_Y_MAX_P1_P2 = 500.0;
        double PAD_L_Y_MIN_P1_P2 = 0.5;
        double PAD_L_Y_MAX_P1_P2 = 2.0;
        double PAD_LL_Y_MIN_P1_P2 = -3.0;
        double PAD_LL_Y_MAX_P1_P2 = 3.0;

        double PAD_U_Y_MIN_P1 = 0.0;
        double PAD_U_Y_MAX_P1 = 500.0;
        double PAD_L_Y_MIN_P1 = 0.5;
        double PAD_L_Y_MAX_P1 = 2.0;
        double PAD_LL_Y_MIN_P1 = -3.0;
        double PAD_LL_Y_MAX_P1 = 3.0;

        double PAD_U_Y_MIN_P2 = 0.0;
        double PAD_U_Y_MAX_P2 = 500.0;
        double PAD_L_Y_MIN_P2 = 0.5;
        double PAD_L_Y_MAX_P2 = 2.0;
        double PAD_LL_Y_MIN_P2 = -3.0;
        double PAD_LL_Y_MAX_P2 = 3.0;

        double PAD_U_Y_MIN_Px = 0.0;
        double PAD_U_Y_MAX_Px = 0.0;
        double PAD_L_Y_MIN_Px = 0.0;
        double PAD_L_Y_MAX_Px = 0.0;
        double PAD_LL_Y_MIN_Px = 0.0;
        double PAD_LL_Y_MAX_Px = 0.0;

        if(channel == 0)
        {
            PAD_U_Y_MAX_P1_P2 = 300.0;

            PAD_U_Y_MAX_P1 = 80;
            PAD_U_Y_MAX_P2 = 250;
        }
        else if(channel == 1)
        {
            PAD_U_Y_MAX_P1_P2 = 1000.0;

            PAD_U_Y_MAX_P1 = 200;
            PAD_U_Y_MAX_P2 = 800;
        }
        else if(channel == 2)
        {
            PAD_U_Y_MAX_P1_P2 = 400.0;

            PAD_U_Y_MAX_P1 = 80;
            PAD_U_Y_MAX_P2 = 350;
        }
        else if(channel == 3)
        {
            PAD_U_Y_MAX_P1_P2 = 900.0;

            PAD_U_Y_MAX_P1 = 20;
            PAD_U_Y_MAX_P2 = 800;
        }
        else if(channel == 4)
        {
            PAD_U_Y_MAX_P1_P2 = 300.0;

            PAD_U_Y_MAX_P1 = 80;
            PAD_U_Y_MAX_P2 = 250;
        }
        else if(channel == 5)
        {
            PAD_U_Y_MAX_P1_P2 = 550.0;

            PAD_U_Y_MAX_P1 = 150;
            PAD_U_Y_MAX_P2 = 450;
        }
        else
        {
            PAD_U_Y_MAX_P1_P2 = 350.0;

            PAD_U_Y_MAX_P1 = 100;
            PAD_U_Y_MAX_P2 = 250;
        }

        if(phase_arg_str == "P1")
        {
            PAD_U_Y_MIN_Px = PAD_U_Y_MIN_P1;
            PAD_U_Y_MAX_Px = PAD_U_Y_MAX_P1;
            PAD_L_Y_MIN_Px = PAD_L_Y_MIN_P1;
            PAD_L_Y_MAX_Px = PAD_L_Y_MAX_P1;
            PAD_LL_Y_MIN_Px = PAD_LL_Y_MIN_P1;
            PAD_LL_Y_MAX_Px = PAD_LL_Y_MAX_P1;
        }
        else if(phase_arg_str == "P2")
        {
            PAD_U_Y_MIN_Px = PAD_U_Y_MIN_P2;
            PAD_U_Y_MAX_Px = PAD_U_Y_MAX_P2;
            PAD_L_Y_MIN_Px = PAD_L_Y_MIN_P2;
            PAD_L_Y_MAX_Px = PAD_L_Y_MAX_P2;
            PAD_LL_Y_MIN_Px = PAD_LL_Y_MIN_P2;
            PAD_LL_Y_MAX_Px = PAD_LL_Y_MAX_P2;
        }
        else if(phase_arg_str == "PBOTH")
        {
            PAD_U_Y_MIN_Px = PAD_U_Y_MIN_P1_P2;
            PAD_U_Y_MAX_Px = PAD_U_Y_MAX_P1_P2;
            PAD_L_Y_MIN_Px = PAD_L_Y_MIN_P1_P2;
            PAD_L_Y_MAX_Px = PAD_L_Y_MAX_P1_P2;
            PAD_LL_Y_MIN_Px = PAD_LL_Y_MIN_P1_P2;
            PAD_LL_Y_MAX_Px = PAD_LL_Y_MAX_P1_P2;
        }



        
        /*
        h_2nubb_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px);
        h_tl208_int_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px);
        h_bi214_int_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px);
        h_bi207_int_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px);
        h_internal_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px);
        h_neighbour_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px);
        h_radon_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px);
        h_external_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px);
        */

        /*
        THStack *stacks1D_major_Px;
        stacks1D_major_Px = new THStack("stacks1D_major_" + channel_str + "_Px", "stacks1D_major_" + channel_str + "_Px");
        stacker_helper_2(stacks1D_major_Px,
                         h_2nubb_Px,
                         h_tl208_int_Px,
                         h_bi214_int_Px,
                         h_bi207_int_Px,
                         h_internal_Px,
                         h_neighbour_Px,
                         h_radon_Px,
                         h_external_Px,
                         h_other_Px);
        */

// TODO: check this works?
////////        stacks1D_major_Px->SetMaximum(PAD_U_Y_MAX_Px);
////////        stacks1D_major_Px->SetMinimum(PAD_U_Y_MIN_Px);


        //stacks1D_major[i]->Draw("hist");
        //stacks1D_major[i]->SetMaximum(PAD_U_Y_MAX);
        //stacks1D_major[i]->SetMinimum(PAD_U_Y_MIN);
    //    stacks1D_major->GetYaxis()->SetRangeUser(PAD_U_Y_MIN, PAD_U_Y_MAX);
    /*
        hAllMC1D_Px->SetMaximum(PAD_U_Y_MAX_Px);
        hAllMC1D_Px->SetMinimum(PAD_U_Y_MIN_Px);
        hAllMC1D_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px);
    */
////////        stacks1D_major_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px);
    //    data1D->SetMaximum(PAD_U_Y_MAX);
    //    data1D->SetMinimum(PAD_U_Y_MIN);
    /*
        if(g_mode_fake_data == false)
        {
            data1D_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px);
        }
        else if(g_mode_fake_data == true)
        {
            fakeData1D_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px);
        }
    */
    /*
        if(g_mode_fake_data == false)
        {
            TString name = TString(data1D_Px->GetName()) + "_ratio";
            hRatio_Px = (TH1D*)data1D_Px->Clone(name);
        }
        else if(g_mode_fake_data == true)
        {
            TString name = TString(fakeData1D_Px->GetName()) + "_ratio";
            hRatio_Px = (TH1D*)fakeData1D_Px->Clone(name);
        }
        //hRatio->Sumw2();
        hRatio_Px->Divide(hAllMC1D_Px);
        hRatio_Px->SetTitle("");
    */








        h_ratio->Add((TH1D*)h_data_fakedata, 1.0);
        //hRatio->Sumw2();
        h_ratio->Divide(h_stack_total_MC);
        for(Int_t i = 1; i <= h_ratio->GetNbinsX(); ++ i)
        {
            Double_t content_M = h_stack_total_MC->GetBinContent(i);
            Double_t content_D = h_data_fakedata->GetBinContent(i);
            //Double_t error_M = h_stack_total_MC->GetBinError(i);
            //Double_t error_M = std::sqrt(h_stack_total_MC->GetBinContent(i));
            //Double_t error_D = h_data_fakedata->GetBinError(i);
            //Double_t content = std::pow(1.0 / content_M, 2.0) * std::pow(error_D, 2.0)
            // error on D is effective zero
            if(content_M <= 0.0)
            {
                h_ratio->SetBinContent(i, 0.0);
                h_ratio->SetBinError(i, 0.0);
            }
            else
            {
                Double_t error_M = std::sqrt(h_stack_total_MC->GetBinContent(i));
                Double_t error_sq = std::pow(content_D / (content_M * content_M), 2.0) * std::pow(error_M, 2.0);
                h_ratio->SetBinError(i, std::sqrt(error_sq));
            }
        }
//        hRatio_Px->SetTitle("");
        h_pull->Add((TH1D*)h_data_fakedata, 1.0);
        h_pull->Add((TH1D*)h_stack_total_MC, -1.0);
        for(Int_t i = 1; i <= h_pull->GetNbinsX(); ++ i)
        {
            Double_t content = h_pull->GetBinContent(i);
            Double_t content_M = h_stack_total_MC->GetBinContent(i);
            Double_t error = std::sqrt(content_M);
            if(error <= 0.0)
            {
                h_pull->SetBinContent(i, 0.0);
            }
            else
            {
                h_pull->SetBinContent(i, content / error);
            }
        }

    /*
        for(Int_t i = 1; i <= hSig_Px->GetNbinsX(); ++ i)
        { 
            Double_t content_data{0.0};
            if(g_mode_fake_data == false)
            {
                content_data = data1D_Px->GetBinContent(i);
            }
            else if(g_mode_fake_data == true)
            {
                content_data = fakeData1D_Px->GetBinContent(i);
            }
            Double_t content_MC{hAllMC1D_Px->GetBinContent(i)};
            Double_t error_data{0.0};
            if(g_mode_fake_data == false)
            {
                error_data = data1D_Px->GetBinError(i);
            }
            else if(g_mode_fake_data == true)
            {
                error_data = fakeData1D_Px->GetBinError(i);
            }
            //Double_t error_MC{hAllMC1D_Px->GetBinError(i)};
            Double_t error_MC{std::sqrt(hAllMC1D_Px->GetBinContent(i))};
            //Double_t error{std::sqrt(error_fixed * error_fixed + error_var * error_var)}; // TODO: this is wrong
            Double_t error{std::sqrt(error_MC * error_MC)}; // TODO: this is wrong
            Double_t content{content_data - content_MC};
            if((error <= 0.0) && (content != 0.0))
            {
                std::cout << "problem" << std::endl;
            }
            else if((error > 0.0) && (content != 0.0))
            {
//                hSig_Px->SetBinContent(i, content / error);
            }
            else
            {
                //std::cout << "content=0.0 -> ignore" << std::endl;
//                hSig_Px->SetBinContent(i, 0.0);
            }
//            hSig_Px->SetBinError(i, 0.0);
        }

        hSig_Px->SetFillColor(kBlack);
        hSig_Px->SetFillStyle(3002);
    */


        ///////////////////////////////////////////////////////////////////////
        // ALL APPEARANCE RELATED FUNCTION CALLS
        // EXCLUDING THOSE WHICH ARE UNIVERSAL AND APPLIED WHEN ALLOCATING
        ///////////////////////////////////////////////////////////////////////


        if(channel == 0)
        {
            h_pull->GetXaxis()->SetTitle("Total Electron Energy [MeV]");
        }
        else if(channel == 1)
        {
            h_pull->GetXaxis()->SetTitle("Single Electron Energy [MeV]");
        }
        else if(channel == 2)
        {
            h_pull->GetXaxis()->SetTitle("High Energy Electron [MeV]");
        }
        else if(channel == 3)
        {
            h_pull->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
        }
        else if(channel == 4)
        {
            h_pull->GetXaxis()->SetTitle("Electron Energy Sum [MeV]");
        }
        else if(channel == 5)
        {
            h_pull->GetXaxis()->SetTitle("Electron Energy Diff [MeV]");
        }

        h_pull->SetFillColor(kBlack);
        h_pull->SetFillStyle(3002);

    /*
        //hRatio[i]->SetMaximum(PAD_L_Y_MAX);
        //hRatio[i]->SetMinimum(PAD_L_Y_MIN);
        hRatio_Px->GetYaxis()->SetRangeUser(PAD_L_Y_MIN_Px, PAD_L_Y_MAX_Px);

//        hRatio_Px->GetXaxis()->SetTickLength(0.03 / (0.2));
//        hRatio_Px->GetYaxis()->SetTickLength(0.01 / (0.2));
//        hRatio_Px->GetXaxis()->SetTickLength(0.02 * 213. / 93.);
        //hRatio_Px->GetXaxis()->SetTickSize(0.1);

        hRatio_Px->GetXaxis()->SetLabelFont(43);
        hRatio_Px->GetXaxis()->SetLabelSize(15);
        hRatio_Px->GetXaxis()->SetTitleFont(43);
        hRatio_Px->GetXaxis()->SetTitleSize(20);
        hRatio_Px->GetXaxis()->SetTitleOffset(3.25);

        hRatio_Px->GetYaxis()->SetLabelFont(43);
        hRatio_Px->GetYaxis()->SetTitleFont(43);
        hRatio_Px->GetYaxis()->SetTitleSize(20);
        hRatio_Px->GetYaxis()->SetTitleOffset(1.0);
        hRatio_Px->GetYaxis()->SetTitle("data / MC");

    */


        // TODO: fix the error bar with for ratio so that it is the same
        // as the line width for the main histogram
        // and remove the horizontal lines from both



        //hSig[i]->SetMaximum(PAD_L_Y_MAX);
        //hSig[i]->SetMinimum(PAD_L_Y_MIN);

/*
        PAD_LL_Y_MIN_Px = 0.0;
        PAD_LL_Y_MAX_Px = 0.0;
        for(Int_t i = 0; i < hSig_Px->GetNbinsX(); ++ i)
        {
            Double_t content = hSig_Px->GetBinContent(i);
            if(content > PAD_LL_Y_MAX_Px)
            {
                PAD_LL_Y_MAX_Px = content;
            }
            if(content < PAD_LL_Y_MIN_Px)
            {
                PAD_LL_Y_MIN_Px = content;
            }
        }
        if(PAD_LL_Y_MAX_Px <= PAD_LL_Y_MIN_Px)
        {
            PAD_LL_Y_MAX_Px = 0.5 * (PAD_LL_Y_MAX_Px + PAD_LL_Y_MIN_Px);
            PAD_LL_Y_MIN_Px = PAD_LL_Y_MAX_Px;
        std::cout << "CATCH" << std::endl;
            PAD_LL_Y_MAX_Px += 0.5;
            PAD_LL_Y_MIN_Px -= 0.5;
        std::cout << "MAX=" << PAD_LL_Y_MAX_Px << std::endl;
        std::cout << "MIN=" << PAD_LL_Y_MIN_Px << std::endl;
        std::cin.get();
        std::cout << std::endl;
        }
        if(PAD_LL_Y_MAX_Px <= 0.0 || PAD_LL_Y_MIN_Px >= 0.0)
        {
        std::cout << "STILL BROKEN" << std::endl;
            PAD_LL_Y_MAX_Px = 1.0;
            PAD_LL_Y_MAX_Px = -1.0;
        }
*/

    /*
    
        hSig_Px->GetYaxis()->SetRangeUser(PAD_LL_Y_MIN_Px, PAD_LL_Y_MAX_Px);
        hSig_Px->SetMaximum(PAD_LL_Y_MAX_Px);
        hSig_Px->SetMinimum(PAD_LL_Y_MIN_Px);

        hSig_Px->GetXaxis()->SetLabelFont(43);
        hSig_Px->GetXaxis()->SetLabelSize(15);
        hSig_Px->GetXaxis()->SetTitleFont(43);
        hSig_Px->GetXaxis()->SetTitleSize(20);
        hSig_Px->GetXaxis()->SetTitleOffset(3.25);

        hSig_Px->GetYaxis()->SetLabelFont(43);
        hSig_Px->GetYaxis()->SetTitleFont(43);
        hSig_Px->GetYaxis()->SetTitleSize(20);
        hSig_Px->GetYaxis()->SetTitleOffset(1.0);
        hSig_Px->GetYaxis()->SetTitle("#chi^{2} pull (#sigma)");

        if(AXISMODE == 2)
        {
            hSig_Px->GetYaxis()->SetLabelSize(0.0);
            hSig_Px->GetYaxis()->SetTickLength(0.0);
        }
        else if(AXISMODE == 1)
        {
            hSig_Px->GetYaxis()->SetLabelSize(15);
            hSig_Px->GetYaxis()->SetTickLength(1.0 / 93.0);
        }

    */

        ///////////////////////////////////////////////////////////////////////
        // canvas construction
        ///////////////////////////////////////////////////////////////////////


        canvas_Px->cd();
        canvas_Px->SetBottomMargin(0.0);

        // in code copying from, canvas alloc. here

        canvas_Px->cd();
        //canvas_Px->SetBottomMargin(0.0);
        p0_Px = new TPad("pad0_" + channel_str + "_" + phase_arg_str, "pad0" + channel_str + "_" + phase_arg_str, 0.0, 0.5, 1.0, 1.0);
        p0_Px->SetBottomMargin(0.0);
        //p0->SetGridx(1);
        //p0->SetGridy(1);
        p0_Px->SetGrid(0, 0);
        p0_Px->SetTicks(2, 2);
        p0_Px->SetRightMargin(0.05);
        p0_Px->Draw();

        canvas_Px->cd();
        p1_Px = new TPad("pad1_" + channel_str + "_" + phase_arg_str, "pad1_" + channel_str + "_" + phase_arg_str, 0.0, 0.3, 1.0, 0.5);
        p1_Px->SetTopMargin(0.0);
        p1_Px->SetBottomMargin(0.0);
        //p1->SetGridx(1);
        //p1->SetGridy(1);
        p1_Px->SetGrid(0, 0);
        p1_Px->SetTicks(2, 2);
        p1_Px->SetRightMargin(0.05); //// CHECK THIS IS WORKING
        p1_Px->Draw();

        canvas_Px->cd();
        p2_Px = new TPad("pad2_" + channel_str + "_" + phase_arg_str, "pad2_" + channel_str + "_" + phase_arg_str, 0.0, 0.0, 1.0, 0.3);
        p2_Px->SetTopMargin(0.0);
        p2_Px->SetBottomMargin(1.0/3.0);
        //p2->SetGridx(1);
        //p2->SetGridy(1);
        p2_Px->SetGrid(0, 0);
        p2_Px->SetTicks(2, 2);
        p2_Px->SetRightMargin(0.05);
        p2_Px->Draw();

        canvas_Px->cd();
        //TPad *p3_Px = new TPad("pad3_" + channel_str + "_" + phase_arg_str, "pad3_" + channel_str + "_" + phase_arg_str, 0.615, 0.4, 0.855, 0.85);
        TPad *p3_Px = new TPad("pad3_" + channel_str + "_" + phase_arg_str, "pad3_" + channel_str + "_" + phase_arg_str, 0.68, 0.53, 0.92, 0.99);
        p3_Px->SetTopMargin(0.05);
        p3_Px->SetBottomMargin(0.05);
        p3_Px->SetLeftMargin(0.05);
        p3_Px->SetRightMargin(0.05);
        p3_Px->SetFillStyle(4000);
        p3_Px->Draw();



        ///////////////////////////////////////////////////////////////////////
        // PAD 0 - TOP
        ///////////////////////////////////////////////////////////////////////


        p0_Px->cd();


        h_stack_2nubb->SetMaximum(PAD_U_Y_MAX_Px);
        h_stack_2nubb->SetMinimum(PAD_U_Y_MIN_Px);

        h_stack_2nubb->GetXaxis()->SetLabelFont(43);
        h_stack_2nubb->GetXaxis()->SetLabelSize(0 * 15);
        h_stack_2nubb->GetXaxis()->SetTitleFont(43);
        h_stack_2nubb->GetXaxis()->SetTitleSize(0);
        h_stack_2nubb->GetXaxis()->SetTitleOffset(1.0);
        h_stack_2nubb->GetXaxis()->SetTitle(0);

        h_stack_2nubb->GetYaxis()->SetLabelFont(43);
        h_stack_2nubb->GetYaxis()->SetLabelSize(0 * 15);
        h_stack_2nubb->GetYaxis()->SetTitleFont(43);
        h_stack_2nubb->GetYaxis()->SetTitleSize(20);
        h_stack_2nubb->GetYaxis()->SetTitleOffset(1.0);
        h_stack_2nubb->GetYaxis()->SetTitle("Events / 0.1 MeV");

        h_stack_2nubb->GetXaxis()->SetTickSize(0.04);
        h_stack_2nubb->GetYaxis()->SetTickSize(0);



        h_stack_2nubb->Draw("hist");
        h_stack_tl208_int->Draw("histsame");
        h_stack_bi214_int->Draw("histsame");
        h_stack_bi207_int->Draw("histsame");
        h_stack_internal->Draw("histsame");
        h_stack_neighbour->Draw("histsame");
        h_stack_radon->Draw("histsame");
        h_stack_external->Draw("histsame");
        h_stack_total_MC->Draw("histsame");

        h_data_fakedata->Draw("Esame");



    /*
        // TODO: does this work
//        stacks1D_major_Px->Draw("histsame"); // this didn't used to be here
        stacks1D_major_Px->SetMaximum(1000.0);
        stacks1D_major_Px->SetMinimum(0.0);
        stacks1D_major_Px->Draw("hist"); // this didn't used to be here
    */
  
// moved above
//        stacks1D_major_Px->SetMaximum(PAD_U_Y_MAX_Px);
//        stacks1D_major_Px->SetMinimum(PAD_U_Y_MIN_Px);
        //stacks1D_major_Px->Draw("hist"); // this used to be uncommented
//        stacks1D_major_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px); // this used to be commented
        //stacks1D_major_Px->Draw("hist");
        //p0->Update();
    //    stacks1D_major_Px->SetMaximum(PAD_U_Y_MAX);
    //    stacks1D_major_Px->SetMinimum(PAD_U_Y_MIN);
    //    stacks1D_major_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN, PAD_U_Y_MAX);

#if 0
        stacks1D_major_Px->SetTitle("");

//        stacks1D_major_Px->GetYaxis()->SetTickSize(0.0);
        //stacks1D_major_Px->GetXaxis()->SetTickSize(0.0);

        stacks1D_major_Px->GetXaxis()->SetLabelFont(43);
        stacks1D_major_Px->GetXaxis()->SetLabelSize(0 * 15);
        stacks1D_major_Px->GetXaxis()->SetTitleFont(43);
        stacks1D_major_Px->GetXaxis()->SetTitleSize(0);
        stacks1D_major_Px->GetXaxis()->SetTitleOffset(1.0);
        stacks1D_major_Px->GetXaxis()->SetTitle("");

        stacks1D_major_Px->GetYaxis()->SetLabelFont(43);
        stacks1D_major_Px->GetYaxis()->SetTitleFont(43);
        stacks1D_major_Px->GetYaxis()->SetTitleSize(20);
        stacks1D_major_Px->GetYaxis()->SetTitleOffset(1.0);
        stacks1D_major_Px->GetYaxis()->SetTitle("Events / 0.1 MeV");

        if(AXISMODE == 2)
        {
            stacks1D_major_Px->GetYaxis()->SetLabelSize(0.0);
            stacks1D_major_Px->GetYaxis()->SetTickSize(0.0);
        }
        else if(AXISMODE == 1)
        {
            stacks1D_major_Px->GetYaxis()->SetLabelSize(15);
            stacks1D_major_Px->GetYaxis()->SetTickLength(1.0 / 213.0);
        }
        //stacks1D_major_Px->GetXaxis()->SetTickSize(0.05);
#endif

    /*
        hAllMC1D_Px->GetXaxis()->SetLabelFont(43);
        hAllMC1D_Px->GetXaxis()->SetLabelSize(0 * 15);
        hAllMC1D_Px->GetXaxis()->SetTitleFont(43);
        hAllMC1D_Px->GetXaxis()->SetTitleSize(0);
        hAllMC1D_Px->GetXaxis()->SetTitleOffset(1.0);
        hAllMC1D_Px->GetXaxis()->SetTitle("");

        hAllMC1D_Px->GetYaxis()->SetLabelFont(43);
        hAllMC1D_Px->GetYaxis()->SetTitleFont(43);
        hAllMC1D_Px->GetYaxis()->SetTitleSize(20);
        hAllMC1D_Px->GetYaxis()->SetTitleOffset(1.0);
        hAllMC1D_Px->GetYaxis()->SetTitle("Events / 0.1 MeV");

        if(AXISMODE == 2)
        {
            hAllMC1D_Px->GetYaxis()->SetLabelSize(0.0);
            hAllMC1D_Px->GetYaxis()->SetTickSize(0.0);
        }
        else if(AXISMODE == 1)
        {
            hAllMC1D_Px->GetYaxis()->SetLabelSize(15);
            hAllMC1D_Px->GetYaxis()->SetTickLength(1.0 / 213.0);
        }
    */
    /*
        hAllMC1D_Px->SetLineWidth(2);
        hAllMC1D_Px->SetLineColor(kBlack);
        hAllMC1D_Px->SetFillColor(kWhite);
        hAllMC1D_Px->SetFillStyle(0);
        //hAllMC1D_Px->Sumw2();
        //hAllMC1D_Px->Draw("hist sames");
    */
        TString Nmc_Px_str;
        //Nmc_Px_str.Form("%i", (int)hAllMC1D_Px->Integral()); // TODO: float?
        Nmc_Px_str.Form("%i", (int)h_stack_total_MC->Integral()); // TODO: float?
        //hAllMC1D_Px->SetTitle("Total MC (" + Nmc_Px_str + ")");
        h_stack_total_MC->SetTitle("Total MC (" + Nmc_Px_str + ")");
    /*
std::cout << "hALLMC1D->Draw()" << std::endl;
////////        hAllMC1D_Px->Draw("hist same");
        hAllMC1D_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px);
        if(g_mode_fake_data == false)
        {
            data1D_Px->SetLineWidth(2);
            //data1D_Px->SetMarkerStyle(20);
            //data1D_Px->SetMarkerSize(1.0);
            data1D_Px->SetLineColor(kBlack); // TODO: not needed? I forget reason for adding
            //data1D_Px->SetMarkerColor(kBlack); // TODO
            data1D_Px->SetFillColor(kBlack); // TODO
        }
        if(g_mode_fake_data == true)
        {
            fakeData1D_Px->SetLineWidth(2);
            //fakeData1D_Px->SetMarkerStyle(20);
            //fakeData1D_Px->SetMarkerSize(1.0);
            fakeData1D_Px->SetLineColor(kBlack);
            //fakeData1D_Px->SetMarkerColor(kBlack);
            fakeData1D_Px->SetFillColor(kBlack);
        }
    */
        //TString Ndata_Px_str;
        TString Ndata_fakedata_str;
        if(g_mode_fake_data == false)
        {
            //Ndata_Px_str.Form("%i", (int)data1D_Px->Integral()); // TODO: float?
            //data1D_Px->SetTitle("Data (" + Ndata_Px_str + ")");
            Ndata_fakedata_str.Form("%i", (int)h_data_fakedata->Integral()); // TODO: float?
            h_data_fakedata->SetTitle("Data (" + Ndata_fakedata_str + ")");
        }
        //TString Nfakedata_Px_str;
        if(g_mode_fake_data == true)
        {
            //Nfakedata_Px_str.Form("%i", (int)fakeData1D_Px->Integral()); // TODO: float?
            //fakeData1D_Px->SetTitle("Fake Data (" + Nfakedata_Px_str + ")");
            Ndata_fakedata_str.Form("%i", (int)h_data_fakedata->Integral()); // TODO: float?
            h_data_fakedata->SetTitle("Fake Data (" + Ndata_fakedata_str + ")");
        }
    /*
        if(g_mode_fake_data == false)
        {
            //data1D_Px[i]->Draw("PEsames");
            //data1D_Px->Draw("PEsame");
////////            data1D_Px->Draw("Esame");
            data1D_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px);
        }
        if(g_mode_fake_data == true)
        {
            //fakeData1D_Px->Draw("PEsame");
////////            fakeData1D_Px->Draw("Esame");
            fakeData1D_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px);
        }
        //data1D_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN, PAD_U_Y_MAX); // TODO???
    */
//std::cout << "data draw" << std::endl;
        //double chi2;
        int nch_Px = -1;
        nch_Px = get_nch_1D(h_stack_total_MC, h_data_fakedata);
        #if 0
        if(g_mode_fake_data == false)
        {
            //nch_Px = get_nch_1D(hAllMC1D_Px, data1D_Px);
            nch_Px = get_nch_1D(h_stack_total_MC, h_data_fakedata);
            /*
            if(phase_arg_str == "P1")
            {
                ndf_Px -= drawinputdata.nfreeparam_P1; 
            }
            else if(phase_arg_str == "P2")
            {
                ndf_Px -= drawinputdata.nfreeparam_P2;
            }
            else if(phase_arg_str == "PBOTH")
            {
                ndf_Px -= drawinputdata.nfreeparam;
            }
            */
        }
        if(g_mode_fake_data == true)
        {
            //nch_Px = get_nch_1D(hAllMC1D_Px, fakeData1D_Px);
            nch_Px = get_nch_1D(h_stack_total_MC, h_data_fakedata);
        }
        #endif
        //nch_Px -= drawinputdata.nfreeparam;
        //int igood;

        // TODO: should chisquare value include the constraints? because at
        // the moment it does not

        // TODO: chi2 value is different from fit_2e code
        //double prob = data1D->Chi2TestX(hAllMC1D, chi2, ndf, igood, "UW");
        //std::cout << "1: prob=" << prob << " chi2=" << chi2 << " igood=" << igood << " ndf=" << ndf << std::endl;
        /*
        double prob = data1D->Chi2TestX(hAllMC1D, chi2, ndf, igood, "UU");
        std::cout << "2: prob=" << prob << " chi2=" << chi2 << " igood=" << igood << " ndf=" << ndf << std::endl;
        */
        //prob = data1D->Chi2TestX(hAllMC1D, chi2, ndf, igood, "WU");
        //std::cout << "3: prob=" << prob << " chi2=" << chi2 << " igood=" << igood << " ndf=" << ndf << std::endl;
        //prob = data1D->Chi2TestX(hAllMC1D, chi2, ndf, igood, "WW");
        //std::cout << "4: prob=" << prob << " chi2=" << chi2 << " igood=" << igood << " ndf=" << ndf << std::endl;

        /*
        double mychi2 = 0.0;
        for(int i = 1; i <= data1D->GetNbinsX(); ++ i)
        {
            double content1 = data1D->GetBinContent(i);
            if(content1 <= 0.0) continue;
            double content2 = hAllMC1D->GetBinContent(i);
            double error1 = data1D->GetBinError(i);
            double error2 = hAllMC1D->GetBinError(i);
            //std::cout << "i=" << i << " " << content1 << " " << content2 << " " << error1 << " " << error2 << std::endl;
            mychi2 += std::pow((content1 - content2) / error1, 2.0);
        }
        std::cout << "mychi2=" << mychi2 << std::endl;
        */

    //    std::cout << "ROOT: chi2=" << chi2

        // TODO: check if I can get fcn value from the minuit fit object
        /*
        chi2_str.Form("%4.3f", chi2);
        */
        double fval_Px = 0.0;
        /*
        if(phase_arg_str == "P1")
        {
            fval_Px = drawinputdata.chi2_P1;
        }
        else if(phase_arg_str == "P2")
        {
            fval_Px = drawinputdata.chi2_P2;
        }
        else if(phase_arg_str == "PBOTH")
        {
            fval_Px = drawinputdata.chi2;
        }
        */
        if(g_mode_fake_data == false)
        {
            //fval_Px = calc_chi2_draw(data1D_Px, hAllMC1D_Px);
            //fval_Px = calc_chi2_draw_metric_CHISQ_ERR_SQRT_MC(data1D_Px, hAllMC1D_Px);
            fval_Px = calc_chi2_draw_metric_CHISQ_ERR_SQRT_MC(h_data_fakedata, h_stack_total_MC);
        }
        if(g_mode_fake_data == true)
        {
            //fval_Px = calc_chi2_draw(fakeData1D_Px, hAllMC1D_Px);
            //fval_Px = calc_chi2_draw_metric_CHISQ_ERR_SQRT_MC(fakeData1D_Px, hAllMC1D_Px);
            fval_Px = calc_chi2_draw_metric_CHISQ_ERR_SQRT_MC(h_data_fakedata, h_stack_total_MC);
        }
        TString fval_Px_str;
        fval_Px_str.Form("%.1f", fval_Px);
        int ndf_Px = nch_Px - drawinputdata.nfp;
        TString ndf_Px_str;
        ndf_Px_str.Form("%i", ndf_Px);
        /*
        mychi2_str.Form("%4.3f", mychi2);
        */

        
        TGaxis *axis = nullptr;
        TGaxis *axis2 = nullptr;
        if(AXISMODE == 2)
        {
            axis = new TGaxis(0.0, PAD_U_Y_MIN_Px,
                                      0.0, PAD_U_Y_MAX_Px,
                                      PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px,
                                      510, "");
            axis->SetLabelFont(43);
            axis->SetLabelSize(15);
            //axis->SetTickLength(0.05 / (0.7 * 0.9)); // TODO: think the margin is 0.1
            //axis->SetTickSize(0.05 / (0.7 * 0.9));
            axis->ChangeLabel(1, -1, 0);
            axis->Draw();

            axis2 = new TGaxis(5.0, PAD_U_Y_MIN_Px,
                                       5.0, PAD_U_Y_MAX_Px,
                                       PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px,
                                       510, "+");
            axis2->SetLabelFont(43);
            axis2->SetLabelSize(0);
            //axis2->SetTickLength(0.05 / (0.3 * 0.6)); // TODO: the margin is 0.4
            //axis2->SetTickSize(0.05 / (0.3 * 0.6));
            axis2->ChangeLabel(1, -1, 0);
            axis2->Draw();
        }
//std::cout << "axis draw" << std::endl;


        double xi_31_Px = 0.0;
        if(g_pg.get_xi_31_int_param_number() != -1)
        {
            int xi_31_ext_param_number = g_pg.get_xi_31_ext_param_number();
            /*double*/ xi_31_Px = params.at(xi_31_ext_param_number);
        }
        TString xilatexstr;
        xilatexstr.Form("#xi_{31}^{2#nu#beta#beta} = %.3f", xi_31_Px);
        TLatex xilatex;
        xilatex.SetNDC();
        xilatex.SetTextFont(63);
        xilatex.SetTextSize(20);
        xilatex.DrawLatex(0.45, 0.70, xilatexstr);

        double chi2_global = drawinputdata.chi2;
        int nch_global = drawinputdata.nch;
        int nfp_global = drawinputdata.nfp;
        int ndf_global = nch_global - nfp_global;
        TString chi2_global_str;
        chi2_global_str.Form("%.1f", chi2_global);
        TString ndf_global_str;
        ndf_global_str.Form("%i", ndf_global);

        TString chilatexstr;
        chilatexstr.Form("#frac{#chi^{2}_{global}}{ndf} = #frac{%s}{%s}", chi2_global_str.Data(), ndf_global_str.Data());

        TLatex chilatex;
        chilatex.SetNDC();
        chilatex.SetTextFont(63);
        chilatex.SetTextSize(20);
        chilatex.DrawLatex(0.45, 0.45, chilatexstr);






        ///////////////////////////////////////////////////////////////////////
        // PAD 2 - MIDDLE
        ///////////////////////////////////////////////////////////////////////





        // second pad
        p1_Px->cd();


        h_ratio->SetMaximum(PAD_L_Y_MAX_Px);
        h_ratio->SetMinimum(PAD_L_Y_MIN_Px);

        h_ratio->GetXaxis()->SetLabelFont(43);
        h_ratio->GetXaxis()->SetLabelSize(0 * 15);
        h_ratio->GetXaxis()->SetTitleFont(43);
        h_ratio->GetXaxis()->SetTitleSize(0);
        h_ratio->GetXaxis()->SetTitleOffset(1.0);
        h_ratio->GetXaxis()->SetTitle(0);

        h_ratio->GetYaxis()->SetLabelFont(43);
        h_ratio->GetYaxis()->SetLabelSize(0 * 15);
        h_ratio->GetYaxis()->SetTitleFont(43);
        h_ratio->GetYaxis()->SetTitleSize(20);
        h_ratio->GetYaxis()->SetTitleOffset(1.0);
        h_ratio->GetYaxis()->SetTitle("data / MC");

        h_ratio->GetXaxis()->SetTickSize(0.06);
        h_ratio->GetYaxis()->SetTickSize(0);



        //hRatio_Px->SetMarkerStyle(20);
        //hRatio_Px->SetMarkerSize(1.0);
        h_ratio->SetLineWidth(2);
        h_ratio->Draw("axis");
        TLine *linezero = new TLine(0.0, 1.0, 5.0, 1.0);
        linezero->SetLineWidth(2);
        linezero->SetLineColorAlpha(kGray + 2, 0.5);
        linezero->Draw();
        //hRatio_Px->Draw("EPsame");
        h_ratio->Draw("Esame");
        //TLine *zeroline = new TLine(0.0, 0.0, 5.0, 0.0);
        //zeroline->Draw();

        std::vector<TLine*> outofbounds_line;
        std::vector<TMarker*> outofbounds_marker;
        for(Int_t i = 1; i <= h_ratio->GetNbinsX(); ++ i)
        {
            Double_t content = h_ratio->GetBinContent(i);
            Double_t error = h_ratio->GetBinError(i);
            Double_t center = h_ratio->GetBinCenter(i);
            if(content == 0.0)
            {
                if(error == 0.0)
                {
                    continue;
                }
            }
            if(content > PAD_L_Y_MAX_Px)
            {
                if(content - error >= PAD_L_Y_MAX_Px)
                {
                    Double_t y_center = 0.5 * (PAD_L_Y_MAX_Px + PAD_L_Y_MIN_Px);
                    Double_t y_delta = 0.5 * (PAD_L_Y_MAX_Px - PAD_L_Y_MIN_Px);
                    Double_t y_value = y_center + 0.95 * y_delta;
                    TMarker *m = new TMarker(center, y_value, 22);
                    m->SetMarkerColor(kBlack);
                    outofbounds_marker.push_back(m);
                    m->Draw();

                    //std::cout << "drawing marker: " << center << " " << y_center << std::endl;
                }
                else
                {
                    TLine *l = new TLine(center, content - error, center, PAD_L_Y_MAX_Px);
                    l->SetLineWidth(2.0);
                    l->SetLineColor(kBlack);
                    outofbounds_line.push_back(l);
                    l->Draw();
                }
            }
            else if(content < PAD_L_Y_MIN_Px)
            {
                if(content + error <= PAD_L_Y_MAX_Px)
                {
                    Double_t y_center = 0.5 * (PAD_L_Y_MAX_Px + PAD_L_Y_MIN_Px);
                    Double_t y_delta = 0.5 * (PAD_L_Y_MAX_Px - PAD_L_Y_MIN_Px);
                    Double_t y_value = y_center - 0.95 * y_delta;
                    TMarker *m = new TMarker(center, y_value, 23);
                    m->SetMarkerColor(kBlack);
                    outofbounds_marker.push_back(m);
                    m->Draw();

                    //std::cout << "drawing marker: " << center << " " << y_center << std::endl;
                }
                else
                {
                    TLine *l = new TLine(center, content + error, center, PAD_L_Y_MIN_Px);
                    l->SetLineWidth(2.0);
                    l->SetLineColor(kBlack);
                    outofbounds_line.push_back(l);
                    l->Draw();
                }
            }
        }


        TGaxis *axis3 = nullptr;
        TGaxis *axis4 = nullptr;
        if(AXISMODE == 2)
        {
            //axis3 = new TGaxis(0.0, PAD_L_Y_MIN_Px, 0.0, PAD_L_Y_MAX_Px, PAD_L_Y_MIN_Px, PAD_L_Y_MAX_Px, 503, "");
            axis3 = new TGaxis(0.0, PAD_L_Y_MIN_Px, 0.0, PAD_L_Y_MAX_Px, PAD_L_Y_MIN_Px, PAD_L_Y_MAX_Px, 505, "");
            axis3->SetLabelFont(43);
            axis3->SetLabelSize(15);
            axis3->ChangeLabel(1, -1, 0);
            axis3->ChangeLabel(-1, -1, 0);
            axis3->Draw();

            //axis4 = new TGaxis(5.0, PAD_L_Y_MIN_Px, 5.0, PAD_L_Y_MAX_Px, PAD_L_Y_MIN_Px, PAD_L_Y_MAX_Px, 503, "+");
            axis4 = new TGaxis(5.0, PAD_L_Y_MIN_Px, 5.0, PAD_L_Y_MAX_Px, PAD_L_Y_MIN_Px, PAD_L_Y_MAX_Px, 505, "+");
            axis4->SetLabelFont(43);
            axis4->SetLabelSize(0);
            axis4->ChangeLabel(1, -1, 0);
            axis4->ChangeLabel(-1, -1, 0);
            axis4->Draw();
        }

        ///////////////////////////////////////////////////////////////////////
        // PAD2 - BOTTOM
        ///////////////////////////////////////////////////////////////////////


        // third pad
        p2_Px->cd();


        h_pull->SetMaximum(PAD_LL_Y_MAX_Px);
        h_pull->SetMinimum(PAD_LL_Y_MIN_Px);

        h_pull->GetXaxis()->SetLabelFont(43);
        h_pull->GetXaxis()->SetLabelSize(15);
        h_pull->GetXaxis()->SetTitleFont(43);
        h_pull->GetXaxis()->SetTitleSize(20);
        h_pull->GetXaxis()->SetTitleOffset(3.0);
        //h_pull->GetXaxis()->SetTitle(0);

        h_pull->GetYaxis()->SetLabelFont(43);
        h_pull->GetYaxis()->SetLabelSize(0 * 15);
        h_pull->GetYaxis()->SetTitleFont(43);
        h_pull->GetYaxis()->SetTitleSize(20);
        h_pull->GetYaxis()->SetTitleOffset(1.0);
        h_pull->GetYaxis()->SetTitle("#chi^{2} pull (#sigma) ");

        h_pull->GetXaxis()->SetTickSize(0.05);
        h_pull->GetYaxis()->SetTickSize(0);

        //hSig_Px->SetMarkerStyle(20);
        //hSig_Px->SetMarkerSize(1.0);
        h_pull->Draw("axis");
        //std::cout << "hSig_Px: "
        //          << hSig_Px->GetMaximum() << " "
        //          << hSig_Px->GetMinimum() << std::endl;
        TLine *linezero2 = new TLine(0.0, 0.0, 5.0, 0.0);
        linezero2->SetLineWidth(2);
        linezero2->SetLineColorAlpha(kGray + 2, 0.5);
        linezero2->Draw();
        h_pull->Draw("histsame");

#if 0
        bool draw_ok = false;
        for(Int_t i = 0; i < hSig_Px->GetNbinsX(); ++ i)
        {
            if(hSig_Px->GetBinContent(i) > 0.0)
            {
                draw_ok = true;
                break;
            }
        }
        if(draw_ok == true)
        {
        hSig_Px->Draw("histsame");
        }
        
        //TODO fillstyle
std::cout << "sigdraw" << std::endl;
#endif

//        std::cout << "PAD_LL_Y_MIN_Px=" << PAD_LL_Y_MIN_Px 
//                  << "PAD_LL_Y_MAX_Px=" << PAD_LL_Y_MAX_Px << std::endl;
        TGaxis *axis5 = nullptr;
        TGaxis *axis6 = nullptr;
        if(AXISMODE == 2 && true)
        {
            //TGaxis *axis5 = new TGaxis(0.0, PAD_LL_Y_MIN_Px, 0.0, PAD_LL_Y_MAX_Px - 0.0001, PAD_LL_Y_MIN_Px, PAD_LL_Y_MAX_Px - 0.0001, 503, "");
            //TGaxis *axis5 = new TGaxis(0.0, PAD_LL_Y_MIN_Px, 0.0 + 0.0001, PAD_LL_Y_MAX_Px, PAD_LL_Y_MIN_Px, PAD_LL_Y_MAX_Px, 503, "");
            axis5 = new TGaxis(0.0, PAD_LL_Y_MIN_Px, 0.0, PAD_LL_Y_MAX_Px, PAD_LL_Y_MIN_Px, PAD_LL_Y_MAX_Px, 505, "");

            axis5->SetLabelFont(43);
            axis5->SetLabelSize(15);
            axis5->ChangeLabel(-1, 0, 0);
            //axis5->ChangeLabel(1, 0, 0);
            axis5->Draw();
        }

        if(AXISMODE == 2 && true)
        {
            //TGaxis *axis6 = new TGaxis(5.0, PAD_LL_Y_MIN_Px, 5.0, PAD_LL_Y_MAX_Px - 0.0001, PAD_LL_Y_MIN_Px, PAD_LL_Y_MAX_Px - 0.0001, 503, "+");
            //axis6 = new TGaxis(5.0, PAD_LL_Y_MIN_Px, 5.0, PAD_LL_Y_MAX_Px, PAD_LL_Y_MIN_Px, PAD_LL_Y_MAX_Px, 503, "+");
            axis6 = new TGaxis(5.0, PAD_LL_Y_MIN_Px, 5.0, PAD_LL_Y_MAX_Px, PAD_LL_Y_MIN_Px, PAD_LL_Y_MAX_Px, 505, "+");
            axis6->SetLabelFont(43);
            axis6->SetLabelSize(0);
            axis6->ChangeLabel(-1, 0, 0);
            //axis6->ChangeLabel(1, 0, 0);
            axis6->Draw();
        }

        Double_t deltapos = 0.0;
        Double_t deltaneg = 0.0;
        for(Int_t i = 0; i <= h_pull->GetNbinsX(); ++ i)
        {
            Double_t content = h_pull->GetBinContent(i);
            if(content >= 0.0)
            {
                deltapos += std::pow(content, 2.0);
            }
            else
            {
                deltaneg += std::pow(content, 2.0);
            }
        }
        TString deltaposlatexstr;
        deltaposlatexstr.Form("#chi^{2}_{+} = %1.1f", deltapos);
        TLatex deltaposlatex;
        deltaposlatex.SetNDC();
        deltaposlatex.SetTextFont(63);
        deltaposlatex.SetTextSize(20);
        //deltaposlatex.SetTextAlign(31);
        deltaposlatex.DrawLatex(0.80, 0.75, deltaposlatexstr);
        TString deltaneglatexstr;
        deltaneglatexstr.Form("#chi^{2}_{-} = %1.1f", deltaneg);
        TLatex deltaneglatex;
        deltaneglatex.SetNDC();
        deltaneglatex.SetTextFont(63);
        deltaneglatex.SetTextSize(20);
        //deltaneglatex.SetTextAlign(31);
        deltaneglatex.DrawLatex(0.80, 0.45, deltaneglatexstr);
//        std::cout << "deltapos=" << deltapos << " deltaneg=" << deltaneg << std::endl;





        ///////////////////////////////////////////////////////////////////////
        // PAD3 - ADDITIONAL PAD FOR LEGEND
        ///////////////////////////////////////////////////////////////////////


        p3_Px->cd();
        TLegend *leg = new TLegend(0.02, 0.02, 0.95, 0.95);
        if(g_mode_fake_data == false)
        {
            //leg->AddEntry(data1D_Px, "Data (" + Ndata_Px_str + ")", "PE"); // TODO PEL ??? works?
            //leg->AddEntry(data1D_Px, "Data (" + Ndata_Px_str + ")", "EL"); // TODO PEL ??? works?
            //leg->AddEntry(h_data_fakedata, "Data (" + Ndata_Px_str + ")", "EL"); // TODO PEL ??? works?
            leg->AddEntry(h_data_fakedata, "Data (" + Ndata_fakedata_str + ")", "EL"); // TODO PEL ??? works?
        }
        if(g_mode_fake_data == true)
        {
            //leg->AddEntry(fakeData1D_Px, "Fake Data (" + Nfakedata_Px_str + ")", "PE");
            //leg->AddEntry(fakeData1D_Px, "Fake Data (" + Nfakedata_Px_str + ")", "EL");
            //leg->AddEntry(h_data_fakedata, "Fake Data (" + Nfakedata_Px_str + ")", "EL");
            leg->AddEntry(h_data_fakedata, "Fake Data (" + Ndata_fakedata_str + ")", "EL");
        }
        leg->AddEntry(h_stack_total_MC, "Total MC (" + Nmc_Px_str + ")", "L");
        leg->AddEntry(h_stack_2nubb, "2#nu#beta#beta", "F");
        leg->AddEntry(h_stack_tl208_int, "^{208}Tl ^{150}Nd Foil", "F");
        leg->AddEntry(h_stack_bi214_int, "^{214}Bi ^{150}Nd Foil", "F");
        leg->AddEntry(h_stack_bi207_int, "^{207}Bi ^{150}Nd Foil", "F");
        leg->AddEntry(h_stack_internal, "Internal Foil", "F");
        leg->AddEntry(h_stack_neighbour, "Neighbour Foil", "F");
        leg->AddEntry(h_stack_radon, "Radon", "F");
        leg->AddEntry(h_stack_external, "External", "F");
        //leg->AddEntry((TObject*)nullptr, "#chi^{2}/ndf=" + chi2_str + "/" + ndf_str, "");
        //leg->AddEntry((TObject*)nullptr, "fval #chi^{2}/ndf=" + fval_str + "/" + ndf_str, "");
        leg->AddEntry((TObject*)nullptr, "#chi^{2}/ndf=" + fval_Px_str + "/" + ndf_Px_str, "");
        //leg->AddEntry((TObject*)nullptr, "my #chi^{2}/ndf=" + mychi2_str + "/" + ndf_str, "");
        //leg->AddEntry(h_other, "other", "f");
        leg->SetBorderSize(0);
        leg->SetFillColor(0);
        //leg->SetTextFont(62);
        leg->SetTextFont(63);
        //leg->SetTextSize(0.035);
        leg->SetTextSize(15);
        leg->SetShadowColor(kGray + 2);
        leg->SetBorderSize(5);
        leg->Draw("BR");
        p0_Px->cd();

//std::cout << "leg draw" << std::endl;
        //TLatex latexlabel;
        //latexlabel.SetNDC();
        //latexlabel.SetTextFont(62);
        //latexlabel.SetTextSize(0.035);
        //latexlabel.DrawLatex(0.63, 0.23, "#frac{#chi^{2}}{ndf} = #frac{" + chi2_str + "}{" + ndf_str + "}");




//        canvas_Px->Show();

        // TODO: other parameters in this class for filename/dir
        std::string saveas_filename = drawinputdata.saveas_filename;
        std::string saveas_dir;
        bool mode_parallel = drawinputdata.mode_parallel;
        if(mode_parallel == false)
        {
            saveas_dir = drawinputdata.serial_dir;
        }
        else if(mode_parallel == true)
        {
            saveas_dir = drawinputdata.parallel_dir;
        }
        if(saveas_filename != "NOSAVE")
        {
            // NOTE: do not delete may be useful later
            //std::string base_name;
            //std::string extension;
            //filename_split_extension(saveas_filename, base_name, extension);
            std::vector<std::string> extensions;
            if(drawinputdata.saveas_png == true)
            {
                extensions.push_back(".png");
            }
            if(drawinputdata.saveas_pdf == true)
            {
                extensions.push_back(".pdf");
            } 
            if(drawinputdata.saveas_eps == true)
            {
                extensions.push_back(".eps");
            } 
            if(drawinputdata.saveas_C == true)
            {
                extensions.push_back(".C");
            } 

            //std::string dir = base_name + "_c" + "_" + std::string(channel_str);
            std::string dir = saveas_dir;
            for(std::vector<std::string>::iterator it{extensions.begin()};
                it != extensions.end(); ++ it)
            {
                std::string extension = *it;
                std::string base_name = saveas_filename;
                if(mode_parallel == true)
                {
                    int number_job_id = drawinputdata.parallel_job_id;
                    std::string name = base_name + "_JID" + std::to_string(number_job_id) + "_c" + "_" + std::string(channel_str) + "_" + phase_arg_str + extension;
                    canvas_saveas_helper(dir, name, canvas_Px);
                }
                else if(mode_parallel == false)
                {
                    std::string name = base_name + "_c" + "_" + std::string(channel_str) + "_" + phase_arg_str + extension;
                    //std::cout << "saving " << name << std::endl;
                    canvas_saveas_helper(dir, name, canvas_Px);
                }
            }
        }

//std::cout << "done with save" << std::endl;

        draw_output_data drawoutputdata_Px;
        //drawoutputdata_Px.stacks1D = stacks1D_Px;
        //drawoutputdata_Px.h_2nubb = h_2nubb_Px;
        //drawoutputdata_Px.h_tl208_int = h_tl208_int_Px;
        //drawoutputdata_Px.h_bi214_int = h_bi214_int_Px;
        //drawoutputdata_Px.h_bi207_int = h_bi207_int_Px;
        //drawoutputdata_Px.h_internal = h_internal_Px;
        //drawoutputdata_Px.h_external = h_external_Px;
        //drawoutputdata_Px.h_radon = h_radon_Px;
        //drawoutputdata_Px.h_neighbour = h_neighbour_Px;
        //drawoutputdata_Px.h_other = h_other_Px;
        drawoutputdata_Px.canvas = canvas_Px;
        //drawoutputdata_Px.hRatio = hRatio_Px;
        //drawoutputdata_Px.hSig = hSig_Px;
        //drawoutputdata_Px.hAllMC1D = hAllMC1D_Px;
        //drawoutputdata_Px.data1D = data1D_Px;
        //drawoutputdata_Px.fakeData1D = fakeData1D_Px;
        if(phase_arg_str == "P1")
        {
            drawOutputDataMap[channel][1] = drawoutputdata_Px;
        }
        else if(phase_arg_str == "P2")
        {
            drawOutputDataMap[channel][2] = drawoutputdata_Px;
        }
        else if(phase_arg_str == "PBOTH")
        {
            drawOutputDataMap[channel][0] = drawoutputdata_Px;
        }
    }



}



#endif // NEWLOGLIKFITTER_DRAW_CHANNEL_PHASE_WITH_PULL_H
