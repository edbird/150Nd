#ifndef NEWLOGLIKFITTER_DRAWCHANNEL_H
#define NEWLOGLIKFITTER_DRAWCHANNEL_H












// how to choose which histograms are drawn
// (channel)
// use channel_enable_draw_1D / channel_enable_draw_2D
class draw_input_data
{
    public:

    draw_input_data()
    {
        draw_P1 = true;
        draw_P2 = true;
        draw_P1P2 = true;

        chi2 = -1.0;
        chi2_P1 = -1.0;
        chi2_P2 = -1.0;
        nch = -1.0;
        nch_P1 = -1.0;
        nch_P2 = -1.0;
        nfp = 0.0;
        //nfreeparam_P1 = 0.0;
        //nfreeparam_P2 = 0.0;

        mode_parallel = false;
        parallel_job_id = -1;
        parallel_dir = "parallel";
        serial_dir = ".";

        saveas_filename = "untitled";

        saveas_png = true;
        saveas_pdf = false;
        saveas_eps = false;
        saveas_C = false;
    }

    bool draw_P1;
    bool draw_P2;
    bool draw_P1P2;

    // chi2 (global)
    double chi2;
    double chi2_P1;
    double chi2_P2;
    // ndf (global)
    //double ndf;
    //double ndf_P1;
    //double ndf_P2;
    int nch;
    int nch_P1;
    int nch_P2; // number of channel = number of bin
    // number of free parameters (global)
    int nfp;
    //double nfreeparam_P1;
    //double nfreeparam_P2;
 
    // obtain from global variable
    //bool g_mode_fake_data;

    bool mode_parallel;
    int parallel_job_id;
    std::string parallel_dir;
    std::string serial_dir;

    std::string saveas_filename;
    bool saveas_png;
    bool saveas_pdf;
    bool saveas_eps;
    bool saveas_C;

};


class draw_output_data
{

    public:

    draw_output_data()
    {
        stacks1D = nullptr;
        h_2nubb = nullptr;
        h_tl208_int = nullptr;
        h_bi214_int = nullptr;
        h_bi207_int = nullptr;
        h_internal = nullptr;
        h_external = nullptr;
        h_radon = nullptr;
        h_neighbour = nullptr;
        h_other = nullptr;

        canvas = nullptr;
        hRatio = nullptr;
        hSig = nullptr;
        hAllMC1D = nullptr;
        data1D = nullptr;
        fakeData1D = nullptr;
    }

    THStack *stacks1D;
    TH1D *h_2nubb;
    TH1D *h_tl208_int;
    TH1D *h_bi214_int;
    TH1D *h_bi207_int;
    TH1D *h_internal;
    TH1D *h_external;
    TH1D *h_radon; 
    TH1D *h_neighbour;
    TH1D *h_other;


    TCanvas *canvas;
    //TPad *p0;
    //TPad *p1;
    TH1D *hRatio;
    TH1D *hSig;
    TH1D *hAllMC1D;
    TH1D *data1D;
    TH1D *fakeData1D;

};


// map channel to map phase to 
// channel index: 0 to number1DHists
// phase index: 0, 1, 2: 0=BOTH, 1=P1, 2=P2
//
std::map<int, std::map<int, draw_output_data>> drawOutputDataMap;



void draw_channel_phase(
    const int channel,
    const std::string& phase_arg_str,
    const draw_input_data &drawinputdata,
    const std::vector<double> &params,
    const std::vector<double> &param_errs
    );



void draw_channel_phase_with_pull(
    const int channel,
    const std::string& phase_arg_str,
    const draw_input_data &drawinputdata,
    const std::vector<double> &params,
    const std::vector<double> &param_errs
    );




void draw_channel(
    const int channel,
    const draw_input_data &drawinputdata,
    const std::vector<double> &params,
    const std::vector<double> &param_errs
    )
{

    if(drawinputdata.draw_P1 == true)
    {
        //std::cout << "draw P1" << std::endl;
        //draw_channel_phase(channel, "P1", drawinputdata, params, param_errs);
        draw_channel_phase_with_pull(channel, "P1", drawinputdata, params, param_errs);
    }

    if(drawinputdata.draw_P2 == true)
    {
        //std::cout << "draw P2" << std::endl;
        //draw_channel_phase(channel, "P2", drawinputdata, params, param_errs);
        draw_channel_phase_with_pull(channel, "P2", drawinputdata, params, param_errs);
    }

    if(drawinputdata.draw_P1P2 == true)
    {
        //std::cout << "draw PBOTH" << std::endl;
        //draw_channel_phase(channel, "PBOTH", drawinputdata, params, param_errs);
        draw_channel_phase_with_pull(channel, "PBOTH", drawinputdata, params, param_errs);
    }

    


}




// thePhase_flag:
// 0 = P1
// 1 = P2
// 2 = P1+P2
void draw_channel_phase(
    const int channel,
    const std::string& phase_arg_str,
    const draw_input_data &drawinputdata,
    const std::vector<double> &params,
    const std::vector<double> &param_errs
    )
{

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

        ///////////////////////////////////////////////////////////////////////////
        // Draw phase as specified by phase_arg_str
        // which can be one of:
        // P1, P2, PBOTH
        ///////////////////////////////////////////////////////////////////////////

        THStack *stacks1D_Px;
        TH1D *h_2nubb_Px = nullptr;
        TH1D *h_tl208_int_Px = nullptr;
        TH1D *h_bi214_int_Px = nullptr;
        TH1D *h_bi207_int_Px = nullptr;
        TH1D *h_internal_Px = nullptr;
        TH1D *h_external_Px = nullptr;
        TH1D *h_radon_Px = nullptr; 
        TH1D *h_neighbour_Px = nullptr;
        TH1D *h_other_Px = nullptr;


        TCanvas *canvas_Px = nullptr;
        TPad *p0_Px = nullptr;
        TPad *p1_Px = nullptr;
        TH1D *hRatio_Px = nullptr;
        TH1D *hAllMC1D_Px = nullptr;
        TH1D *data1D_Px = nullptr;
        TH1D *fakeData1D_Px = nullptr;

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
            if(g_mode_fake_data == false)
            {
                search_object_Px = histname + std::string(DataFile) + "_" + (*it);
            }
            else if(g_mode_fake_data == true)
            {
                search_object_Px = histname + std::string("fakedata") + "_" + (*it);
            }
            //std::cout << "search_object_Px=" << search_object_Px << std::endl;
            //std::string search_object_P2 = histname + std::string(DataFile) + "_P2";
            if(g_mode_fake_data == false)
            {
                if(data1D_Px == nullptr)
                {
                    TString name = TString(allDataSamples1D->FindObject(search_object_Px.c_str())->GetName()) + "_drawclone";
                    data1D_Px = (TH1D*)allDataSamples1D->FindObject(search_object_Px.c_str())->Clone(name);
                    //tmpDataHistP2 = (TH1D*)allDataSamples1D->FindObject(search_object_P2.c_str());
                }
                else
                {
                    data1D_Px->Add((TH1D*)allDataSamples1D->FindObject(search_object_Px.c_str()));
                }
            }
            else if(g_mode_fake_data == true)
            {
                if(fakeData1D_Px == nullptr)
                {
                    TString name = TString(allFakeDataSamples1D->FindObject(search_object_Px.c_str())->GetName()) + "_drawclone";
                    fakeData1D_Px = (TH1D*)allFakeDataSamples1D->FindObject(search_object_Px.c_str())->Clone(name);
                    //tmpDataHistP2 = (TH1D*)allFakeDataSamples1D->FindObject(search_object_P2.c_str());
                    // TODO: will not work if logLikelihood not called before
                    // because LL function calls function to construct fakedata
                }
                else
                {
                    fakeData1D_Px->Add((TH1D*)allFakeDataSamples1D->FindObject(search_object_Px.c_str()));
                }
            }
        }

        TString channel_str;
        channel_str.Form("ch%i", channel);

        canvas_Px = new TCanvas("c_" + channel_str + "_" + phase_arg_str, "c_" + channel_str + "_" + phase_arg_str);
        canvas_Px->SetFillColor(kWhite);

        
        stacks1D_Px = new THStack("stacks1D_" + channel_str + "_" + phase_arg_str, "stacks1D_" + channel_str + "_" + phase_arg_str);


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
                    TString name = TString(allMCSamples1D[channel]->FindObject(search_object_P1.c_str())->GetName()) + "_tmp";
                    tmpHist1D_P1 = (TH1D*)allMCSamples1D[channel]->FindObject(search_object_P1.c_str())->Clone(name); // need to clone to scale
                    
                    if(tmpHist1D_P1 == nullptr)
                    {
                        std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
                        throw "problem";
                    }


                    Double_t scale_factor_P1 = params.at(paramNumberInt);
                    tmpHist1D_P1->Scale(scale_factor_P1);
                    if(debuglevel >= 6)
                    {
                        if(TString(tmpHist1D_P1->GetName()).Contains("hSingleEnergy_"))
                        {
                            std::cout << "P1: name=" << tmpHist1D_P1->GetName() << std::endl;
                            std::cout << "scale_factor_P1=" << scale_factor_P1 << std::endl;
                            std::cout << "bin 10: " << tmpHist1D_P1->GetBinContent(10) << std::endl;
                        }
                    }


                    if(tmpHist1D_P1->Integral() > 0)
                    {
                        stacks1D_Px->Add((TH1D*)tmpHist1D_P1->Clone()); // probably no clone here? TODO
                        stacker_helper(tmpHist1D_P1, // TODO Clone names in this function and above
                                       h_2nubb_Px,
                                       h_tl208_int_Px,
                                       h_bi214_int_Px,
                                       h_bi207_int_Px,
                                       h_internal_Px,
                                       h_neighbour_Px,
                                       h_radon_Px,
                                       h_external_Px,
                                       h_other_Px);

                        if(hAllMC1D_Px == nullptr)
                        {
                            hAllMC1D_Px = (TH1D*)tmpHist1D_P1->Clone("Total MC");
                        }
                        else
                        {
                            hAllMC1D_Px->Add((TH1D*)tmpHist1D_P1);
                        }
                    }
                    else
                    {
                        //std::cout << "not adding to stack, Integral() <= 0: " << tmpHist_draw1D->GetName() << std::endl;
                    }
                }

                
                if(paramEnabledP2 == true)
                {
                    TString name = TString(allMCSamples1D[channel]->FindObject(search_object_P2.c_str())->GetName()) + "_tmp";
                    tmpHist1D_P2 = (TH1D*)allMCSamples1D[channel]->FindObject(search_object_P2.c_str())->Clone(name); // need clone to scale
                
                    if(tmpHist1D_P2 == nullptr)
                    {
                        std::cout << "ERROR: Could not find object " << search_object_P2 << std::endl;
                        throw "problem";
                    }


                    Double_t scale_factor_P2 = params.at(paramNumberInt);
                    tmpHist1D_P2->Scale(scale_factor_P2);
                    if(debuglevel >= 6)
                    {
                        if(TString(tmpHist1D_P2->GetName()).Contains("hSingleEnergy_"))
                        {
                            std::cout << "P2: name=" << tmpHist1D_P2->GetName() << std::endl;
                            std::cout << "scale_factor_P2=" << scale_factor_P2 << std::endl;
                            std::cout << "bin 10: " << tmpHist1D_P2->GetBinContent(10) << std::endl;
                        }
                    }

                    
                    if(tmpHist1D_P2->Integral() > 0)
                    {
                        stacks1D_Px->Add((TH1D*)tmpHist1D_P2->Clone());
                        stacker_helper(tmpHist1D_P2,
                                       h_2nubb_Px,
                                       h_tl208_int_Px,
                                       h_bi214_int_Px,
                                       h_bi207_int_Px,
                                       h_internal_Px,
                                       h_neighbour_Px,
                                       h_radon_Px,
                                       h_external_Px,
                                       h_other_Px);

                        if(hAllMC1D_Px == nullptr)
                        {
                            hAllMC1D_Px = (TH1D*)tmpHist1D_P2->Clone("Total MC");
                        }
                        else
                        {
                            hAllMC1D_Px->Add((TH1D*)tmpHist1D_P2);
                        }
                    }
                    else
                    {
                        //std::cout << "not adding to stack, Integral() <= 0: " << tmpHist_draw1D->GetName() << std::endl;
                    }
                }

            }
        }


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

        double PAD_U_Y_MIN_P1_P2 = 0.0;
        double PAD_U_Y_MAX_P1_P2 = 500.0;
        double PAD_L_Y_MIN_P1_P2 = 0.0;
        double PAD_L_Y_MAX_P1_P2 = 3.0;

        double PAD_U_Y_MIN_P1 = 0.0;
        double PAD_U_Y_MAX_P1 = 500.0;
        double PAD_L_Y_MIN_P1 = 0.0;
        double PAD_L_Y_MAX_P1 = 3.0;

        double PAD_U_Y_MIN_P2 = 0.0;
        double PAD_U_Y_MAX_P2 = 500.0;
        double PAD_L_Y_MIN_P2 = 0.0;
        double PAD_L_Y_MAX_P2 = 3.0;

        double PAD_U_Y_MIN_Px = 0.0;
        double PAD_U_Y_MAX_Px = 0.0;
        double PAD_L_Y_MIN_Px = 0.0;
        double PAD_L_Y_MAX_Px = 0.0;

        if(channel == 0)
        {
            PAD_U_Y_MAX_P1_P2 = 300.0;

            PAD_U_Y_MAX_P1 = 70;
            PAD_U_Y_MAX_P2 = 240;
        }
        else if(channel == 1)
        {
            PAD_U_Y_MAX_P1_P2 = 1000.0;

            PAD_U_Y_MAX_P1 = 200;
            PAD_U_Y_MAX_P2 = 750;
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

            PAD_U_Y_MAX_P1 = 180;
            PAD_U_Y_MAX_P2 = 700;
        }
        else if(channel == 4)
        {
            PAD_U_Y_MAX_P1_P2 = 300.0;

            PAD_U_Y_MAX_P1 = 70;
            PAD_U_Y_MAX_P2 = 250;
        }
        else if(channel == 5)
        {
            PAD_U_Y_MAX_P1_P2 = 550.0;

            PAD_U_Y_MAX_P1 = 140;
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
        }
        else if(phase_arg_str == "P2")
        {
            PAD_U_Y_MIN_Px = PAD_U_Y_MIN_P2;
            PAD_U_Y_MAX_Px = PAD_U_Y_MAX_P2;
            PAD_L_Y_MIN_Px = PAD_L_Y_MIN_P2;
            PAD_L_Y_MAX_Px = PAD_L_Y_MAX_P2;
        }
        else if(phase_arg_str == "PBOTH")
        {
            PAD_U_Y_MIN_Px = PAD_U_Y_MIN_P1_P2;
            PAD_U_Y_MAX_Px = PAD_U_Y_MAX_P1_P2;
            PAD_L_Y_MIN_Px = PAD_L_Y_MIN_P1_P2;
            PAD_L_Y_MAX_Px = PAD_L_Y_MAX_P1_P2;
        }


        //stacks1D_major[i]->Draw("hist");
        //stacks1D_major[i]->SetMaximum(PAD_U_Y_MAX);
        //stacks1D_major[i]->SetMinimum(PAD_U_Y_MIN);
    //    stacks1D_major->GetYaxis()->SetRangeUser(PAD_U_Y_MIN, PAD_U_Y_MAX);
    //    hAllMC1D->SetMaximum(PAD_U_Y_MAX);
    //    hAllMC1D->SetMinimum(PAD_U_Y_MIN);
        hAllMC1D_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px);
    //    data1D->SetMaximum(PAD_U_Y_MAX);
    //    data1D->SetMinimum(PAD_U_Y_MIN);
        if(g_mode_fake_data == false)
        {
            data1D_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px);
        }
        else if(g_mode_fake_data == true)
        {
            fakeData1D_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px);
        }

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

        if(channel == 0)
        {
            hRatio_Px->GetXaxis()->SetTitle("2e Electron Energy [MeV]");
        }
        else if(channel == 1)
        {
            hRatio_Px->GetXaxis()->SetTitle("Single Electron Energy [MeV]");
        }
        else if(channel == 2)
        {
            hRatio_Px->GetXaxis()->SetTitle("High Energy Electron [MeV]");
        }
        else if(channel == 3)
        {
            hRatio_Px->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
        }
        else if(channel == 4)
        {
            hRatio_Px->GetXaxis()->SetTitle("Energy Sum [MeV]");
        }
        else if(channel == 5)
        {
            hRatio_Px->GetXaxis()->SetTitle("Energy Diff [MeV]");
        }

        //hRatio[i]->SetMaximum(PAD_L_Y_MAX);
        //hRatio[i]->SetMinimum(PAD_L_Y_MIN);
        hRatio_Px->GetYaxis()->SetRangeUser(PAD_L_Y_MIN_Px, PAD_L_Y_MAX_Px);

        hRatio_Px->GetXaxis()->SetTickSize(0.1);

        //hRatio->GetXaxis()->SetTitle("Electron Energy [MeV]");
        hRatio_Px->GetXaxis()->SetTitleFont(43);
        hRatio_Px->GetXaxis()->SetTitleSize(20);
        hRatio_Px->GetXaxis()->SetTitleOffset(3.25);

        hRatio_Px->GetYaxis()->SetTitle("data / MC");
        hRatio_Px->GetYaxis()->SetTitleFont(43);
        hRatio_Px->GetYaxis()->SetTitleSize(20);
        hRatio_Px->GetYaxis()->SetTitleOffset(1.0);

        hRatio_Px->GetYaxis()->SetLabelFont(43);
        hRatio_Px->GetYaxis()->SetLabelSize(0.0 * 31);

        hRatio_Px->GetYaxis()->SetTickSize(0.0);

        hRatio_Px->GetXaxis()->SetLabelFont(43);
        hRatio_Px->GetXaxis()->SetLabelSize(15);

        // TODO: fix the error bar with for ratio so that it is the same
        // as the line width for the main histogram
        // and remove the horizontal lines from both

        // in code copying from, canvas alloc. here
        canvas_Px->cd();
        canvas_Px->SetBottomMargin(0.0);
        p0_Px = new TPad("pad0_" + channel_str + "_" + phase_arg_str, "pad0" + channel_str + "_" + phase_arg_str, 0.0, 0.3, 1.0, 1.0);
        p0_Px->SetBottomMargin(0.0);
        //p0->SetGridx(1);
        //p0->SetGridy(1);
        p0_Px->SetGrid(0, 0);
        p0_Px->SetTicks(2, 2);
        p0_Px->Draw();

        canvas_Px->cd();
        p1_Px = new TPad("pad1_" + channel_str + "_" + phase_arg_str, "pad1_" + channel_str + "_" + phase_arg_str, 0.0, 0.0, 1.0, 0.3);
        p1_Px->SetTopMargin(0.0);
        p1_Px->SetBottomMargin(0.4);
        //p1->SetGridx(1);
        //p1->SetGridy(1);
        p1_Px->SetGrid(0, 0);
        p1_Px->SetTicks(2, 2);
        p1_Px->Draw();

        p0_Px->cd();

        // draw regular pad1
        // leave alone for now
        // TODO: axis
        // copy code from other file
        hAllMC1D_Px->SetTitle("");

        // TODO

        stacks1D_major_Px->Draw("hist"); // this didn't used to be here

        stacks1D_major_Px->GetXaxis()->SetTickLength(0.05);
        stacks1D_major_Px->GetXaxis()->SetTickSize(0.05);
        
        stacks1D_major_Px->SetMaximum(PAD_U_Y_MAX_Px);
        stacks1D_major_Px->SetMinimum(PAD_U_Y_MIN_Px);
        //stacks1D_major_Px->Draw("hist"); // this used to be uncommented
        stacks1D_major_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px); // this used to be commented
        //stacks1D_major_Px->Draw("hist");
        //p0->Update();
    //    stacks1D_major_Px->SetMaximum(PAD_U_Y_MAX);
    //    stacks1D_major_Px->SetMinimum(PAD_U_Y_MIN);
    //    stacks1D_major_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN, PAD_U_Y_MAX);

        stacks1D_major_Px->SetTitle("");

        stacks1D_major_Px->GetYaxis()->SetTickSize(0.0);
        //stacks1D_major_Px->GetXaxis()->SetTickSize(0.0);

        stacks1D_major_Px->GetYaxis()->SetTitle("Events / 0.1 MeV");
        stacks1D_major_Px->GetYaxis()->SetTitleSize(20);
        stacks1D_major_Px->GetYaxis()->SetTitleFont(43);
        stacks1D_major_Px->GetYaxis()->SetTitleOffset(1.0);

        stacks1D_major_Px->GetYaxis()->SetLabelFont(43);
        stacks1D_major_Px->GetYaxis()->SetLabelSize(0.0 * 31);

        stacks1D_major_Px->GetXaxis()->SetTitle("");
        stacks1D_major_Px->GetXaxis()->SetTitleSize(0);
        stacks1D_major_Px->GetXaxis()->SetTitleFont(43);
        stacks1D_major_Px->GetXaxis()->SetTitleOffset(1.0);

        stacks1D_major_Px->GetXaxis()->SetLabelSize(0);
        stacks1D_major_Px->GetXaxis()->SetLabelFont(63);


        hAllMC1D_Px->SetLineWidth(2);
        hAllMC1D_Px->SetLineColor(kBlack);
        hAllMC1D_Px->SetFillColor(kWhite);
        hAllMC1D_Px->SetFillStyle(0);
        //hAllMC1D_Px->Sumw2();
        //hAllMC1D_Px->Draw("hist sames");
        TString Nmc_Px_str;
        Nmc_Px_str.Form("%i", (int)hAllMC1D_Px->Integral()); // TODO: float?
        hAllMC1D_Px->SetTitle("Total MC (" + Nmc_Px_str + ")");
        hAllMC1D_Px->Draw("hist same");
        hAllMC1D_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px);
        if(g_mode_fake_data == false)
        {
            data1D_Px->SetLineWidth(2);
            data1D_Px->SetMarkerStyle(20);
            data1D_Px->SetMarkerSize(1.0);
            data1D_Px->SetLineColor(kBlack); // TODO: not needed? I forget reason for adding
            data1D_Px->SetMarkerColor(kBlack); // TODO
            data1D_Px->SetFillColor(kBlack); // TODO
        }
        if(g_mode_fake_data == true)
        {
            fakeData1D_Px->SetLineWidth(2);
            fakeData1D_Px->SetMarkerStyle(20);
            fakeData1D_Px->SetMarkerSize(1.0);
            fakeData1D_Px->SetLineColor(kBlack);
            fakeData1D_Px->SetMarkerColor(kBlack);
            fakeData1D_Px->SetFillColor(kBlack);
        }
        TString Ndata_Px_str;
        if(g_mode_fake_data == false)
        {
            Ndata_Px_str.Form("%i", (int)data1D_Px->Integral()); // TODO: float?
            data1D_Px->SetTitle("Data (" + Ndata_Px_str + ")");
        }
        TString Nfakedata_Px_str;
        if(g_mode_fake_data == true)
        {
            Nfakedata_Px_str.Form("%i", (int)fakeData1D_Px->Integral()); // TODO: float?
            fakeData1D_Px->SetTitle("Fake Data (" + Nfakedata_Px_str + ")");
        }
        if(g_mode_fake_data == false)
        {
            //data1D_Px[i]->Draw("PEsames");
            data1D_Px->Draw("PEsame");
            data1D_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px);
        }
        if(g_mode_fake_data == true)
        {
            fakeData1D_Px->Draw("PEsame");
            fakeData1D_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN_Px, PAD_U_Y_MAX_Px);
        }
        //data1D_Px->GetYaxis()->SetRangeUser(PAD_U_Y_MIN, PAD_U_Y_MAX); // TODO???

        //double chi2;
        int nch_Px = -1;
        if(g_mode_fake_data == false)
        {
            nch_Px = get_nch_1D(hAllMC1D_Px, data1D_Px);
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
            nch_Px = get_nch_1D(hAllMC1D_Px, fakeData1D_Px);
        }
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
            fval_Px = calc_chi2_draw_metric_CHISQ_ERR_SQRT_MC(data1D_Px, hAllMC1D_Px);
        }
        if(g_mode_fake_data == true)
        {
            //fval_Px = calc_chi2_draw(fakeData1D_Px, hAllMC1D_Px);
            fval_Px = calc_chi2_draw_metric_CHISQ_ERR_SQRT_MC(fakeData1D_Px, hAllMC1D_Px);
        }
        TString fval_Px_str;
        fval_Px_str.Form("%.1f", fval_Px);
        int ndf_Px = nch_Px - drawinputdata.nfp;
        TString ndf_Px_str;
        ndf_Px_str.Form("%i", ndf_Px);
        /*
        mychi2_str.Form("%4.3f", mychi2);
        */

        
        TGaxis *axis = new TGaxis(0.0, PAD_U_Y_MIN_Px + 0.01,
                                  0.0, PAD_U_Y_MAX_Px,
                                  PAD_U_Y_MIN_Px + 0.01, PAD_U_Y_MAX_Px,
                                  510, "");
        axis->SetLabelFont(43);
        axis->SetLabelSize(15);
        //axis->SetTickLength(0.05 / (0.7 * 0.9)); // TODO: think the margin is 0.1
        //axis->SetTickSize(0.05 / (0.7 * 0.9));
        axis->Draw();

        TGaxis *axis2 = new TGaxis(5.0, PAD_U_Y_MIN_Px + 0.01,
                                   5.0, PAD_U_Y_MAX_Px,
                                   PAD_U_Y_MIN_Px + 0.01, PAD_U_Y_MAX_Px,
                                   510, "+");
        axis2->SetLabelFont(43);
        axis2->SetLabelSize(0);
        //axis2->SetTickLength(0.05 / (0.3 * 0.6)); // TODO: the margin is 0.4
        //axis2->SetTickSize(0.05 / (0.3 * 0.6));
        axis2->Draw();

        TLegend *leg = new TLegend(0.615, 0.1, 0.855, 0.85);
        if(g_mode_fake_data == false)
        {
            leg->AddEntry(data1D_Px, "Data (" + Ndata_Px_str + ")", "PE"); // TODO PEL ??? works?
        }
        if(g_mode_fake_data == true)
        {
            leg->AddEntry(fakeData1D_Px, "Fake Data (" + Nfakedata_Px_str + ")", "PE");
        }
        leg->AddEntry(hAllMC1D_Px, "Total MC (" + Nmc_Px_str + ")", "L");
        leg->AddEntry(h_2nubb_Px, "2#nu#beta#beta", "F");
        leg->AddEntry(h_tl208_int_Px, "^{208}Tl Int", "F");
        leg->AddEntry(h_bi214_int_Px, "^{214}Bi Int", "F");
        leg->AddEntry(h_bi207_int_Px, "^{207}Bi Int", "F");
        leg->AddEntry(h_internal_Px, "Internal", "F");
        leg->AddEntry(h_neighbour_Px, "Neighbour Foil", "F");
        leg->AddEntry(h_radon_Px, "Radon", "F");
        leg->AddEntry(h_external_Px, "External", "F");
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

        //TLatex latexlabel;
        //latexlabel.SetNDC();
        //latexlabel.SetTextFont(62);
        //latexlabel.SetTextSize(0.035);
        //latexlabel.DrawLatex(0.63, 0.23, "#frac{#chi^{2}}{ndf} = #frac{" + chi2_str + "}{" + ndf_str + "}");
        

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
        xilatex.DrawLatex(0.425, 0.70, xilatexstr);

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
        chilatex.DrawLatex(0.425, 0.45, chilatexstr);


        // second pad
        p1_Px->cd();
        hRatio_Px->SetMarkerStyle(20);
        hRatio_Px->SetMarkerSize(1.0);
        hRatio_Px->SetLineWidth(2);
        hRatio_Px->Draw("AXIS");
        TLine *linezero = new TLine(0.0, 1.0, 5.0, 1.0);
        linezero->SetLineWidth(2);
        linezero->SetLineColor(kGray + 2);
        linezero->Draw();
        hRatio_Px->Draw("EPsame");
        //TLine *zeroline = new TLine(0.0, 0.0, 5.0, 0.0);
        //zeroline->Draw();



        TGaxis *axis3 = new TGaxis(0.0, PAD_L_Y_MIN_Px, 0.0, PAD_L_Y_MAX_Px - 0.0001, PAD_L_Y_MIN_Px, PAD_L_Y_MAX_Px - 0.0001, 503, "");
        axis3->SetLabelFont(43);
        axis3->SetLabelSize(15);
        axis3->Draw();

        TGaxis *axis4 = new TGaxis(5.0, PAD_L_Y_MIN_Px, 5.0, PAD_L_Y_MAX_Px - 0.0001, PAD_L_Y_MIN_Px, PAD_L_Y_MAX_Px - 0.0001, 503, "+");
        axis4->SetLabelFont(43);
        axis4->SetLabelSize(0);
        axis4->Draw();

        



        canvas_Px->Show();

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
                    canvas_saveas_helper(dir, name, canvas_Px);
                }
            }
        }


        draw_output_data drawoutputdata_Px;
        drawoutputdata_Px.stacks1D = stacks1D_Px;
        drawoutputdata_Px.h_2nubb = h_2nubb_Px;
        drawoutputdata_Px.h_tl208_int = h_tl208_int_Px;
        drawoutputdata_Px.h_bi214_int = h_bi214_int_Px;
        drawoutputdata_Px.h_bi207_int = h_bi207_int_Px;
        drawoutputdata_Px.h_internal = h_internal_Px;
        drawoutputdata_Px.h_external = h_external_Px;
        drawoutputdata_Px.h_radon = h_radon_Px;
        drawoutputdata_Px.h_neighbour = h_neighbour_Px;
        drawoutputdata_Px.h_other = h_other_Px;
        drawoutputdata_Px.canvas = canvas_Px;
        drawoutputdata_Px.hRatio = hRatio_Px;
        drawoutputdata_Px.hAllMC1D = hAllMC1D_Px;
        drawoutputdata_Px.data1D = data1D_Px;
        drawoutputdata_Px.fakeData1D = fakeData1D_Px;
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

#endif //NEWLOGLIKFITTER_DRAWCHANNEL_H
