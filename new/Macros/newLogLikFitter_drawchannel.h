#ifndef NEWLOGLIKFITTER_DRAWCHANNEL_H
#define NEWLOGLIKFITTER_DRAWCHANNEL_H


void draw_channel(const int channel,
                  const Double_t *const p,
                  const double fval,
                  const std::string& saveas_filename,
                  const std::string& saveas_dir = ".")
{

    THStack *stacks1D;
    TH1F *h_2nubb = nullptr;
    TH1F *h_tl208_int = nullptr;
    TH1F *h_bi214_int = nullptr;
    TH1F *h_bi207_int = nullptr;
    TH1F *h_internal = nullptr;
    TH1F *h_external = nullptr;
    TH1F *h_radon = nullptr; 
    TH1F *h_neighbours = nullptr;
    TH1F *h_other = nullptr;


    TCanvas *c;
    TPad *p0;
    TPad *p1;
    TH1F *hRatio;
    TH1F *hAllMC1D;
    TH1F *data1D;

    std::cout << "debug: number of data samples: " << allDataSamples1D->GetEntries() << std::endl;
    std::cout << "debug: number of MC samples: " << allMCSamples1D[0]->GetEntries() << std::endl;


    // additional array index
    data1D = (TH1F*)allDataSamples1D->At(channel)->Clone();

    TString channel_str;
    channel_str.Form("%i", channel);

    c = new TCanvas("c" + channel_str, "c" + channel_str);
    c->SetFillColor(kWhite);

    
    stacks1D = new THStack("stacks1D" + channel_str, channel_str);


    TH1F *tmpHist_draw1D;


    // TODO i should be channel here?
    for(int j = 0; j < allMCSamples1D[channel]->GetEntries(); j++)
    {

        TString j_str;
        j_str.Form("%i", j);

        tmpHist_draw1D = (TH1F*)allMCSamples1D[channel]->At(j)->Clone();
        TString tmpName = tmpHist_draw1D->GetName();

        //std::cout << "looking for " << tmpName << std::endl;
        int which_param = -1;
        bool found_param = false;

        // used later
        //double activity_scale_branching_ratio = 1.0;

        // get index for parameter
        {
            std::string tmp_hist_name(tmpName);
            auto i_start = tmp_hist_name.find('_') + 1;
            auto i_end = tmp_hist_name.rfind('_');
            if(i_end - i_start > 0)
            {
                std::string tmp_sample_name = tmp_hist_name.substr(i_start, i_end - i_start);

                // set branching ratio fraction
                /*
                if(tmp_sample_name == std::string("tl208_int_rot") ||
                   tmp_sample_name == std::string("tl208_feShield") ||
                   tmp_sample_name == std::string("tl208_pmt"))
                {
                    activity_scale_branching_ratio = 0.36;
                }
                */

                if(MCNameToParamNumberMap.count(tmp_sample_name) > 0)
                {
                    int paramNumber = MCNameToParamNumberMap.at(tmp_sample_name);
                    // TODO: removed std::string, change tmpName type to be std::string from TString
                
                    //std::cout << "paramNumber=" << paramNumber << " -> tmp_sample_name=" << tmp_sample_name << " ~> tmpName=" << tmpName << std::endl;                    
                    which_param = paramNumberToMinuitParamNumberMap.at(paramNumber);
                    found_param = true;
                }
                else
                {
                   std::cout << "ERROR: could not find " << tmp_sample_name << " in MCNameToParamNumberMap" << std::endl;
                }
            }
        }

        if(found_param == true)
        {
            //std::cout << "found histogram: tmpName=" << tmpName << " which_param=" << which_param << std::endl;

            // scale histogram to correct size using output parameter
            // from fit
            if(which_param >= numberEnabledParams)
            {
                std::cout << "throwing exception, which_param=" << which_param << std::endl;
                throw std::runtime_error("which_param invalid value");
            }

            // no error thrown, which_param is presumably the correct index
            //Double_t activity_scale = AdjustActs[which_param] * activity_scale_branching_ratio;
            Double_t activity_scale = p[which_param]; // * activity_scale_branching_ratio; // TODO: caution: changed from AdjustActs to p
            tmpHist_draw1D->Scale(activity_scale);

            if(tmpHist_draw1D->Integral() > 0)
            {

                stacks1D->Add((TH1F*)tmpHist_draw1D->Clone());
                stacker_helper(tmpHist_draw1D,
                               h_2nubb,
                               h_tl208_int,
                               h_bi214_int,
                               h_bi207_int,
                               h_internal,
                               h_neighbours,
                               h_radon,
                               h_external,
                               h_other);


                if(j == 0)
                {
                    hAllMC1D = (TH1F*)tmpHist_draw1D->Clone("Total MC");
                    //hAllMC1D[i] = (TH1F*)tmpHist_drawpointer->Clone("Total MC");
                }
                else
                {
                    hAllMC1D->Add((TH1F*)tmpHist_draw1D);
                    //hAllMC1D[i]->Add((TH1F*)tmpHist_drawpointer);
                }
            }
            else
            {
                //std::cout << "not adding to stack, Integral() <= 0: " << tmpHist_draw1D->GetName() << std::endl;
            }
        }
        else
        {
            std::cout << "error could not find histogram: tmpName=" << tmpName << std::endl;
        } 

    }


    //stacks1D[i]->SetMaximum(350.);
    //stacks1D[i]->Draw("hist");

    //stacks1D_2nubb[i]->SetMaximum(350.0);
    //stacks1D_2nubb[i]->Draw("hist");
    //stacks1D_tl208_int[i]->Draw("histsame");
    //stacks1D_bi214_int[i]->Draw("histsame");
    //stacks1D_bi207_int[i]->Draw("histsame");
    //stacks1D_internal[i]->Draw("histsame");
    //stacks1D_external[i]->Draw("histsame");
    //stacks1D_radon[i]->Draw("histsame");
    //stacks1D_neighbours[i]->Draw("histsame");
    //stacks1D_other[i]->Draw("histsame");


    THStack *stacks1D_major;
    stacks1D_major = new THStack("stacks1D_major" + channel_str, channel_str);
    stacker_helper_2(stacks1D_major,
                     h_2nubb,
                     h_tl208_int,
                     h_bi214_int,
                     h_bi207_int,
                     h_internal,
                     h_neighbours,
                     h_radon,
                     h_external,
                     h_other);


    double PAD_U_Y_MIN = 0.0;
    double PAD_U_Y_MAX = 500.0;
    double PAD_L_Y_MAX = 3.0;
    double PAD_L_Y_MIN = 0.0;

    if(channel == 0)
    {
        PAD_U_Y_MAX = 350.0;
    }
    else if(channel == 1)
    {
        PAD_U_Y_MAX = 1000.0;
        PAD_U_Y_MAX = 1200.0;
    }
    else if(channel == 2)
    {
        PAD_U_Y_MAX = 450.0;
    }
    else if(channel == 3)
    {
        PAD_U_Y_MAX = 1000.0;
    }
    else if(channel == 4)
    {
        PAD_U_Y_MAX = 400.0;
    }
    else if(channel == 5)
    {
        PAD_U_Y_MAX = 600.0;
    }
    else
    {
        PAD_U_Y_MAX = 350.0;
    }
    //stacks1D_major[i]->Draw("hist");
    //stacks1D_major[i]->SetMaximum(PAD_U_Y_MAX);
    //stacks1D_major[i]->SetMinimum(PAD_U_Y_MIN);
//    stacks1D_major->GetYaxis()->SetRangeUser(PAD_U_Y_MIN, PAD_U_Y_MAX);
    // TODO moved
//    hAllMC1D->SetMaximum(PAD_U_Y_MAX);
//    hAllMC1D->SetMinimum(PAD_U_Y_MIN);
    hAllMC1D->GetYaxis()->SetRangeUser(PAD_U_Y_MIN, PAD_U_Y_MAX);
//    data1D->SetMaximum(PAD_U_Y_MAX);
//    data1D->SetMinimum(PAD_U_Y_MIN);
    data1D->GetYaxis()->SetRangeUser(PAD_U_Y_MIN, PAD_U_Y_MAX);


    hRatio = (TH1F*)data1D->Clone();
    hRatio->Sumw2();
    hRatio->Divide(hAllMC1D);
    hRatio->SetTitle("");


    if(channel == 0)
    {
        hRatio->GetXaxis()->SetTitle("2e Electron Energy [MeV]");
    }
    else if(channel == 1)
    {
        hRatio->GetXaxis()->SetTitle("Single Electron Energy [MeV]");
    }
    else if(channel == 2)
    {
        hRatio->GetXaxis()->SetTitle("High Energy Electron [MeV]");
    }
    else if(channel == 3)
    {
        hRatio->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
    }
    else if(channel == 4)
    {
        hRatio->GetXaxis()->SetTitle("Energy Sum [MeV]");
    }
    else if(channel == 5)
    {
        hRatio->GetXaxis()->SetTitle("Energy Diff [MeV]");
    }

    //hRatio[i]->SetMaximum(PAD_L_Y_MAX);
    //hRatio[i]->SetMinimum(PAD_L_Y_MIN);
    hRatio->GetYaxis()->SetRangeUser(PAD_L_Y_MIN, PAD_L_Y_MAX);

    hRatio->GetXaxis()->SetTickSize(0.1);

    //hRatio->GetXaxis()->SetTitle("Electron Energy [MeV]");
    hRatio->GetXaxis()->SetTitleFont(43);
    hRatio->GetXaxis()->SetTitleSize(20);
    hRatio->GetXaxis()->SetTitleOffset(3.0);

    hRatio->GetYaxis()->SetTitle("data / MC");
    hRatio->GetYaxis()->SetTitleFont(43);
    hRatio->GetYaxis()->SetTitleSize(20);
    hRatio->GetYaxis()->SetTitleOffset(1.0);

    hRatio->GetYaxis()->SetLabelFont(43);
    hRatio->GetYaxis()->SetLabelSize(0.0 * 31);

    hRatio->GetYaxis()->SetTickSize(0.0);

    hRatio->GetXaxis()->SetLabelFont(43);
    hRatio->GetXaxis()->SetLabelSize(15);

    // in code copying from, canvas alloc. here
    c->cd();
    c->SetBottomMargin(0.0);
    p0 = new TPad("pad0", "pad0", 0.0, 0.3, 1.0, 1.0);
    p0->SetBottomMargin(0.0);
    //p0->SetGridx(1);
    //p0->SetGridy(1);
    p0->SetGrid(0, 0);
    p0->SetTicks(2, 2);
    p0->Draw();

    c->cd();
    p1 = new TPad("pad1", "pad1", 0.0, 0.0, 1.0, 0.3);
    p1->SetTopMargin(0.0);
    p1->SetBottomMargin(0.4);
    //p1->SetGridx(1);
    //p1->SetGridy(1);
    p1->SetGrid(1, 1);
    p1->SetTicks(2, 2);
    p1->Draw();

    p0->cd();

    // draw regular pad1
    // leave alone for now
    // TODO: axis
    // copy code from other file
    hAllMC1D->SetTitle("");

    stacks1D_major->Draw("hist"); // this didn't used to be here
    stacks1D_major->SetMaximum(PAD_U_Y_MAX);
    stacks1D_major->SetMinimum(PAD_U_Y_MIN);
    //stacks1D_major->Draw("hist"); // this used to be uncommented
    stacks1D_major->GetYaxis()->SetRangeUser(PAD_U_Y_MIN, PAD_U_Y_MAX); // this used to be commented
    //stacks1D_major->Draw("hist");
    //p0->Update();
//    stacks1D_major->SetMaximum(PAD_U_Y_MAX);
//    stacks1D_major->SetMinimum(PAD_U_Y_MIN);
//    stacks1D_major->GetYaxis()->SetRangeUser(PAD_U_Y_MIN, PAD_U_Y_MAX);

    stacks1D_major->SetTitle("");

    stacks1D_major->GetYaxis()->SetTickSize(0.0);
    //stacks1D_major->GetXaxis()->SetTickSize(0.0);

    stacks1D_major->GetYaxis()->SetTitle("Spectrum");
    stacks1D_major->GetYaxis()->SetTitleSize(20);
    stacks1D_major->GetYaxis()->SetTitleFont(43);
    stacks1D_major->GetYaxis()->SetTitleOffset(1.0);

    stacks1D_major->GetYaxis()->SetLabelFont(43);
    stacks1D_major->GetYaxis()->SetLabelSize(0.0 * 31);

    stacks1D_major->GetXaxis()->SetTitle("");
    stacks1D_major->GetXaxis()->SetTitleSize(0);
    stacks1D_major->GetXaxis()->SetTitleFont(43);
    stacks1D_major->GetXaxis()->SetTitleOffset(1.0);

    stacks1D_major->GetXaxis()->SetLabelSize(0);
    stacks1D_major->GetXaxis()->SetLabelFont(63);


    hAllMC1D->SetLineWidth(2);
    hAllMC1D->SetLineColor(kBlack);
    hAllMC1D->SetFillColor(kWhite);
    hAllMC1D->SetFillStyle(0);
    //hAllMC1D->Sumw2();
    //hAllMC1D->Draw("hist sames");
    TString Nmc_str;
    Nmc_str.Form("%i", (int)hAllMC1D->Integral()); // TODO: float?
    hAllMC1D->SetTitle("Total MC (" + Nmc_str + ")");
    hAllMC1D->Draw("hist same");
    hAllMC1D->GetYaxis()->SetRangeUser(PAD_U_Y_MIN, PAD_U_Y_MAX);
    data1D->SetLineWidth(2);
    data1D->SetMarkerStyle(20);
    data1D->SetMarkerSize(1.0);
    TString Ndata_str;
    Ndata_str.Form("%i", (int)data1D->Integral()); // TODO: float?
    data1D->SetTitle("Data (" + Ndata_str + ")");
    //data1D->Draw("PEsames");
    data1D->Draw("PEsame");
    data1D->GetYaxis()->SetRangeUser(PAD_U_Y_MIN, PAD_U_Y_MAX);

    double chi2;
    int ndf;
    int igood;
    TString chi2_str;
    TString ndf_str;
    TString fval_str;
    TString mychi2_str;

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
    ndf_str.Form("%i", ndf);
    fval_str.Form("%4.3f", fval);
    /*
    mychi2_str.Form("%4.3f", mychi2);
    */

    
    TGaxis *axis = new TGaxis(0.0, PAD_U_Y_MIN + 0.01, 0.0, PAD_U_Y_MAX, PAD_U_Y_MIN + 0.01, PAD_U_Y_MAX, 510, "");
    axis->SetLabelFont(43);
    axis->SetLabelSize(15);
    axis->Draw();

    TGaxis *axis2 = new TGaxis(5.0, PAD_U_Y_MIN + 0.01, 5.0, PAD_U_Y_MAX, PAD_U_Y_MIN + 0.01, PAD_U_Y_MAX, 510, "+");
    axis2->SetLabelFont(43);
    axis2->SetLabelSize(0);
    axis2->Draw();

    TLegend *leg = new TLegend(0.6, 0.1, 0.85, 0.85);
    leg->AddEntry(data1D, "Data (" + Ndata_str + ")", "PE");
    leg->AddEntry(hAllMC1D, "Total MC (" + Nmc_str + ")", "L");
    leg->AddEntry(h_2nubb, "2#nu#beta#beta", "F");
    leg->AddEntry(h_tl208_int, "^{208}Tl Int", "F");
    leg->AddEntry(h_bi214_int, "^{214}Bi Int", "F");
    leg->AddEntry(h_bi207_int, "^{207}Bi Int", "F");
    leg->AddEntry(h_internal, "Internal", "F");
    leg->AddEntry(h_neighbours, "Neighbour Foil", "F");
    leg->AddEntry(h_radon, "Radon", "F");
    leg->AddEntry(h_external, "External", "F");
    //leg->AddEntry((TObject*)nullptr, "#chi^{2}/ndf=" + chi2_str + "/" + ndf_str, "");
    //leg->AddEntry((TObject*)nullptr, "fval #chi^{2}/ndf=" + fval_str + "/" + ndf_str, "");
    leg->AddEntry((TObject*)nullptr, "#chi^{2}/ndf=" + fval_str + "/" + ndf_str, "");
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


    // second pad
    p1->cd();
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(1.0);
    hRatio->Draw("EP");



    TGaxis *axis3 = new TGaxis(0.0, PAD_L_Y_MIN, 0.0, PAD_L_Y_MAX - 0.0001, PAD_L_Y_MIN, PAD_L_Y_MAX - 0.0001, 503, "");
    axis3->SetLabelFont(43);
    axis3->SetLabelSize(15);
    axis3->Draw();

    TGaxis *axis4 = new TGaxis(5.0, PAD_L_Y_MIN, 5.0, PAD_L_Y_MAX - 0.0001, PAD_L_Y_MIN, PAD_L_Y_MAX - 0.0001, 503, "+");
    axis4->SetLabelFont(43);
    axis4->SetLabelSize(0);
    axis4->Draw();




    c->Show();

    if(saveas_filename != "NOSAVE")
    {
        std::string base_name;
        std::string extension;
        filename_split_extension(saveas_filename, base_name, extension);

        //std::string dir = base_name + "_c" + "_" + std::string(channel_str);
        std::string dir = saveas_dir;
        std::string name = base_name + "_c" + "_" + std::string(channel_str) + extension;
        canvas_saveas_helper(dir, name, c);
    }



}

#endif //NEWLOGLIKFITTER_DRAWCHANNEL_H
