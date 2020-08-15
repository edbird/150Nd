#ifndef NEWLOGLIKFITTER_DRAW_2D_H
#define NEWLOGLIKFITTER_DRAW_2D_H

// TODO: draw functions for other fit histograms, particularly hHighLowEnergy
void draw_2D(
    const std::vector<double> &params,
    const std::vector<double> &param_errs,
    const std::string &saveas_filename,
    TH1D* const hHighEnergy_allMC,
    TH1D* const hLowEnergy_allMC,
    TH1D* const hHighEnergy_data,
    TH1D* const hLowEnergy_data)
{
    
#if 0
    TH2D *h_data[number2DHists] = { nullptr };
    TH2D *h_MC[number2DHists] = { nullptr };
    TH2D *h_diff[number2DHists] = { nullptr };
    TH2D *h_ratio[number2DHists] = { nullptr };
    TH2D *h_err[number2DHists] = { nullptr };
    TH2D *h_pull[number2DHists] = { nullptr };

    TCanvas *c_data[number2DHists] = { nullptr };
    TCanvas *c_MC[number2DHists] = { nullptr };
    TCanvas *c_diff[number2DHists] = { nullptr };
    TCanvas *c_ratio[number2DHists] = { nullptr };
    TCanvas *c_err[number2DHists] = { nullptr };
    TCanvas *c_pull[number2DHists] = { nullptr };

    TH1D *h_data_low_e_projection[number2DHists] = { nullptr };
    TH1D *h_data_high_e_projection[number2DHists] = { nullptr };
    TH1D *h_MC_low_e_projection[number2DHists] = { nullptr };
    TH1D *h_MC_high_e_projection[number2DHists] = { nullptr };

    TCanvas *c_low_e_projection[number2DHists] = { nullptr };
    TCanvas *c_high_e_projection[number2DHists] = { nullptr };
    //TCanvas c_MC_low_e_projection[number2DHists] = { nullptr };
    //TCanvas c_MC_high_e_projection[number2DHists] = { nullptr };

    TH2D *hAllMC2D[number2DHists] = { nullptr };
    TH2D *data2D[number2DHists] = { nullptr };



    // TODO NOTE: added some names in clone call functions, might need to
    // remove these as they may be breaking code? check against 1D draw
    // function

    // each channel 1D hists
    // this is for(i = 0; i < 1; ++ i)
    for(int channel = 0; channel < allDataSamples2D->GetEntries(); ++ channel)
    {

        TString i_str;
        i_str.Form("%i", channel);

        data2D[channel] = (TH2D*)allDataSamples2D->At(channel)->Clone("data2D_" + i_str);
        h_ratio[channel] = (TH2D*)allDataSamples2D->At(channel)->Clone("h_ratio" + i_str);
        h_diff[channel] = (TH2D*)allDataSamples2D->At(channel)->Clone("h_diff" + i_str);
    
        c_data[channel] = new TCanvas("c_data_2D" + i_str, "c_data_2D" + i_str);
        c_data[channel]->SetFillColor(kWhite);

        c_MC[channel] = new TCanvas("c_MC_2D" + i_str, "c_MC_2D" + i_str);
        c_MC[channel]->SetFillColor(kWhite);
        
        c_diff[channel] = new TCanvas("c_diff" + i_str, "c_diff" + i_str);
        c_diff[channel]->SetFillColor(kWhite);

        c_ratio[channel] = new TCanvas("c_ratio_2D" + i_str, "c_ratio_2D" + i_str);
        c_ratio[channel]->SetFillColor(kWhite);

        c_err[channel] = new TCanvas("c_err_2D" + i_str, "c_err_2D" + i_str);
        c_err[channel]->SetFillColor(kWhite);

        c_pull[channel] = new TCanvas("c_pull_2D" + i_str, "c_pull_2D" + i_str);
        c_pull[channel]->SetFillColor(kWhite);

        c_low_e_projection[channel] = new TCanvas("c_low_e_projection" + i_str, "c_low_e_projection" + i_str);
        c_low_e_projection[channel]->SetFillColor(kWhite);

        c_high_e_projection[channel] = new TCanvas("c_high_e_projection" + i_str, "c_high_e_projection" + i_str);
        c_high_e_projection[channel]->SetFillColor(kWhite);


        // add each MC sample
        for(int j = 0; j < allMCSamples2D[channel]->GetEntries(); ++ j)
        {

            TString j_str;
            j_str.Form("%i", j);

            TH2D *tmpHist_draw2D = (TH2D*)allMCSamples2D[channel]->At(j)->Clone();
            TString tmpName = tmpHist_draw2D->GetName();

            //std::cout << "looking for " << tmpName << std::endl;
            int which_param = -1;
            bool found_param = false;

            // used later
            double activity_scale_branching_ratio = 1.0;

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
                        //std::cout << "tmp_sample_name=" << tmp_sample_name << std::endl;
                        //std::cin.get();
                        activity_scale_branching_ratio = 0.36;
                    }
                    */

                    //std::cout << "tmp_sample_name=" << tmp_sample_name << std::endl;
                    if(MCNameToParamNumberMap.count(tmp_sample_name) > 0)
                    {
                        int paramNumber = MCNameToParamNumberMap.at(tmp_sample_name);
                        // TODO: removed std::string, change tmpName type to be std::string from TString
                    
                        //std::cout << "paramNumber=" << paramNumber << " -> tmp_sample_name=" << tmp_sample_name << " ~> tmpName=" << tmpName << std::endl;                    
                        //which_param = paramNumber;
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
                // scale histogram to correct size using output parameter
                // from fit
                if(which_param >= numberEnabledParams)
                {
                    std::cout << "throwing exception, which_param=" << which_param << std::endl;
                    throw std::runtime_error("which_param invalid value");
                }


                // no error thrown, which_param is presumably the correct index
                Double_t activity_scale = params[which_param]; // * activity_scale_branching_ratio;
                tmpHist_draw2D->Scale(activity_scale);
                // TODO: fix this

                
                if(tmpHist_draw2D->Integral() > 0)
                {
                    //h_data_div_MC[i]->Add((TH2D*)tmpHist_draw2D->Clone(), -1.0);
                    
                    if(j == 0)
                    {
                        hAllMC2D[channel] = (TH2D*)tmpHist_draw2D->Clone("Total_MC_" + i_str);
                    }
                    else
                    {
                        hAllMC2D[channel]->Add((TH2D*)tmpHist_draw2D);
                    }
	            }
                else
                {
                    //std::cout << "not adding to stack, Integral() <= 0: " << tmpHist_draw2D->GetName() << std::endl;
                }
            }
            else
            {
                std::cout << "error could not find histogram: tmpName=" << tmpName << std::endl;
            } 

        }


        h_data_low_e_projection[channel] = data2D[channel]->ProjectionX("_px", 0, -1, "e");
        h_data_high_e_projection[channel] = data2D[channel]->ProjectionY("_py", 0, -1, "e");
        h_MC_low_e_projection[channel] = hAllMC2D[channel]->ProjectionX("_px", 0, -1, "e");
        h_MC_high_e_projection[channel] = hAllMC2D[channel]->ProjectionY("_py", 0, -1, "e");

        c_low_e_projection[channel]->cd();
        h_MC_low_e_projection[channel]->SetTitle("");
        h_MC_low_e_projection[channel]->SetFillColor(kBlue);
        h_MC_low_e_projection[channel]->SetFillStyle(1001);
        h_MC_low_e_projection[channel]->SetLineColor(kBlue);
        h_MC_low_e_projection[channel]->SetLineWidth(0);
        h_MC_low_e_projection[channel]->Divide(hLowEnergy_allMC);
        h_MC_low_e_projection[channel]->Draw("hist");
        if(hLowEnergy_allMC == nullptr)
        {
            std::cout << "nullptr!" << std::endl;
            std::cin.get();
        }
        else
        {
            hLowEnergy_allMC->SetFillColor(kRed);
            hLowEnergy_allMC->SetFillStyle(1001);
            hLowEnergy_allMC->SetLineColor(kRed);
            hLowEnergy_allMC->SetLineWidth(0);
//            hLowEnergy_allMC->Draw("histsame");
        }
        h_data_low_e_projection[channel]->SetMarkerStyle(20);
        h_data_low_e_projection[channel]->Divide(hLowEnergy_data);
        h_data_low_e_projection[channel]->Draw("PEsame");

        c_high_e_projection[channel]->cd();
        h_MC_high_e_projection[channel]->SetTitle("");
        h_MC_high_e_projection[channel]->SetFillColor(kBlue);
        h_MC_high_e_projection[channel]->SetFillStyle(1001);
        h_MC_high_e_projection[channel]->SetLineColor(kBlue);
        h_MC_high_e_projection[channel]->SetLineWidth(0);
        h_MC_high_e_projection[channel]->Divide(hHighEnergy_allMC);
        h_MC_high_e_projection[channel]->Draw("hist");
        if(hHighEnergy_allMC == nullptr)
        {
            std::cout << "nullptr!" << std::endl;
            std::cin.get();
        }
        else
        {
            hHighEnergy_allMC->SetFillColor(kRed);
            hHighEnergy_allMC->SetFillStyle(1001);
            hHighEnergy_allMC->SetLineColor(kRed);
            hHighEnergy_allMC->SetLineWidth(0);
//            hHighEnergy_allMC->Draw("histsame");
            std::cout << "----->>>>> Integral()=" << hHighEnergy_allMC->Integral() << std::endl;
        }
        h_data_high_e_projection[channel]->SetMarkerStyle(20);
        h_data_high_e_projection[channel]->Divide(hHighEnergy_data);
        h_data_high_e_projection[channel]->Draw("PEsame");





        h_data[channel] = (TH2D*)data2D[channel]->Clone("h_data_2D_" + i_str);
        h_MC[channel] = (TH2D*)hAllMC2D[channel]->Clone("h_MC_2D_" + i_str);

        h_err[channel] = (TH2D*)h_data[channel]->Clone("h_err_2D_" + i_str);
        //h_err[channel] = new TH2D("h_err_2D_" + i_str, "h_err_2D_" + i_str,
        //                          50, 0.0, 5.0, 50, 0.0, 5.0);
        for(int i = 1; i <= h_data[channel]->GetNbinsX(); ++ i)
        {
            for(int j = 1; j <= h_data[channel]->GetNbinsY(); ++ j)
            {

                //Double_t content = h_data[channel]->GetBinContent(i, j);
                Double_t error = h_data[channel]->GetBinError(i, j);

                h_err[channel]->SetBinContent(i, j, error);
                //h_err[channel]->SetBinContent(i, j, content);
                //Double_t ix = h_err[channel]->GetXaxis()->GetBinCenter(i);
                //Double_t jy = h_err[channel]->GetYaxis()->GetBinCenter(j);
                //h_err[channel]->Fill(ix, jy, content);
                //h_err[channel]->Fill(ix, jy, error);
                //h_err[channel]->SetBinError(0.0, i, j);
            }
        }

        ///////////////////////////////////////////////////////////////////////

        c_data[channel]->cd();
        h_data[channel]->SetContour(1000);
        h_data[channel]->SetTitle("");
        h_data[channel]->Draw("colz");

        /*
        {
            std::string filename = std::string(c_data[channel]->GetName()) + std::string(".png");
            c_data[channel]->SaveAs(filename.c_str());
        }
        */

        ///////////////////////////////////////////////////////////////////////

        c_MC[channel]->cd();
        h_MC[channel]->SetTitle("");
        h_MC[channel]->SetContour(1000);
        h_MC[channel]->Draw("colz");

        /*
        {
            std::string filename = std::string(c_MC[channel]->GetName()) + std::string(".png");
            c_MC[channel]->SaveAs(filename.c_str());
        }
        */

        ///////////////////////////////////////////////////////////////////////


        h_diff[channel]->Add(hAllMC2D[channel], -1.0);

        c_diff[channel]->cd();
        h_diff[channel]->SetContour(1000);
        h_diff[channel]->SetTitle("");
        h_diff[channel]->Draw("colz");

        // change to division
        //h_data_div_MC[channel] = (TH2D*)data2D[channel]->Clone("h_data_div_MC_2D");
        h_ratio[channel]->Divide(hAllMC2D[channel]);

        h_pull[channel] = (TH2D*)h_diff[channel]->Clone("h_pull_2D_" + i_str);
        h_pull[channel]->Divide(h_err[channel]);

        c_ratio[channel]->cd();
        h_ratio[channel]->SetContour(1000);
        h_ratio[channel]->SetTitle("");
        h_ratio[channel]->Draw("colz");
        //hAllMC1D[i]->Draw("hist sames");
        //TString Nmc_str;
        //Nmc_str.Form("%i", (int)h_data_div_MC[i]->Integral()); // TODO: float?
        //hAllMC1D[i]->SetTitle("Total MC (" + Nmc_str + ")");
        //hAllMC1D[i]->Draw("hist same");
        //data1D[i]->SetLineWidth(2);
        //data1D[i]->SetMarkerStyle(20);
        //data1D[i]->SetMarkerSize(1.0);
        //TString Ndata_str;
        //Ndata_str.Form("%i", (int)data1D[i]->Integral()); // TODO: float?
        //data1D[i]->SetTitle("Data (" + Ndata_str + ")");
        //data1D[i]->Draw("PEsames");
        //data1D[i]->Draw("PEsame");

        double chi2;
        int ndf;
        int igood;
        TString chi2_str;
        TString ndf_str;

        // TODO: should chisquare value include the constraints? because at
        // the moment it does not

        // TODO: chi2 value is different from fit_2e code
        double prob = data2D[channel]->Chi2TestX(hAllMC2D[channel], chi2, ndf, igood, "UW");
        // TODO: check if I can get fcn value from the minuit fit object
        chi2_str.Form("%4.3f", chi2);
        ndf_str.Form("%i", ndf);


        TLegend *leg = new TLegend(0.6, 0.88, 0.88, 0.7);
        //leg->AddEntry(data1D[i], "Data (" + Ndata_str + ")", "PE");
        //leg->AddEntry(hAllMC1D[i], "Total MC (" + Nmc_str + ")", "L");
        leg->AddEntry(h_ratio[channel], "Data / MC", "F");
        leg->AddEntry((TObject*)nullptr, "#chi^{2}/ndf=" + chi2_str + "/" + ndf_str, "");
        //leg->AddEntry(h_other[i], "other", "f");
        leg->SetBorderSize(5);
        leg->SetFillColor(0);
        leg->SetTextFont(62);
        leg->SetTextSize(0.035);
        leg->SetShadowColor(kBlack);
        leg->Draw("BR");


        TLatex latexlabel;
        latexlabel.SetNDC();
        latexlabel.SetTextFont(62);
        latexlabel.SetTextSize(0.035);
        latexlabel.DrawLatex(0.63, 0.23, "#frac{#chi^{2}}{ndf} = #frac{" + chi2_str + "}{" + ndf_str + "}");



        c_err[channel]->cd();
        h_err[channel]->SetContour(1000);
        h_err[channel]->SetTitle("");
        h_err[channel]->Draw("colz");

        c_pull[channel]->cd();
        h_pull[channel]->SetContour(1000);
        h_pull[channel]->SetTitle("");
        h_pull[channel]->Draw("colz");
        
        Double_t chi2_2 = 0.0;
        for(int i = 1; i <= h_pull[channel]->GetNbinsX(); ++ i)
        {
            for(int j = 1; j <= h_pull[channel]->GetNbinsY(); ++ j)
            {
                chi2_2 += std::pow(h_pull[channel]->GetBinContent(i, j), 2.0);
            }
        }
        TString chi2_str2;
        chi2_str2.Form("%.3f", chi2_2);
        TString ndf_str2 = ndf_str;
        TLatex llpull0;
        llpull0.SetNDC();
        llpull0.SetTextFont(43);
        llpull0.SetTextSize(25);
        llpull0.DrawLatex(0.6, 0.2, "#frac{#chi^{2}}{ndf} = #frac{" + chi2_str2 + "}{" + ndf_str2 + "}");
    

        std::string base_name;
        std::string extension;
        filename_split_extension(saveas_filename, base_name, extension);

        {
            std::string dir = base_name + "_c_MC";
            std::string name = base_name + "_c_MC" + "_" + std::string(i_str);
            std::string name_plus_ext = name + extension;
            canvas_saveas_helper(dir, name_plus_ext, c_MC[channel]);
        }

        {
            std::string dir = base_name + "_c_data";
            std::string name = base_name + "_c_data" + "_" + std::string(i_str);
            std::string name_plus_ext = name + extension;
            canvas_saveas_helper(dir, name_plus_ext, c_data[channel]);
        }

        {
            std::string dir = base_name + "_c_ratio";
            std::string name = base_name + "_c_ratio" + "_" + std::string(i_str);
            std::string name_plus_ext = name + extension;
            canvas_saveas_helper(dir, name_plus_ext, c_ratio[channel]);
        }

        {
            std::string dir = base_name + "_c_err";
            std::string name = base_name + "_c_err" + "_" + std::string(i_str);
            std::string name_plus_ext = name + extension;
            canvas_saveas_helper(dir, name_plus_ext, c_err[channel]);
        }

        {
            std::string dir = base_name + "_c_pull";
            std::string name = base_name + "_c_pull" + "_" + std::string(i_str);
            std::string name_plus_ext = name + extension;
            canvas_saveas_helper(dir, name_plus_ext, c_pull[channel]);
        }


    }

#endif
}


#endif //NEWLOGLIKFITTER_DRAW_2D_H
