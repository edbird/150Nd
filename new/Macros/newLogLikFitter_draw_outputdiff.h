#ifndef NEWLOGLIKFITTER_DRAW_OUTPUTDIFF_H
#define NEWLOGLIKFITTER_DRAW_OUTPUTDIFF_H


// NOTE: obsolete function

// h2 must be orig
void ratioplotOld(TH1D* h1, TH1D* h2, const double xi_31, const std::string &saveas_filename, int flag)
{


    // Define the ratio plot
    TH1D *h3 = (TH1D*)h1->Clone("hratio");
    h3->Sumw2();
    h3->Divide(h2);

   // Define two gaussian histograms. Note the X and Y title are defined
   // at booking time using the convention "Hist_title ; X_title ; Y_title"
//   TH1D *h1 = new TH1D("h1", "Two gaussian plots and their ratio;x title; h1 and h2 gaussian histograms", 100, -5, 5);
//   TH1D *h2 = new TH1D("h2", "h2", 100, -5, 5);
//   h1->FillRandom("gaus");
//   h2->FillRandom("gaus");
    h1->SetTitle("");
    //h1->SetMinimum(0.0);
    //h1->SetMaximum(0.16);
    h1->SetStats(0);          // No statistics on upper plot

    double PAD_U_Y_MIN = 0.0; // -0.1; // 0.0
    double PAD_U_Y_MAX = 0.16; //0.3; // 0.16
    //const double PAD_U_Y_MIN = -0.2; // -0.1; // 0.0
    //const double PAD_U_Y_MAX = 0.4; //0.3; // 0.16
    //double PAD_U_Y_MIN = -5000; // -0.1; // 0.0
    //double PAD_U_Y_MAX = 5000; //0.3; // 0.16
    //const double PAD_U_Y_MIN = 0.0;
    //const double PAD_U_Y_MAX = 0.16;

    double PAD_L_Y_MIN = 0.5; //-1.0; // 0.5
    double PAD_L_Y_MAX = 1.5; //3.0; // 1.5
    //const double PAD_L_Y_MIN = -1.0; // 0.5
    //const double PAD_L_Y_MAX = 3.0; // 1.5
    //double PAD_L_Y_MIN = -1000.0; // 0.5
    //double PAD_L_Y_MAX = 100.0; // 1.5
    //const double PAD_L_Y_MIN = 0.5;
    //const double PAD_L_Y_MAX = 1.5;

    if(flag == 0)
    {
        PAD_U_Y_MIN = -3000.0;
        PAD_U_Y_MAX = 1000.0;
        PAD_L_Y_MIN = -10.0;
        PAD_L_Y_MAX = 10.0;
    }
    else if(flag == 1)
    {
        PAD_U_Y_MIN = 0.0;
        PAD_U_Y_MAX = 3000.0;
        PAD_L_Y_MIN = -10.0;
        PAD_L_Y_MAX = 10.0;
    }
    else if(flag == -1)
    {
        // do nothing
    }
    else if(flag == -2)
    {
        PAD_U_Y_MIN = 0.0;
        PAD_U_Y_MAX = 800.0;
        PAD_L_Y_MAX = 1.0;
        PAD_L_Y_MIN = 0.0;
    }
    else if(flag == -3)
    {
        PAD_U_Y_MIN = 0.0;
        PAD_U_Y_MAX = 800.0;
        PAD_L_Y_MAX = 1.1;
        PAD_L_Y_MIN = 0.9;
    }
    else if(flag == 4)
    {
        PAD_U_Y_MIN = 0.0;
        PAD_U_Y_MAX = 700.0;
        PAD_L_Y_MAX = 2.0;
        PAD_L_Y_MIN = 0.0;
    }

    h2->SetTitle("");
    //h2->SetMaximum(PAD_U_Y_MAX); // TODO: i re-enabled these 2 lines but were not needed
    //h2->SetMinimum(PAD_U_Y_MIN);
    h2->GetYaxis()->SetRangeUser(PAD_U_Y_MIN, PAD_U_Y_MAX);
    h2->SetStats(0);

    h3->SetTitle("");
    //h3->SetMaximum(PAD_L_Y_MAX);
    //h3->SetMinimum(PAD_L_Y_MIN);
    h3->GetYaxis()->SetRangeUser(PAD_L_Y_MIN, PAD_L_Y_MAX);
    h3->SetStats(0);          // No statistics on upper plot

    h3->SetLineColor(kBlack);
    h3->SetMarkerColor(kBlack);
    h3->SetMarkerStyle(20);
    //h3->Sumw2();



    ///////////////////////////////////////////////////////////////////////////
    // Y axis h1 plot settings
//    h1->GetYaxis()->SetNdivisions(505);
//    h1->GetYaxis()->SetTicks("+");
    h2->GetYaxis()->SetTickSize(0.0);
    
    h2->GetYaxis()->SetTitle("Normalized Spectrum");
    h2->GetYaxis()->SetTitleSize(20);
    h2->GetYaxis()->SetTitleFont(43);
    h2->GetYaxis()->SetTitleOffset(1.2);
   
    h2->GetYaxis()->SetLabelFont(63);
    h2->GetYaxis()->SetLabelSize(0);

    ///////////////////////////////////////////////////////////////////////////
    // X axis h1 plot settings
//    h1->GetXaxis()->SetTicks("+");
    h2->GetXaxis()->SetTickSize(0.05);

    h2->GetXaxis()->SetTitle("");
    h2->GetXaxis()->SetTitleSize(0);
    h2->GetXaxis()->SetTitleFont(43);
    h2->GetXaxis()->SetTitleOffset(1.0);

    h2->GetXaxis()->SetLabelFont(63);
    h2->GetXaxis()->SetLabelSize(0);

    ///////////////////////////////////////////////////////////////////////////
    // Y axis h3 plot settings
    h3->GetYaxis()->SetNdivisions(502);
//    h3->GetYaxis()->SetTicks("+");
    h3->GetYaxis()->SetTickSize(0.0);
//    h3->GetYaxis()->SetSize();

    TString xi_31_str;
    xi_31_str.Form("%3.3f", xi_31);
    //h3->GetYaxis()->SetTitle("ratio #xi^{31}_{"+xi_31_str+"}/#xi^{31}_{0.296}");
    h3->GetYaxis()->SetTitle("ratio #xi^{31}_{"+xi_31_str+"}/#xi^{31}_{0.0}");
    h3->GetYaxis()->SetTitleSize(20);
    h3->GetYaxis()->SetTitleFont(43);
    h3->GetYaxis()->SetTitleOffset(1.2);
   
    h3->GetYaxis()->SetLabelFont(63);
    h3->GetYaxis()->SetLabelSize(0);

    ///////////////////////////////////////////////////////////////////////////
    // X axis h3 plot settings
//    h3->GetXaxis()->SetTicks("+");
    h3->GetXaxis()->SetTickSize(0.1);
    
    h3->GetXaxis()->SetTitle("Single Electron Energy (T_{e}) [MeV]");
    h3->GetXaxis()->SetTitleSize(20);
    h3->GetXaxis()->SetTitleFont(43);
    h3->GetXaxis()->SetTitleOffset(4.0);

    h3->GetXaxis()->SetLabelFont(63);
    h3->GetXaxis()->SetLabelSize(15);
    h3->GetXaxis()->SetLabelOffset(0.02);



    // Define the Canvas
    TCanvas *c = new TCanvas("c", "canvas", 800, 600);
    c->SetBottomMargin(0);

    // Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1.0, 1.0);
    pad1->SetBottomMargin(0.0); // Upper and lower plot are joined
    pad1->SetGridx();         // Vertical grid
    pad1->SetGridy();         // Horizontal grid
    pad1->SetTicks(2, 2);
    pad1->Draw();             // Draw the upper pad: pad1

    // lower plot will be in pad2
    c->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.5);
    pad2->SetGridx(); // vertical grid
    pad2->SetGridy(); // horizontal grid
    pad2->SetTicks(2, 2); // ticks on both sides
    pad2->Draw();

    // Do not draw the Y axis label on the upper plot and redraw a small
    // axis instead, in order to avoid the first label (0) to be clipped.
    //h1->GetYaxis()->SetLabelSize(0);

    pad1->cd();               // pad1 becomes the current pad


//    TGaxis *axistmp = new TGaxis(2.0, 0.01, 3.0, 0.1, 0.0, 1.0, 510, "");
//    axistmp->Draw();

    

    // Y axis ratio plot settings


    h2->Draw("hist");
    h1->Draw("histsame");

    pad2->cd();       // pad2 becomes the current pad
    h3->Draw("EP");       // Draw the ratio plot

    pad1->cd();
    TLegend *leg = new TLegend(0.6, 0.4, 0.85, 0.85);
    //leg->SetNColumns(2);
    leg->AddEntry(h1, "#xi^{31}_{"+xi_31_str+"}", "F");
    //leg->AddEntry(h2, "#xi^{31}_{0.296}", "F");
    //leg->AddEntry(h2, "#xi^{31}_{0.0} (HSD)", "F");
    leg->AddEntry(h2, "#xi^{31}_{0.0}", "F");
    leg->SetBorderSize(5);
    leg->SetShadowColor(kGray + 2);
    leg->SetTextFont(62);
    leg->SetTextSize(0.07);
    leg->SetMargin(0.5);
    leg->Draw("BR");

    TGaxis *axis = new TGaxis(0.0, PAD_U_Y_MIN + 0.01, 0.0, PAD_U_Y_MAX, PAD_U_Y_MIN + 0.01, PAD_U_Y_MAX, 510, "");
    axis->SetLabelFont(43);
    axis->SetLabelSize(15);
    axis->Draw();

    TGaxis *axis2 = new TGaxis(5.0, PAD_U_Y_MIN + 0.01, 5.0, PAD_U_Y_MAX, PAD_U_Y_MIN + 0.01, PAD_U_Y_MAX, 510, "+");
    axis2->SetLabelFont(43);
    axis2->SetLabelSize(0);
    axis2->Draw();

    pad2->cd();

    TGaxis *axis3 = new TGaxis(0.0, PAD_L_Y_MIN, 0.0, PAD_L_Y_MAX - 0.01, PAD_L_Y_MIN, PAD_L_Y_MAX - 0.01, 502, "");
    axis3->SetLabelFont(43);
    axis3->SetLabelSize(15);
    axis3->Draw();

    TGaxis *axis4 = new TGaxis(5.0, PAD_L_Y_MIN, 5.0, PAD_L_Y_MAX - 0.01, PAD_L_Y_MIN, PAD_L_Y_MAX - 0.01, 502, "+");
    axis4->SetLabelFont(43);
    axis4->SetLabelSize(0);
    axis4->Draw();



    if(saveas_filename.size() > 0)
    {
        //TODO:
        std::size_t length = saveas_filename.size();
        if(length >= 2)
        {
            if((saveas_filename[length - 1] == '*') && (saveas_filename[length - 2] == '.'))
            {
                std::string base_name = saveas_filename.substr(0, length - 2);
                std::vector<std::string> extensions;
                extensions.push_back("C");
                extensions.push_back("png");
                extensions.push_back("eps");
                extensions.push_back("pdf");
                for(auto it{extensions.begin()}; it != extensions.end(); ++ it)
                {
                    //std::string name = base_name + "_" + std::string(i_str) + "." + *it;
                    std::string name = base_name + "." + *it;
                    std::cout << "saving as " << name << std::endl;
                    //c[channel]->SaveAs(name.c_str());
                    c->SaveAs(name.c_str());
                }
            }
            else
            {
                c->SaveAs(saveas_filename.c_str());
            }
        }
        else
        {
            // TODO: bug: add i_str
            c->SaveAs(saveas_filename.c_str());
        }
    }
    else
    {
        std::string filename = std::string("c_ratio");
        if(xi_31 >= 0.0)
        {
            filename += "+";
        }
        filename += std::string(xi_31_str) + std::string(".png");
        c->SaveAs(filename.c_str());
        //c->SaveAs("c_outputdiff.eps");
        //c->SaveAs("c_outputdiff.pdf");
        //c->SaveAs("c_outputdiff.C");

        //TFile *fcout = new TFile("c_outputdiff.root", "recreate");
        //c->Write();
    }
}



// NOTE: obsolete function
void draw_outputdiff(const double *const AdjustActs,
                     const double xi_31_orig,
                     const std::string &saveas_filename,
                     int flag)
{
#if 0
    // 2020-06-17
    const double xi_31_init_P1 = paramInitValueP1Map[1];
    const double xi_31_init_P2 = paramInitValueP2Map[1];
    double xi_31_init = 0.0;
    if(thePhase == 0)
    {
        xi_31_init = xi_31_init_P1;
    }
    else if(thePhase == 1)
    {
        xi_31_init = xi_31_init_P2;
    }
    else
    {
        throw "invalid value for thePhase";
    }
    //const double xi_31_fit = AdjustActs[1] * xi_31_init;
    const double xi_31_fit = AdjustActs[1];
    //const double xi_31_baseline{0.296};
    //const double xi_31_baseline{0.0};

    double xi_31 = xi_31_orig;
    TH1D *hTotalE_orig = nullptr;
    TH1D *hSingleEnergy_orig = nullptr;
    TH1D *hHighEnergy_orig = nullptr;
    TH1D *hLowEnergy_orig = nullptr;
    TH1D *hEnergySum_orig = nullptr;
    TH1D *hEnergyDiff_orig = nullptr;
    TH2D *hHighLowEnergy_orig = nullptr;
    TH1D *hWeight_orig = nullptr;
    reweight_apply(hTotalE_orig,
                   hSingleEnergy_orig,
                   hHighEnergy_orig,
                   hLowEnergy_orig,
                   hEnergySum_orig,
                   hEnergyDiff_orig,
                   hHighLowEnergy_orig,
                   hWeight_orig,
                   "Nd150_2eNg_output_truth_postprocessed_small.root",
                   xi_31,
                   xi_31_baseline,
                   h_nEqNull,
                   h_nEqTwo,
                   psiN0,
                   psiN2,
                   bb_Q);

    xi_31 = xi_31_fit;
    TH1D *hTotalE_fit = nullptr;
    TH1D *hSingleEnergy_fit = nullptr;
    TH1D *hHighEnergy_fit = nullptr;
    TH1D *hLowEnergy_fit = nullptr;
    TH1D *hEnergySum_fit = nullptr;
    TH1D *hEnergyDiff_fit = nullptr;
    TH2D *hHighLowEnergy_fit = nullptr;
    TH1D *hWeight_fit = nullptr;
    reweight_apply(hTotalE_fit,
                   hSingleEnergy_fit,
                   hHighEnergy_fit,
                   hLowEnergy_fit,
                   hEnergySum_fit,
                   hEnergyDiff_fit,
                   hHighLowEnergy_fit,
                   hWeight_fit,
                   "Nd150_2eNg_output_truth_postprocessed_small.root",
                   xi_31,
                   xi_31_baseline,
                   h_nEqNull,
                   h_nEqTwo,
                   psiN0,
                   psiN2,
                   bb_Q);

/*
//    TCanvas *cWeight_orig = new TCanvas("cWeight_orig", "cWeight_orig");
    TCanvas *cWeight = new TCanvas("cWeight", "cWeight");
    hWeight_orig->SetLineColor(kBlue);
    hWeight_orig->Draw("hist");
//    TCanvas *cWeight_fit = new TCanvas("cWeight_fit", "cWeight_fit");
    hWeight_fit->SetLineColor(kRed);
    hWeight_fit->Draw("histsame");
*/



    hTotalE_fit->Scale(AdjustActs[0]);
    hSingleEnergy_fit->Scale(AdjustActs[0]);
    hLowEnergy_fit->Scale(AdjustActs[0]);
    hHighEnergy_fit->Scale(AdjustActs[0]);
    hHighLowEnergy_fit->Scale(AdjustActs[0]);

    // TODO: 150 Nd is scaled by AdjustActs[0] later in this code.
    // check!
    
/*    hTotalE_orig->Sumw2();
    hSingleEnergy_orig->Sumw2();
    hLowEnergy_orig->Sumw2();
    hHighEnergy_orig->Sumw2();
    hHighLowEnergy_orig->Sumw2();
    
    hTotalE_fit->Sumw2();
    hSingleEnergy_fit->Sumw2();
    hLowEnergy_fit->Sumw2();
    hHighEnergy_fit->Sumw2();
    hHighLowEnergy_fit->Sumw2();*/

    hTotalE_orig->SetLineWidth(2);
    hSingleEnergy_orig->SetLineWidth(2);
    hLowEnergy_orig->SetLineWidth(2);
    hHighEnergy_orig->SetLineWidth(2);
    hHighLowEnergy_orig->SetLineWidth(2);

    hTotalE_fit->SetLineWidth(2);
    hSingleEnergy_fit->SetLineWidth(2);
    hLowEnergy_fit->SetLineWidth(2);
    hHighEnergy_fit->SetLineWidth(2);
    hHighLowEnergy_fit->SetLineWidth(2);

    //hTotalE_orig->SetHatchesLineWidth(2);
    //hSingleEnergy_orig->SetHatchesLineWidth(2);
    //hLowEnergy_orig->SetHatchesLineWidth(2);
    //hHighEnergy_orig->SetHatchesLineWidth(2);
    //hHighLowEnergy_orig->SetHatchesLineWidth(2);

    //hTotalE_fit->SetHatchesLineWidth(2);
    //hSingleEnergy_fit->SetHatchesLineWidth(2);
    //hLowEnergy_fit->SetHatchesLineWidth(2);
    //hHighEnergy_fit->SetHatchesLineWidth(2);
    //hHighLowEnergy_fit->SetHatchesLineWidth(2);
    
    hTotalE_orig->SetFillStyle(3002);
    hSingleEnergy_orig->SetFillStyle(3002);
    hHighEnergy_orig->SetFillStyle(3002);
    hLowEnergy_orig->SetFillStyle(3002);
    hEnergySum_orig->SetFillStyle(3002);
    hEnergyDiff_orig->SetFillStyle(3002);
    hHighLowEnergy_orig->SetFillStyle(3002);

    hTotalE_fit->SetFillStyle(3453);
    hSingleEnergy_fit->SetFillStyle(3453);
    hHighEnergy_fit->SetFillStyle(3453);
    hLowEnergy_fit->SetFillStyle(3453);
    hEnergySum_fit->SetFillStyle(3453);
    hEnergyDiff_fit->SetFillStyle(3453);
    hHighLowEnergy_fit->SetFillStyle(3453);

    hTotalE_orig->SetFillColor(kBlue);
    hSingleEnergy_orig->SetFillColor(kBlue);
    hHighEnergy_orig->SetFillColor(kBlue);
    hLowEnergy_orig->SetFillColor(kBlue);
    hEnergySum_orig->SetFillColor(kBlue);
    hEnergyDiff_orig->SetFillColor(kBlue);
    hHighLowEnergy_orig->SetFillColor(kBlue);

    hTotalE_fit->SetFillColor(kMagenta);
    hSingleEnergy_fit->SetFillColor(kMagenta);
    hHighEnergy_fit->SetFillColor(kMagenta);
    hLowEnergy_fit->SetFillColor(kMagenta);
    hEnergySum_fit->SetFillColor(kMagenta);
    hEnergyDiff_fit->SetFillColor(kMagenta);
    hHighLowEnergy_fit->SetFillColor(kMagenta);

    hTotalE_orig->SetLineColor(kBlue);
    hSingleEnergy_orig->SetLineColor(kBlue);
    hHighEnergy_orig->SetLineColor(kBlue);
    hLowEnergy_orig->SetLineColor(kBlue);
    hEnergySum_orig->SetLineColor(kBlue);
    hEnergyDiff_orig->SetLineColor(kBlue);
    hHighLowEnergy_orig->SetLineColor(kBlue);

    hTotalE_fit->SetLineColor(kMagenta);
    hSingleEnergy_fit->SetLineColor(kMagenta);
    hHighEnergy_fit->SetLineColor(kMagenta);
    hLowEnergy_fit->SetLineColor(kMagenta);
    hEnergySum_fit->SetLineColor(kMagenta);
    hEnergyDiff_fit->SetLineColor(kMagenta);
    hHighLowEnergy_fit->SetLineColor(kMagenta);

    //hSingleEnergy_orig->SetMaximum(800.0);
//    hSingleEnergy_orig->SetMaximum(1.0);

    // normalize
    //Double_t hSingleEnergy_orig_integral = hSingleEnergy_orig->Integral();
    //Double_t hSingleEnergy_fit_integral = hSingleEnergy_fit->Integral();
    Double_t hSingleEnergy_orig_integral = 0.0;
    Double_t hSingleEnergy_fit_integral = 0.0;
    for(Int_t ix{1}; ix < hSingleEnergy_orig->GetNbinsX(); ++ ix)
    {
        hSingleEnergy_orig_integral += std::abs(hSingleEnergy_orig->GetBinContent(ix));
        hSingleEnergy_fit_integral += std::abs(hSingleEnergy_fit->GetBinContent(ix));
    }
    //hSingleEnergy_orig->Scale(1.0 / std::abs(hSingleEnergy_orig_integral));
    //hSingleEnergy_fit->Scale(1.0 / std::abs(hSingleEnergy_fit_integral));
//    hSingleEnergy_orig->Scale(1.0 / hSingleEnergy_orig_integral);
//    hSingleEnergy_fit->Scale(1.0 / hSingleEnergy_fit_integral);
    hSingleEnergy_fit->Scale(hSingleEnergy_orig_integral / hSingleEnergy_fit_integral);

//    TH1D *hSingleEnergy_ratio = (TH1D*)hSingleEnergy_orig->Clone("hSingleEnergy_ratio");
//    hSingleEnergy_ratio->Divide(hSingleEnergy_fit);

    ratioplotOld(hSingleEnergy_fit, hSingleEnergy_orig, xi_31_fit, "", flag);

/*
    // plot title
    hSingleEnergy_orig->SetTitle("");
    hSingleEnergy_ratio->SetTitle("");


    // orig X title
    hSingleEnergy_orig->GetXaxis()->SetTitleSize(15);
    hSingleEnergy_orig->GetXaxis()->SetTitleFont(63);
    hSingleEnergy_orig->GetXaxis()->SetTitle("");
    hSingleEnergy_orig->GetXaxis()->SetTitleOffset(1.0);

    // orig Y title
    hSingleEnergy_orig->GetYaxis()->SetTitleSize(15);
    hSingleEnergy_orig->GetYaxis()->SetTitleFont(63);
    hSingleEnergy_orig->GetYaxis()->SetTitle("some title here");
    hSingleEnergy_orig->GetYaxis()->SetTitleOffset(1.55);
//    hSingleEnergy_orig->GetYaxis()->SetLabelSize(0.);


    // orig X ticks
    hSingleEnergy_orig->GetXaxis()->SetTicks("-");
    hSingleEnergy_orig->GetXaxis()->SetTickSize(0.05);

    // orig Y ticks
    hSingleEnergy_orig->GetYaxis()->SetTicks("-");
    hSingleEnergy_orig->GetYaxis()->SetTickSize(0.05);
    
    // orig X label

    // orig Y label

    //TGaxis *axis = new TGaxis(0.0, 0.0, 5.0, 0.0, 0.0, 5.0, 5, "");
    //axis->SetLabelFont(63);
    //axis->SetLabelSize(15);
    //axis->Draw();



    // ratio line/marker min/max
    hSingleEnergy_ratio->SetLineColor(kBlack);
    hSingleEnergy_ratio->SetMarkerStyle(20);
    hSingleEnergy_ratio->SetMinimum(0.0);
    hSingleEnergy_ratio->SetMaximum(2.0);
    hSingleEnergy_ratio->Sumw2();
//    hSingleEnergy_ratio->SetStats(0);

    //hSingleEnergy_ratio->SetNdivisions(505);
//    hSingleEnergy_ratio->SetTitleSize(20);
//    hSingleEnergy_ratio->SetTitleFont(63);
    //hSingleEnergy_ratio->SetTitleOffset(1.55);
//    hSingleEnergy_ratio->SetLabelSize(15);


    // ratio X title
    //hSingleEnergy_ratio->GetXaxis()->SetTitleSize(15);
    hSingleEnergy_ratio->GetXaxis()->SetTitleFont(63);
    hSingleEnergy_ratio->GetXaxis()->SetTitleOffset(0.8);
    hSingleEnergy_ratio->GetXaxis()->SetTitle("Single Electron Energy (T_e) [MeV]");

    // ratio X label
    hSingleEnergy_ratio->GetXaxis()->SetLabelFont(63);
    hSingleEnergy_ratio->GetXaxis()->SetLabelSize(15);
    hSingleEnergy_ratio->GetXaxis()->SetLabelOffset(-1.0);
    
    // ratio X ticks
    hSingleEnergy_ratio->GetXaxis()->SetTicks("+");

    // ratio Y title
    //hSingleEnergy_ratio->GetYaxis()->SetTitleSize(15);
    hSingleEnergy_ratio->GetYaxis()->SetTitleFont(63);
    hSingleEnergy_ratio->GetYaxis()->SetTitleOffset(1.0);
    hSingleEnergy_ratio->GetYaxis()->SetTitle("ratio orig/fit");

    // ratio Y label
    hSingleEnergy_ratio->GetYaxis()->SetLabelFont(63);
    hSingleEnergy_ratio->GetYaxis()->SetLabelSize(15);
    hSingleEnergy_ratio->GetYaxis()->SetLabelOffset(1.0);

    // ratio Y ticks
    hSingleEnergy_ratio->GetYaxis()->SetTicks("-");
    hSingleEnergy_ratio->GetYaxis()->SetNdivisions(203);




    TCanvas *c = new TCanvas("c_outputdiff", "c_outputdiff");
    c->SetTopMargin(0);
    c->SetBottomMargin(0);
    c->SetLeftMargin(0);
    c->SetRightMargin(0);

    TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
    pad1->SetBottomMargin(0.0);
    pad1->SetTopMargin(0.05);
    pad1->SetLeftMargin(0.10);
    pad1->SetRightMargin(0.05);
    pad1->SetTicks(1, 1);
    //pad1->SetGridx();
    pad1->Draw();
    pad1->cd();

    hSingleEnergy_orig->Draw("hist");
    hSingleEnergy_fit->Draw("histsame");

    c->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.3);
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.4);
    pad2->SetLeftMargin(0.10);
    pad2->SetRightMargin(0.05);
    //pad2->SetBottomMargin(0.5);
    //pad2->SetGridx();
    pad2->SetTicks(1, 1);
    pad2->Draw();
    pad2->cd();

    hSingleEnergy_ratio->Draw("PE");

*/
    
    const Int_t channel = 1; // SingleEnergy
    TH1D *data1D = (TH1D*)allDataSamples1D->At(channel)->Clone();
    
    // iterate through all backgrounds except 150Nd and subtract
    // value from data1D
    for(int j = 0; j < allMCSamples1D[channel]->GetEntries(); ++ j)
    {
        TH1D *tmpHist = (TH1D*)allMCSamples1D[channel]->At(j)->Clone();

        TString tmpName = tmpHist->GetName();

        int which_param = -1;
        bool found_param = false;

        
        //double activity_scale_branching_ratio = 1.0;

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
            Double_t activity_scale = AdjustActs[which_param]; // * activity_scale_branching_ratio;
            tmpHist->Scale(activity_scale);

            // TODO: already scaled by AdjustActs[0] in the case of 150 Nd ?

                
            if(tmpHist->Integral() > 0)
            {
                TString hname = tmpHist->GetName();

                // 150 Nd
                if(hname.Contains("150"))
                {
                    continue;
                }

                for(int bin_ix{1}; bin_ix <= tmpHist->GetNbinsX(); ++ bin_ix)
                {
                    double_t new_content = data1D->GetBinContent(bin_ix);
                    new_content -= tmpHist->GetBinContent(bin_ix);
                    data1D->SetBinContent(bin_ix, new_content);
                }
            }
            else
            {
                //std::cout << "not adding to stack, Integral() <= 0: " << tmpHist->GetName() << std::endl;
            }
        }
        else
        {
            std::cout << "error could not find histogram: tmpName=" << tmpName << std::endl;
        } 
    }

    data1D->SetLineWidth(2);
    data1D->SetMarkerStyle(20);
    data1D->SetMarkerSize(1.0);
    //data1D->Draw("PEsame");
   

    // save to file
    /*
    if(saveas_filename.size() > 0)
    {
        std::size_t length = saveas_filename.size();
        if(length >= 2)
        {
            if((saveas_filename[length - 1] == '*') && (saveas_filename[length - 2] == '.'))
            {
                std::string base_name = saveas_filename.substr(0, length - 2);
                std::vector<std::string> extensions;
                extensions.push_back("C");
                extensions.push_back("png");
                extensions.push_back("eps");
                extensions.push_back("pdf");
                for(auto it{extensions.begin()}; it != extensions.end(); ++ it)
                {
                    //std::string name = base_name + "_" + std::string(i_str) + "." + *it;
                    std::string name = base_name + "." + *it;
                    std::cout << "saving as " << name << std::endl;
                    //c[channel]->SaveAs(name.c_str());
                    c->SaveAs(name.c_str());
                }
            }
            else
            {
                c->SaveAs(saveas_filename.c_str());
            }
        }
        else
        {
            // TODO: bug: add i_str
            c->SaveAs(saveas_filename.c_str());
        }
    }
    else
    {
        c->SaveAs("c_outputdiff.png");
        //c->SaveAs("c_outputdiff.eps");
        //c->SaveAs("c_outputdiff.pdf");
        //c->SaveAs("c_outputdiff.C");

        //TFile *fcout = new TFile("c_outputdiff.root", "recreate");
        //c->Write();
    }
    */

    TString xi_31_str_plus;
    xi_31_str_plus.Form("%+3.3f", xi_31_fit);
    TString xi_31_str;
    xi_31_str.Form("%3.3f", xi_31_fit);

    //hHighLowEnergy_orig->Scale(0.1);
    //std::cout << "E(10,10)=" << hHighLowEnergy_orig->GetBinError(10, 10) << std::endl;
    //hHighLowEnergy_orig->Scale(1.0 / 0.1);
    //std::cout << "E(10,10)=" << hHighLowEnergy_orig->GetBinError(10, 10) << std::endl;


    for(Int_t ix{1}; ix <= hHighLowEnergy_orig->GetNbinsX(); ++ ix)
    {
        for(Int_t iy{1}; iy <= hHighLowEnergy_orig->GetNbinsY(); ++ iy)
        {
            if(hHighLowEnergy_orig->GetBinContent(ix, iy) == 0.0)
            {
                hHighLowEnergy_orig->SetBinContent(ix, iy, -std::numeric_limits<double>::infinity());
            }

            if(hHighLowEnergy_fit->GetBinContent(ix, iy) == 0.0)
            {
                hHighLowEnergy_fit->SetBinContent(ix, iy, -std::numeric_limits<double>::infinity());
            }
        }
    }


    TH2D *hHighLowEnergy_err = (TH2D*)hHighLowEnergy_orig->Clone();
    for(Int_t ix{1}; ix <= hHighLowEnergy_err->GetNbinsX(); ++ ix)
    {
        for(Int_t iy{1}; iy <= hHighLowEnergy_err->GetNbinsY(); ++ iy)
        {
            double content = hHighLowEnergy_orig->GetBinContent(ix, iy);
            double error = hHighLowEnergy_orig->GetBinError(ix, iy);

            if(content > 0.0)
            {
                hHighLowEnergy_err->SetBinContent(ix, iy, error);
                //std::cout << "ix=" << ix << " iy=" << iy << " content=" << content << " sqrt(content)=" << std::sqrt(content) << " error=" << error << std::endl;
            }
            else
            {
                hHighLowEnergy_err->SetBinContent(ix, iy, -std::numeric_limits<double>::infinity());
            }
        }
    }


    TH2D *hHighLowEnergy_div = (TH2D*)hHighLowEnergy_orig->Clone();
    //hHighLowEnergy_div->Sumw2();
    //hHighLowEnergy_div->Divide(hHighLowEnergy_orig);
    //hHighLowEnergy_div->Scale(hHighLowEnergy_div->Integral());
    //hHighLowEnergy_div->Scale(0.1);
    double hHighLowEnergy_div_min = std::numeric_limits<double>::infinity();
    double hHighLowEnergy_div_max = std::numeric_limits<double>::infinity();
    for(Int_t ix{1}; ix <= hHighLowEnergy_div->GetNbinsX(); ++ ix)
    {
        for(Int_t iy{1}; iy <= hHighLowEnergy_div->GetNbinsY(); ++ iy)
        {
            double content_fit = hHighLowEnergy_fit->GetBinContent(ix, iy);
            double content_orig = hHighLowEnergy_orig->GetBinContent(ix, iy);

            //if(content_orig != 0.0)
            if(content_orig != -std::numeric_limits<double>::infinity())
            {
                double new_content = content_fit / content_orig;

                hHighLowEnergy_div->SetBinContent(ix, iy, new_content);
                //std::cout << "ix=" << ix << " iy=" << iy << " fit=" << content_fit << " orig=" << content_orig << " div=" << content_fit / content_orig << std::endl;

                if(hHighLowEnergy_div_min == std::numeric_limits<double>::infinity())
                {
                    hHighLowEnergy_div_min = new_content;
                    hHighLowEnergy_div_max = new_content;
                }
                else
                {
                    if(new_content > hHighLowEnergy_div_max)
                    {
                        hHighLowEnergy_div_max = new_content;
                    }
                    
                    if(new_content < hHighLowEnergy_div_min)
                    {
                        hHighLowEnergy_div_min = new_content;
                    }
                }
            }
            else
            {
                hHighLowEnergy_div->SetBinContent(ix, iy, -std::numeric_limits<double>::infinity());
            }
        }
    }


    TH2D *hHighLowEnergy_sig = (TH2D*)hHighLowEnergy_orig->Clone();
    double hHighLowEnergy_sig_min = std::numeric_limits<double>::infinity();
    double hHighLowEnergy_sig_max = std::numeric_limits<double>::infinity();
    for(Int_t ix{1}; ix <= hHighLowEnergy_sig->GetNbinsX(); ++ ix)
    {
        for(Int_t iy{1}; iy <= hHighLowEnergy_sig->GetNbinsY(); ++ iy)
        {
            Double_t content1 = hHighLowEnergy_orig->GetBinContent(ix, iy);
            Double_t error1 = hHighLowEnergy_orig->GetBinError(ix, iy);
            Double_t content2 = hHighLowEnergy_fit->GetBinContent(ix, iy);

            //if(error != 0.0)
            if(content1 != -std::numeric_limits<double>::infinity())
            {
                //double new_content = std::pow((content1 - content2) / error2, 2.0);
                double new_content = (content2 - content1) / error1;
                hHighLowEnergy_sig->SetBinContent(ix, iy, new_content);


                /*
                std::cout << "ix=" << ix << " iy=" << iy
                          << " content=" << content << " error=" << error
                          << " content1=" << content1 << " content2=" << content2 << " error2=" << error2
                          << " -> " << (content1 - content2) / error2
                          << " : ^2-> " << new_content << std::endl;
                */


                if(hHighLowEnergy_sig_min == std::numeric_limits<double>::infinity())
                {
                    hHighLowEnergy_sig_min = new_content;
                    hHighLowEnergy_sig_max = new_content;
                }
                else
                {
                    if(new_content > hHighLowEnergy_sig_max)
                    {
                        hHighLowEnergy_sig_max = new_content;
                    }
                    
                    if(new_content < hHighLowEnergy_sig_min)
                    {
                        hHighLowEnergy_sig_min = new_content;
                    }
                }
            }
            else
            {
                hHighLowEnergy_sig->SetBinContent(ix, iy, -std::numeric_limits<double>::infinity());
            }
        }
    }

    //std::cout << "make sig_2" << std::endl;

    TH2D *hHighLowEnergy_sig_2 = (TH2D*)hHighLowEnergy_orig->Clone();
    double hHighLowEnergy_sig_2_min = std::numeric_limits<double>::infinity();
    double hHighLowEnergy_sig_2_max = std::numeric_limits<double>::infinity();
    for(Int_t ix{1}; ix <= hHighLowEnergy_sig_2->GetNbinsX(); ++ ix)
    {
        for(Int_t iy{1}; iy <= hHighLowEnergy_sig_2->GetNbinsY(); ++ iy)
        {
            //Double_t content = hHighLowEnergy_sig->GetBinContent(ix, iy);
            //Double_t error = hHighLowEnergy_orig->GetBinError(ix, iy);

            Double_t content1 = hHighLowEnergy_orig->GetBinContent(ix, iy);
            Double_t error1 = hHighLowEnergy_orig->GetBinError(ix, iy);
            Double_t content2 = hHighLowEnergy_fit->GetBinContent(ix, iy);

            //if(error != 0.0)
            if(content1 != -std::numeric_limits<double>::infinity())
            {
                double new_content = std::pow((content2 - content1) / error1, 2.0);
                hHighLowEnergy_sig_2->SetBinContent(ix, iy, new_content);

                if(hHighLowEnergy_sig_2_min == std::numeric_limits<double>::infinity())
                {
                    hHighLowEnergy_sig_2_min = new_content;
                    hHighLowEnergy_sig_2_max = new_content;
                }
                else
                {
                    if(new_content > hHighLowEnergy_sig_2_max)
                    {
                        hHighLowEnergy_sig_2_max = new_content;
                    }
                    
                    if(new_content < hHighLowEnergy_sig_2_min)
                    {
                        hHighLowEnergy_sig_2_min = new_content;
                    }
                }
            }
            else
            {
                hHighLowEnergy_sig_2->SetBinContent(ix, iy, -std::numeric_limits<double>::infinity());
            }
        }
    }

    TH2D *hHighLowEnergy_diff = (TH2D*)hHighLowEnergy_orig->Clone(); 
    for(Int_t ix{1}; ix <= hHighLowEnergy_diff->GetNbinsX(); ++ ix)
    {
        for(Int_t iy{1}; iy <= hHighLowEnergy_diff->GetNbinsY(); ++ iy)
        {
            Double_t content1 = hHighLowEnergy_orig->GetBinContent(ix, iy);
            Double_t content2 = hHighLowEnergy_fit->GetBinContent(ix, iy);
            //Double_t error = hHighLowEnergy_orig->GetBinError(ix, iy);
            //if(error != 0.0)
            if(content1 != -std::numeric_limits<double>::infinity())
            {
                hHighLowEnergy_diff->SetBinContent(ix, iy, content2 - content1);

            }
            else
            {
                hHighLowEnergy_diff->SetBinContent(ix, iy, -std::numeric_limits<double>::infinity());
                //hHighLowEnergy_diff->SetBinContent(ix, iy, 0.0);
            }
            /*
            if(content1 != 0.0 && error != 0.0)
            {
                hHighLowEnergy_diff->SetBinContent(ix, iy, content2 - content1);
            }
            else if(content2 != 0.0 && error != 0.0)
            {
                hHighLowEnergy_diff->SetBinContent(ix, iy, content2);
            }
            */
        }
    }
    //hHighLowEnergy_diff->Add(hHighLowEnergy_orig, -1.0);

    //double hHighLowEnergy_div_min = hHighLowEnergy_div->GetMinimumStored();
    //double hHighLowEnergy_div_max = hHighLowEnergy_div->GetMaximumStored();
    /*
    double hHighLowEnergy_div_min = std::numeric_limits<double>::infinity();
    double hHighLowEnergy_div_max = std::numeric_limits<double>::infinity();
    for(Int_t ix{1}; ix <= hHighLowEnergy_div->GetNbinsX(); ++ ix)
    {
        for(Int_t iy{1}; iy <= hHighLowEnergy_div->GetNbinsY(); ++ iy)
        {
            Double_t content = hHighLowEnergy_div->GetBinContent(ix, iy);
            Double_t error = hHighLowEnergy_div->GetBinError(ix, iy);
            if((1.0e-05 < content) && (content < 1.0e+05))
            {
                if(hHighLowEnergy_div_min == std::numeric_limits<double>::infinity())
                {
                    hHighLowEnergy_div_min = content;
                    hHighLowEnergy_div_max = content;
                }
                else
                {
                    if(content > hHighLowEnergy_div_max)
                    {
                        hHighLowEnergy_div_max = content;
                    }
                    
                    if(content < hHighLowEnergy_div_min)
                    {
                        hHighLowEnergy_div_min = content;
                    }
                }
            }
            else
            {
                //hHighLowEnergy_diff->SetBinContent(ix, iy, -std::numeric_limits<double>::infinity());
            }
        }
    }
    */

    TCanvas *cHighLowEnergy_diff = new TCanvas("cHighLowEnergy_diff", "cHighLowEnergy_diff", 800, std::sqrt(2.0) * 800);

    // top left
    cHighLowEnergy_diff->cd();
    TPad *pad1 = new TPad("pad1", "pad1", 0.0, 3.0/4.0, 0.5, 4.0/4.0);
    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.2);
    pad1->SetBottomMargin(0.15);
    pad1->Draw();
    
    // mid left
    cHighLowEnergy_diff->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.0, 2.0/4.0, 0.5, 3.0/4.0);
    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.2);
    pad2->SetBottomMargin(0.15);
    pad2->Draw();
    
    // top right
    cHighLowEnergy_diff->cd();
    TPad *pad3 = new TPad("pad3", "pad3", 0.5, 3.0/4.0, 1.0, 4.0/4.0);
    pad3->SetLeftMargin(0.15);
    pad3->SetRightMargin(0.25);
    pad3->SetBottomMargin(0.15);
    pad3->Draw();

    // mid right
    cHighLowEnergy_diff->cd();
    TPad *pad4 = new TPad("pad4", "pad4", 0.5, 2.0/4.0, 1.0, 3.0/4.0);
    pad4->SetLeftMargin(0.15);
    pad4->SetRightMargin(0.25);
    pad4->SetBottomMargin(0.15);
    pad4->Draw();

    cHighLowEnergy_diff->cd();
    //TPad *pad5 = new TPad("pad5", "pad5", 0.0, 0.0, 1.0, 0.5);
    TPad *pad5 = new TPad("pad5", "pad5", 0.0, 0.0, 1.0, 0.5);
    pad5->SetTopMargin(0.05);
    pad5->SetRightMargin(0.15);
    pad5->Draw();

    //pad5->cd();
    cHighLowEnergy_diff->cd();
    TPad *pad5_2 = new TPad("pad5_2", "pad5_2", 0.3, 0.475 - (0.85 - 0.3) / 2.0, 0.85, 0.475);
    pad5_2->SetTopMargin(0.05);
    pad5_2->SetLeftMargin(0.15);
    pad5_2->SetRightMargin(0.25);
    pad5_2->SetBottomMargin(0.2);
    //pad5_2->SetBorderMode(-1);
    //pad5_2->SetBorderSize(3);
    pad5_2->Draw();

    cHighLowEnergy_diff->cd();
    TPad *pad3_2 = new TPad("pad3_2", "pad3_2", 0.68, 0.79, 0.87, 0.94);
    pad3_2->SetTopMargin(0.15);
    pad3_2->SetLeftMargin(0.25);
    pad3_2->SetRightMargin(0.1);
    pad3_2->SetBottomMargin(0.2);
    pad3_2->Draw();


    // TODO
    // fix these problems
    //
    // empty bins still being counted in chisquare [FIXED]
    // not sure if -inf values are counted in chisquare [THEY ARE NOT]
    // make "empty bins" -INF [DONE]
    // ensure consistency between histograms by tagging empty
    // bins with -INF [DONE]


    const int FONT = 43;
    const double LABEL_SIZE_S = 17.5;
    const double LABEL_SIZE_L = 20.0;
    const double TITLE_SIZE_S = 17.5;
    const double TITLE_SIZE_L = 20.0;
    const double TITLE_OFFSET_S_Z_SMALL = 3.5;
    const double TITLE_OFFSET_S_Z_MEDIUM = 4.5;
    const double TITLE_OFFSET_S_Z_LARGE = 5.5;
    const double TITLE_OFFSET_S_XY = 4.0;
    const double LATEX_X_LEFT = 0.75;
    const double LATEX_X_RIGHT = 0.7;
    const double LATEX_Y = 0.85;
    const int LATEX_FONT = FONT;
    const double LATEX_FONT_SIZE = 15.0;

    int draworder[] = {5, 4, 3, 2, 1};
    for(int i = 0; i < sizeof(draworder) / sizeof(draworder[0]); ++ i)
    {
        int order = draworder[i];

        if(order == 1)
        {
            // top left
            pad1->cd();
            hHighLowEnergy_orig->SetTitle("");
            hHighLowEnergy_orig->SetContour(10000);
            hHighLowEnergy_orig->Draw("colz");
            hHighLowEnergy_orig->GetXaxis()->SetTitleFont(FONT);
            hHighLowEnergy_orig->GetXaxis()->SetTitleSize(TITLE_SIZE_S);
            hHighLowEnergy_orig->GetXaxis()->SetTitleOffset(TITLE_OFFSET_S_XY);
            hHighLowEnergy_orig->GetXaxis()->SetLabelFont(FONT);
            hHighLowEnergy_orig->GetXaxis()->SetLabelSize(LABEL_SIZE_S);
            hHighLowEnergy_orig->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
            hHighLowEnergy_orig->GetXaxis()->SetNdivisions(205);

            hHighLowEnergy_orig->GetYaxis()->SetTitleFont(FONT);
            hHighLowEnergy_orig->GetYaxis()->SetTitleSize(TITLE_SIZE_S);
            hHighLowEnergy_orig->GetYaxis()->SetTitleOffset(TITLE_OFFSET_S_XY);
            hHighLowEnergy_orig->GetYaxis()->SetLabelFont(FONT);
            hHighLowEnergy_orig->GetYaxis()->SetLabelSize(LABEL_SIZE_S);
            hHighLowEnergy_orig->GetYaxis()->SetTitle("High Energy Electron [MeV]");
            hHighLowEnergy_orig->GetYaxis()->SetNdivisions(205);

            hHighLowEnergy_orig->GetZaxis()->SetTitleFont(FONT);
            hHighLowEnergy_orig->GetZaxis()->SetTitleSize(TITLE_SIZE_S);
            hHighLowEnergy_orig->GetZaxis()->SetTitleOffset(TITLE_OFFSET_S_Z_MEDIUM);
            hHighLowEnergy_orig->GetZaxis()->SetLabelFont(FONT);
            hHighLowEnergy_orig->GetZaxis()->SetLabelSize(LABEL_SIZE_S);
            hHighLowEnergy_orig->GetZaxis()->SetTitle("Weight");

            TString st;
            //st = "h(#xi_{31} = 0.296)";
            st = "h(#xi_{31} = 0.0) (HSD)";
            TLatex tl;
            tl.SetNDC();
            tl.SetTextAlign(33);
            tl.SetTextFont(LATEX_FONT);
            tl.SetTextSize(LATEX_FONT_SIZE);
            tl.DrawLatex(LATEX_X_LEFT, LATEX_Y, st);
            tl.DrawLatex(0.25, LATEX_Y, "(a)");
        }
        else if(order == 2)
        {
            // top right
            pad3->cd();
            hHighLowEnergy_fit->SetTitle("");
            hHighLowEnergy_fit->SetContour(10000);
            hHighLowEnergy_fit->Draw("colz");
            hHighLowEnergy_fit->GetXaxis()->SetTitleFont(FONT);
            hHighLowEnergy_fit->GetXaxis()->SetTitleSize(TITLE_SIZE_S);
            hHighLowEnergy_fit->GetXaxis()->SetTitleOffset(TITLE_OFFSET_S_XY);
            hHighLowEnergy_fit->GetXaxis()->SetLabelFont(FONT);
            hHighLowEnergy_fit->GetXaxis()->SetLabelSize(LABEL_SIZE_S);
            hHighLowEnergy_fit->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
            hHighLowEnergy_fit->GetXaxis()->SetNdivisions(205);

            hHighLowEnergy_fit->GetYaxis()->SetTitleFont(FONT);
            hHighLowEnergy_fit->GetYaxis()->SetTitleSize(TITLE_SIZE_S);
            hHighLowEnergy_fit->GetYaxis()->SetTitleOffset(TITLE_OFFSET_S_XY);
            hHighLowEnergy_fit->GetYaxis()->SetLabelFont(FONT);
            hHighLowEnergy_fit->GetYaxis()->SetLabelSize(LABEL_SIZE_S);
            hHighLowEnergy_fit->GetYaxis()->SetTitle("High Energy Electron [MeV]");
            hHighLowEnergy_fit->GetYaxis()->SetNdivisions(205);

            hHighLowEnergy_fit->GetZaxis()->SetTitleFont(FONT);
            hHighLowEnergy_fit->GetZaxis()->SetTitleSize(TITLE_SIZE_S);
            hHighLowEnergy_fit->GetZaxis()->SetTitleOffset(TITLE_OFFSET_S_Z_MEDIUM);
            hHighLowEnergy_fit->GetZaxis()->SetLabelFont(FONT);
            hHighLowEnergy_fit->GetZaxis()->SetLabelSize(LABEL_SIZE_S);
            hHighLowEnergy_fit->GetZaxis()->SetTitle("Weight");

            TString st;
            st = "h(#xi_{31} = " + xi_31_str + ")";
            TLatex tl;
            tl.SetNDC();
            tl.SetTextAlign(33);
            tl.SetTextFont(LATEX_FONT);
            tl.SetTextSize(LATEX_FONT_SIZE);
            tl.DrawLatex(LATEX_X_RIGHT, LATEX_Y, st);
            tl.DrawLatex(0.25, LATEX_Y, "(b)");

            pad3_2->cd();
            hWeight_fit->Scale(1.0 / hWeight_fit->Integral());
            hWeight_fit->Draw("hist");
            hWeight_fit->SetTitle("");

            hWeight_fit->GetXaxis()->SetTitle("Weight");
            hWeight_fit->GetXaxis()->SetTitleOffset(7.0);
            hWeight_fit->GetYaxis()->SetTitle("");
            hWeight_fit->GetYaxis()->SetTitleOffset(10.0);
            tl.SetTextAlign(13);
            tl.DrawLatex(0.1, 0.95, "Distribution (Norm)");
            tl.SetTextAlign(33);

            hWeight_fit->GetXaxis()->SetTitleFont(FONT);
            hWeight_fit->GetXaxis()->SetTitleSize(0.8 * TITLE_SIZE_S);
            hWeight_fit->GetXaxis()->SetLabelFont(FONT);
            hWeight_fit->GetXaxis()->SetLabelSize(0.8 * TITLE_SIZE_S);
            hWeight_fit->GetXaxis()->SetNdivisions(402);
            hWeight_fit->GetXaxis()->SetTickLength(0.1);
            
            hWeight_fit->GetYaxis()->SetTitleFont(FONT);
            hWeight_fit->GetYaxis()->SetTitleSize(0.8 * TITLE_SIZE_S);
            hWeight_fit->GetYaxis()->SetLabelFont(FONT);
            hWeight_fit->GetYaxis()->SetLabelSize(0.8 * TITLE_SIZE_S);
            hWeight_fit->GetYaxis()->SetNdivisions(402);
            hWeight_fit->GetYaxis()->SetTickLength(0.1);
        }
        else if(order == 3)
        {
            // mid left
            pad2->cd();
            hHighLowEnergy_diff->SetTitle("");
            hHighLowEnergy_diff->SetContour(10000);
            hHighLowEnergy_diff->Draw("colz");
            hHighLowEnergy_diff->GetXaxis()->SetTitleFont(FONT);
            hHighLowEnergy_diff->GetXaxis()->SetTitleSize(TITLE_SIZE_S);
            hHighLowEnergy_diff->GetXaxis()->SetTitleOffset(TITLE_OFFSET_S_XY);
            hHighLowEnergy_diff->GetXaxis()->SetLabelFont(FONT);
            hHighLowEnergy_diff->GetXaxis()->SetLabelSize(LABEL_SIZE_S);
            hHighLowEnergy_diff->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
            hHighLowEnergy_diff->GetXaxis()->SetNdivisions(205);

            hHighLowEnergy_diff->GetYaxis()->SetTitleFont(FONT);
            hHighLowEnergy_diff->GetYaxis()->SetTitleSize(TITLE_SIZE_S);
            hHighLowEnergy_diff->GetYaxis()->SetTitleOffset(TITLE_OFFSET_S_XY);
            hHighLowEnergy_diff->GetYaxis()->SetLabelFont(FONT);
            hHighLowEnergy_diff->GetYaxis()->SetLabelSize(LABEL_SIZE_S);
            hHighLowEnergy_diff->GetYaxis()->SetTitle("High Energy Electron [MeV]");
            hHighLowEnergy_diff->GetYaxis()->SetNdivisions(205);

            hHighLowEnergy_diff->GetZaxis()->SetTitleFont(FONT);
            hHighLowEnergy_diff->GetZaxis()->SetTitleSize(TITLE_SIZE_S);
            hHighLowEnergy_diff->GetZaxis()->SetTitleOffset(4.0);
            hHighLowEnergy_diff->GetZaxis()->SetLabelFont(FONT);
            hHighLowEnergy_diff->GetZaxis()->SetLabelSize(LABEL_SIZE_S);
            hHighLowEnergy_diff->GetZaxis()->SetTitle("Weight");

            //hHighLowEnergy_diff->SetMaximum(-1.0e-2);

            TString st;
            //st = "h(#xi_{31} = " + xi_31_str + ") - h(#xi_{31} = 0.296)";
            st = "h(#xi_{31} = " + xi_31_str + ") - h(#xi_{31} = 0.0)";
            TLatex tl;
            tl.SetNDC();
            tl.SetTextAlign(33);
            tl.SetTextFont(LATEX_FONT);
            tl.SetTextSize(LATEX_FONT_SIZE);
            tl.DrawLatex(LATEX_X_LEFT, LATEX_Y, st);
            tl.DrawLatex(0.25, LATEX_Y, "(c)");
        }
        else if(order == 4)
        {
            // mid right
            pad4->cd();
            hHighLowEnergy_div->SetTitle("");
            hHighLowEnergy_div->SetContour(10000);
            hHighLowEnergy_div->Draw("colz");
            hHighLowEnergy_div->GetXaxis()->SetTitleFont(FONT);
            hHighLowEnergy_div->GetXaxis()->SetTitleSize(TITLE_SIZE_S);
            hHighLowEnergy_div->GetXaxis()->SetTitleOffset(TITLE_OFFSET_S_XY);
            hHighLowEnergy_div->GetXaxis()->SetLabelFont(FONT);
            hHighLowEnergy_div->GetXaxis()->SetLabelSize(LABEL_SIZE_S);
            hHighLowEnergy_div->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
            hHighLowEnergy_div->GetXaxis()->SetNdivisions(205);

            hHighLowEnergy_div->GetYaxis()->SetTitleFont(FONT);
            hHighLowEnergy_div->GetYaxis()->SetTitleSize(TITLE_SIZE_S);
            hHighLowEnergy_div->GetYaxis()->SetTitleOffset(TITLE_OFFSET_S_XY);
            hHighLowEnergy_div->GetYaxis()->SetLabelFont(FONT);
            hHighLowEnergy_div->GetYaxis()->SetLabelSize(LABEL_SIZE_S);
            hHighLowEnergy_div->GetYaxis()->SetTitle("High Energy Electron [MeV]");
            hHighLowEnergy_div->GetYaxis()->SetNdivisions(205);

            hHighLowEnergy_div->GetZaxis()->SetTitleFont(FONT);
            hHighLowEnergy_div->GetZaxis()->SetTitleSize(TITLE_SIZE_S);
            hHighLowEnergy_div->GetZaxis()->SetTitleOffset(TITLE_OFFSET_S_Z_MEDIUM);
            hHighLowEnergy_div->GetZaxis()->SetLabelFont(FONT);
            hHighLowEnergy_div->GetZaxis()->SetLabelSize(LABEL_SIZE_S);
            hHighLowEnergy_div->GetZaxis()->SetTitle("Weight");
            
            hHighLowEnergy_div->GetZaxis()->SetRangeUser(hHighLowEnergy_div_min, hHighLowEnergy_div_max);

            TString st;
            //st = "#frac{h(#xi_{31} = " + xi_31_str + ")}{h(#xi_{31} = 0.296)}";
            st = "#frac{h(#xi_{31} = " + xi_31_str + ")}{h(#xi_{31} = 0.0)}";
            TLatex tl;
            tl.SetNDC();
            tl.SetTextAlign(33);
            tl.SetTextFont(LATEX_FONT);
            tl.SetTextSize(LATEX_FONT_SIZE);
            tl.DrawLatex(LATEX_X_RIGHT, LATEX_Y, st);
            tl.DrawLatex(0.25, LATEX_Y, "(d)");
        }
        else if(order == 5)
        {
            // bottom
            pad5->cd();
            hHighLowEnergy_sig_2->SetTitle("");
            hHighLowEnergy_sig_2->SetContour(10000);
            hHighLowEnergy_sig_2->Draw("colz");
            hHighLowEnergy_sig_2->GetXaxis()->SetTitleFont(FONT);
            hHighLowEnergy_sig_2->GetXaxis()->SetTitleSize(TITLE_SIZE_L);
            hHighLowEnergy_sig_2->GetXaxis()->SetTitleOffset(2.5);
            hHighLowEnergy_sig_2->GetXaxis()->SetLabelFont(FONT);
            hHighLowEnergy_sig_2->GetXaxis()->SetLabelSize(LABEL_SIZE_L);
            hHighLowEnergy_sig_2->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
            //hHighLowEnergy_sig_2->GetXaxis()->SetNdivisions(205);

            hHighLowEnergy_sig_2->GetYaxis()->SetTitleFont(FONT);
            hHighLowEnergy_sig_2->GetYaxis()->SetTitleSize(TITLE_SIZE_L);
            hHighLowEnergy_sig_2->GetYaxis()->SetTitleOffset(2.5);
            hHighLowEnergy_sig_2->GetYaxis()->SetLabelFont(FONT);
            hHighLowEnergy_sig_2->GetYaxis()->SetLabelSize(LABEL_SIZE_L);
            hHighLowEnergy_sig_2->GetYaxis()->SetTitle("High Energy Electron [MeV]");
            //hHighLowEnergy_sig_2->GetYaxis()->SetNdivisions(205);

            hHighLowEnergy_sig_2->GetZaxis()->SetTitleFont(FONT);
            hHighLowEnergy_sig_2->GetZaxis()->SetTitleSize(TITLE_SIZE_L);
            hHighLowEnergy_sig_2->GetZaxis()->SetTitleOffset(2.5);
            hHighLowEnergy_sig_2->GetZaxis()->SetLabelFont(FONT);
            hHighLowEnergy_sig_2->GetZaxis()->SetLabelSize(LABEL_SIZE_L);
            hHighLowEnergy_sig_2->GetZaxis()->SetTitle("Significance #sigma");
            
            //std::cout << "hHighLowEnergy_sig_2_min=" << hHighLowEnergy_sig_2_min << std::endl;
            //hHighLowEnergy_sig_2->GetZaxis()->SetRangeUser(hHighLowEnergy_sig_2_min * 0.99, hHighLowEnergy_sig_2_max);
            //hHighLowEnergy_sig_2->GetZaxis()->SetRangeUser(0.0, hHighLowEnergy_sig_2_max);

            TString st;
            //st = "Significance #sigma = #frac{h(#xi_{31} = " + xi_31_str + ") - h(#xi_{31} = 0.296)}{e(#xi_{31} = 0.296)}";
            st = "Significance #sigma = #frac{h(#xi_{31} = " + xi_31_str + ") - h(#xi_{31} = 0.0)}{e(#xi_{31} = 0.0)}";
            TLatex tl;
            tl.SetNDC();
            tl.SetTextAlign(33);
            tl.SetTextFont(LATEX_FONT);
            tl.SetTextSize(1.1 * LATEX_FONT_SIZE);
            //tl.SetTextSize(0.03);
            tl.DrawLatex(0.8, 0.32, st);
            tl.SetTextSize(1.5 * LATEX_FONT_SIZE);
            tl.DrawLatex(0.18, 0.9, "(e)");

            double chi2 = 0.0;
            //double chi2_2 = 0.0;
            int nbins = 0;
            for(Int_t ix{1}; ix <= hHighLowEnergy_sig_2->GetNbinsX(); ++ ix)
            {
                for(Int_t iy{1}; iy <= hHighLowEnergy_sig_2->GetNbinsY(); ++ iy)
                {
                    double content = hHighLowEnergy_sig_2->GetBinContent(ix, iy);
                    if(content != -std::numeric_limits<double>::infinity())
                    {
                        //chi2 += std::pow(content, 2.0);
                        chi2 += content;
                        //chi2_2 += std::pow(content, 2.0);
                        ++ nbins;

                        //std::cout << "ix=" << ix << " iy=" << iy << " content=" << content << std::endl;
                    }
                }
            }
            TString chi2str;
            //TString chi2str_2;
            //chi2str.Form("#chi^{2}=#frac{%3.1f}{%d}", chi2, nbins);
            chi2str.Form("#chi^{2} = #frac{%d}{%d} = %3.1f", (int)chi2, nbins, chi2 / (double)nbins);
            //chi2str.Form("#chi^{2}_2 = #frac{%d}{%d} = %3.1f", (int)chi2_2, nbins, chi2_2 / (double)nbins);
            tl.SetTextSize(1.1 * LATEX_FONT_SIZE);
            tl.DrawLatex(0.8, 0.2, chi2str);

            if(0)
            {
                pad5_2->cd();
                
                hHighLowEnergy_sig->SetTitle("");
                hHighLowEnergy_sig->SetContour(10000);
                hHighLowEnergy_sig->Draw("colz");
                hHighLowEnergy_sig->GetXaxis()->SetTitleFont(FONT);
                hHighLowEnergy_sig->GetXaxis()->SetTitleSize(TITLE_SIZE_S);
                hHighLowEnergy_sig->GetXaxis()->SetTitleOffset(TITLE_OFFSET_S_XY);
                hHighLowEnergy_sig->GetXaxis()->SetLabelFont(FONT);
                hHighLowEnergy_sig->GetXaxis()->SetLabelSize(LABEL_SIZE_S);
                hHighLowEnergy_sig->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
                hHighLowEnergy_sig->GetXaxis()->SetNdivisions(205);

                hHighLowEnergy_sig->GetYaxis()->SetTitleFont(FONT);
                hHighLowEnergy_sig->GetYaxis()->SetTitleSize(TITLE_SIZE_S);
                hHighLowEnergy_sig->GetYaxis()->SetTitleOffset(TITLE_OFFSET_S_XY);
                hHighLowEnergy_sig->GetYaxis()->SetLabelFont(FONT);
                hHighLowEnergy_sig->GetYaxis()->SetLabelSize(LABEL_SIZE_S);
                hHighLowEnergy_sig->GetYaxis()->SetTitle("High Energy Electron [MeV]");
                hHighLowEnergy_sig->GetYaxis()->SetNdivisions(205);

                hHighLowEnergy_sig->GetZaxis()->SetTitleFont(FONT);
                hHighLowEnergy_sig->GetZaxis()->SetTitleSize(TITLE_SIZE_S);
                hHighLowEnergy_sig->GetZaxis()->SetTitleOffset(TITLE_OFFSET_S_Z_MEDIUM);
                hHighLowEnergy_sig->GetZaxis()->SetLabelFont(FONT);
                hHighLowEnergy_sig->GetZaxis()->SetLabelSize(LABEL_SIZE_S);
                hHighLowEnergy_sig->GetZaxis()->SetTitle("Significance #sigma"); 
                //std::cout << "hHighLowEnergy_sig_min=" << hHighLowEnergy_sig_min << std::endl;
                //hHighLowEnergy_sig->GetZaxis()->SetRangeUser(hHighLowEnergy_sig_min * 0.99, hHighLowEnergy_sig_max);
                //hHighLowEnergy_sig->GetZaxis()->SetRangeUser(0.0, hHighLowEnergy_sig_max);

/*
                TString st;
                st = "Significance #sigma = #frac{h(#xi_{31} = " + xi_31_str + ") - h(#xi_{31} = 0.296)}{e(#xi_{31} = 0.296)}";
                TLatex tl;
                tl.SetNDC();
                tl.SetTextAlign(33);
                tl.SetTextFont(LATEX_FONT);
                tl.SetTextSize(1.2 * LATEX_FONT_SIZE);
                //tl.SetTextSize(0.03);
                tl.DrawLatex(0.8, 0.9, st);
                tl.SetTextSize(1.5 * LATEX_FONT_SIZE);
                tl.DrawLatex(0.18, 0.9, "(e)");

                double chi2 = 0.0;
                //double chi2_2 = 0.0;
                int nbins = 0;
                for(Int_t ix{1}; ix < hHighLowEnergy_sig->GetNbinsX(); ++ ix)
                {
                    for(Int_t iy{1}; iy < hHighLowEnergy_sig->GetNbinsY(); ++ iy)
                    {
                        double content = hHighLowEnergy_sig->GetBinContent(ix, iy);
                        if(content != -std::numeric_limits<double>::infinity())
                        {
                            chi2 += std::pow(content, 2.0);
                            //chi2 += content;
                            //chi2_2 += std::pow(content, 2.0);
                            ++ nbins;
                        }
                    }
                }
                TString chi2str;
                //TString chi2str_2;
                //chi2str.Form("#chi^{2}=#frac{%3.1f}{%d}", chi2, nbins);
                chi2str.Form("#chi^{2} = #frac{%d}{%d} = %3.1f", (int)chi2, nbins, chi2 / (double)nbins);
                //chi2str.Form("#chi^{2}_2 = #frac{%d}{%d} = %3.1f", (int)chi2_2, nbins, chi2_2 / (double)nbins);
                tl.DrawLatex(0.8, 0.6, chi2str);
  */                  
            }
            else
            {
                pad5_2->cd();
                
                hHighLowEnergy_err->SetTitle("");
                hHighLowEnergy_err->SetContour(10000);
                hHighLowEnergy_err->Draw("colz");
                hHighLowEnergy_err->GetXaxis()->SetTitleFont(FONT);
                hHighLowEnergy_err->GetXaxis()->SetTitleSize(TITLE_SIZE_S);
                hHighLowEnergy_err->GetXaxis()->SetTitleOffset(TITLE_OFFSET_S_XY);
                hHighLowEnergy_err->GetXaxis()->SetLabelFont(FONT);
                hHighLowEnergy_err->GetXaxis()->SetLabelSize(LABEL_SIZE_S);
                hHighLowEnergy_err->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
                hHighLowEnergy_err->GetXaxis()->SetNdivisions(205);

                hHighLowEnergy_err->GetYaxis()->SetTitleFont(FONT);
                hHighLowEnergy_err->GetYaxis()->SetTitleSize(TITLE_SIZE_S);
                hHighLowEnergy_err->GetYaxis()->SetTitleOffset(TITLE_OFFSET_S_XY);
                hHighLowEnergy_err->GetYaxis()->SetLabelFont(FONT);
                hHighLowEnergy_err->GetYaxis()->SetLabelSize(LABEL_SIZE_S);
                hHighLowEnergy_err->GetYaxis()->SetTitle("High Energy Electron [MeV]");
                hHighLowEnergy_err->GetYaxis()->SetNdivisions(205);

                hHighLowEnergy_err->GetZaxis()->SetTitleFont(FONT);
                hHighLowEnergy_err->GetZaxis()->SetTitleSize(TITLE_SIZE_S);
                hHighLowEnergy_err->GetZaxis()->SetTitleOffset(TITLE_OFFSET_S_Z_SMALL);
                hHighLowEnergy_err->GetZaxis()->SetLabelFont(FONT);
                hHighLowEnergy_err->GetZaxis()->SetLabelSize(LABEL_SIZE_S);
                hHighLowEnergy_err->GetZaxis()->SetTitle("Error"); 

                TString st;
                //st = "e(#xi_{31} = 0.296)";
                st = "e(#xi_{31} = 0.0) (HSD)";
                TLatex tl;
                tl.SetNDC();
                tl.SetTextAlign(33);
                tl.SetTextFont(LATEX_FONT);
                tl.SetTextSize(LATEX_FONT_SIZE);
                //tl.SetTextSize(0.03);
                tl.DrawLatex(LATEX_X_RIGHT, LATEX_Y, st);
                tl.SetTextSize(1.5 * LATEX_FONT_SIZE);
                tl.DrawLatex(0.25, LATEX_Y, "(f)");

            }

        }
    }

    cHighLowEnergy_diff->SaveAs("./cHighLowEnergy_diff/cHighLowEnergy_diff_" + xi_31_str_plus + ".png");
    cHighLowEnergy_diff->SaveAs("./cHighLowEnergy_diff/cHighLowEnergy_diff_" + xi_31_str_plus + ".eps");
    cHighLowEnergy_diff->SaveAs("./cHighLowEnergy_diff/cHighLowEnergy_diff_" + xi_31_str_plus + ".pdf");
    cHighLowEnergy_diff->SaveAs("./cHighLowEnergy_diff/cHighLowEnergy_diff_" + xi_31_str_plus + ".C");
#endif
}



#endif // NEWLOGLIKFITTER_DRAW_OUTPUTDIFF_H
