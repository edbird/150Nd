#ifndef NEWLOGLIKFITTER_DRAW_MPS_HELPER_H
#define NEWLOGLIKFITTER_DRAW_MPS_HELPER_H




void newloglikfitter_mps_draw_helper_colz
(
    TCanvas *&c_mps,
    mpsdrawdata &mps_draw_data_sysall,
    mpsdrawdata &mps_draw_data_sysnone
)
{

    std::cout << __func__ << std::endl;

    ///////////////////////////////////////////////////////////////////
    // c_mps before/after
    ///////////////////////////////////////////////////////////////////

    if(0 || (MODE_PARALLEL == 0))
    {
        TString c_mps_name_base = "c_mps";
        TString c_mps_name = c_mps_name_base + "_" + "after";

        /*TCanvas **/c_mps = new TCanvas(c_mps_name, c_mps_name);
        c_mps->SetTicks(2, 2);
        c_mps->SetRightMargin(0.15);
        c_mps->SetBottomMargin(0.15);
        //c_mps->SetLogz();
        TVirtualPad *padret = c_mps->cd();
        if(padret == nullptr)
        {
            std::cout << "PAD FAIL" << std::endl;
            std::cin.get();
        }
        //c_mps->GetPad()->cd();
        //c_mps_v.push_back(c_mps);
        //c_mps = nullptr;
        //c_mps->cd();


        //TH2D *h_mps = mps_draw_data_sysnone.h_mps;
        TH2D *h_mps = mps_draw_data_sysall.h_mps;
        
        /*
        TH2D *h_mps_both = new TH2D("h_mps_both", "h_mps_both",
                                    2 * h_mps->GetNbinsX(),
                                    h_mps->GetXaxis()->GetBinLowEdge(1),
                                    h_mps->GetXaxis()->GetBinUpEdge(h_mps->GetNbinsX()),
                                    h_mps->GetNbinsY(),
                                    h_mps->GetYaxis()->GetBinLowEdge(1),
                                    h_mps->GetYaxis()->GetBinUpEdge(h_mps->GetNbinsY()));
        h_mps_both->SetStats(0);
        h_mps_both->SetTitle(0);
        h_mps_both->SetContour(1000);
        c_mps->SetLogz();
        h_mps_both->SetTitle("");
        h_mps_both->SetStats(0);
        h_mps_both->GetZaxis()->SetLabelOffset(0.005);
        h_mps_both->GetXaxis()->SetLabelSize(17.0);
        h_mps_both->GetXaxis()->SetLabelFont(43);
        h_mps_both->GetYaxis()->SetLabelSize(17.0);
        h_mps_both->GetYaxis()->SetLabelFont(43);
        h_mps_both->GetZaxis()->SetLabelSize(17.0);
        h_mps_both->GetZaxis()->SetLabelFont(43);
        h_mps_both->GetXaxis()->SetTitleSize(18.0);
        h_mps_both->GetXaxis()->SetTitleFont(43);
        h_mps_both->GetYaxis()->SetTitleSize(18.0);
        h_mps_both->GetYaxis()->SetTitleFont(43);
        h_mps_both->GetYaxis()->SetTitle("^{150}Nd Amplitude Scale Factor");
        h_mps_both->GetXaxis()->SetTitle("#xi^{2#nu#beta#beta}_{31}");
        h_mps_both->GetXaxis()->SetTitleOffset(1.5);
        h_mps_both->GetYaxis()->SetTitleOffset(1.2);
        h_mps_both->GetXaxis()->SetLabelOffset(0.01);
        h_mps_both->GetYaxis()->SetLabelOffset(0.01);

        for(Int_t j = 1; j <= h_mps->GetNbinsY(); ++ j)
        {
            for(Int_t i = 1; i <= h_mps->GetNbinsX(); ++ i)
            {
                TH2D *h_mps_sysall = mps_draw_data_sysall.h_mps;
                Double_t content_sysnone = h_mps->GetBinContent(i, j);
                Double_t content_sysall = h_mps_sysall->GetBinContent(i, j);
                h_mps_both->SetBinContent(2 * (i - 1) + 0 + 1, j, content_sysnone);
                h_mps_both->SetBinContent(2 * (i - 1) + 1 + 1, j, content_sysall);
            }
        }
        h_mps_both->Draw("colz");
        */

        double min = mps_draw_data_sysnone.min;
        double min_x = mps_draw_data_sysnone.min_x;
        double min_y = mps_draw_data_sysnone.min_y;
        double min_point[2] = {mps_draw_data_sysnone.min_point[0],
                               mps_draw_data_sysnone.min_point[1]};
        double param_1_min = mps_draw_data_sysnone.param_1_min;
        double param_1_max = mps_draw_data_sysnone.param_1_max;
        double param_2_min = mps_draw_data_sysnone.param_2_min;
        double param_2_max = mps_draw_data_sysnone.param_2_max;

        double min_fval = mps_draw_data_sysnone.min_fval;
        std::cout << "min_fval=" << min_fval << std::endl;

        h_mps->SetTitle("");
        h_mps->SetStats(0);
        h_mps->GetZaxis()->SetLabelOffset(0.005);
        h_mps->GetXaxis()->SetLabelSize(17.0);
        h_mps->GetXaxis()->SetLabelFont(43);
        h_mps->GetYaxis()->SetLabelSize(17.0);
        h_mps->GetYaxis()->SetLabelFont(43);
        h_mps->GetZaxis()->SetLabelSize(17.0);
        h_mps->GetZaxis()->SetLabelFont(43);
        h_mps->GetXaxis()->SetTitleSize(18.0);
        h_mps->GetXaxis()->SetTitleFont(43);
        h_mps->GetYaxis()->SetTitleSize(18.0);
        h_mps->GetYaxis()->SetTitleFont(43);
        h_mps->GetYaxis()->SetTitle("^{150}Nd Amplitude Scale Factor                    ");
        h_mps->GetXaxis()->SetTitle("#xi^{2#nu#beta#beta}_{31}");
        h_mps->GetXaxis()->SetTitleOffset(1.5);
        h_mps->GetYaxis()->SetTitleOffset(1.2);
        h_mps->GetXaxis()->SetLabelOffset(0.01);
        h_mps->GetYaxis()->SetLabelOffset(0.01);
        TH2D *h_mps_contour = (TH2D*)h_mps->Clone("h_mps_1_0_clone");
        h_mps->Draw("colz");
        //h_mps->Draw("AXIS");


        std::cout << "min=" << min << " min_x=" << min_x << " min_y=" << min_y << std::endl;
        //double clevels[3] = {min + 1.0, min + 2.0, min + 3.0};
        double clevels[3] = {min_fval + 2.30, min_fval + 4.61, min_fval + 9.21};
        //double clevels[3] = {2.30, 4.61, 9.21}; // true minimum is 0.0 for HSD
        //h_mps_contour->SetLineColor(kGreen); //kRed
        h_mps_contour->SetLineColor(kCyan); //kRed
        //h_mps_contour->SetLineColor(kGray + 2); //kRed
        h_mps_contour->SetLineWidth(2);
//        h_mps_contour->SetLineStyle(0); // in drawmps.C but not used
        h_mps_contour->SetContour(3, clevels);

        c_mps->Update();
        TPaletteAxis *palette = (TPaletteAxis*)h_mps->GetListOfFunctions()->FindObject("palette");
        //TPaletteAxis *palette = (TPaletteAxis*)h_mps_both->GetListOfFunctions()->FindObject("palette");
        palette->SetX1NDC(0.88 + 0.03);
        palette->SetX2NDC(0.92 + 0.03);
        palette->SetY1NDC(0.15);
        palette->SetY2NDC(0.9);
        palette->Draw();
        gPad->Modified();
        gPad->Update();
        c_mps->Modified();

        TLine *lineHSD = new TLine(0.0, param_2_min, 0.0, param_2_max);
        TLine *lineSSD = new TLine(0.296, param_2_min, 0.296, param_2_max);
        TLine *lineY = new TLine(param_1_min, 1.0, param_1_max, 1.0);
        //TLine *lineXc = new TLine(param_1_min, min_y, param_1_max, min_y);
        //TLine *lineYc = new TLine(min_x, param_2_min, min_x, param_2_max);
        TLine *lineXc = new TLine(param_1_min, min_point[1], param_1_max, min_point[1]);
        TLine *lineYc = new TLine(min_point[0], param_2_min, min_point[0], param_2_max);
        //lineHSD->SetLineColor(kWhite);
        //lineSSD->SetLineColor(kWhite);
        //lineY->SetLineColor(kWhite);
        lineHSD->SetLineColorAlpha(kBlack, 1.0);
        lineSSD->SetLineColorAlpha(kBlack, 1.0);
        lineSSD->SetLineStyle(1);

        lineY->SetLineColorAlpha(kBlack, 1.1);
        lineXc->SetLineColorAlpha(kMagenta, 1.0);
        lineYc->SetLineColorAlpha(kMagenta, 1.0);
        Int_t min_ix = h_mps->GetXaxis()->FindBin(min_x);
        Int_t min_iy = h_mps->GetXaxis()->FindBin(min_y);
        Int_t ix_0 = h_mps->GetXaxis()->FindBin(0.0);
        Int_t iy_1 = h_mps->GetXaxis()->FindBin(1.0);
        if(mps_draw_data_sysnone.mode_fake_data_flag == false)
        {
            lineXc->Draw();
            lineYc->Draw();
        }
        if(std::abs(0.0 - min_point[0]) > 1.0e-2)
        {
            lineHSD->Draw();
        }
        if(std::abs(0.296 - min_point[0]) > 1.0e-2)
        {
            lineSSD->Draw();
        }
        if(std::abs(1.0 - min_point[1]) > 1.0e-2)
        {
            lineY->Draw();
        }

        // draw sysnone contour
        h_mps_contour->Draw("cont3same");

        // draw sysnone contour
        TH2D *h_mps_sysnone = mps_draw_data_sysnone.h_mps;
        TH2D *h_mps_contour_sysnone = (TH2D*)h_mps_sysnone->Clone("h_mps_1_0_clone_sysnone");
        h_mps_contour_sysnone->SetLineColor(kGreen);
        //h_mps_contour_sysnone->SetLineColor(kGray);
        h_mps_contour_sysnone->SetLineWidth(2);
        h_mps_contour_sysnone->SetContour(3, clevels);
        h_mps_contour_sysnone->Draw("cont3same");


        ///////////////////////////////////////////////////////////////////////
        // PLOT MARKERS FOR SYSTEMATICS BEST FIT POINTS
        ///////////////////////////////////////////////////////////////////////

        for(int i = 0; i < N_SYSTEMATICS; ++ i)
        {
            if(ENABLE_MIN_POINT_SYSn[i] == true)
            {
                std::cout << "MIN_POINT: i=" << i << " enabled" << std::endl;
                if((mps_draw_data_sysall.min_point_sysn_l[i][0] != 0.0) &&
                   (mps_draw_data_sysall.min_point_sysn_l[i][1] != 0.0))
                {
                    std::cout << "draw! (l) "
                              << mps_draw_data_sysall.min_point_sysn_l[i][0] << " "
                              << mps_draw_data_sysall.min_point_sysn_l[i][1]
                              << std::endl;
//                    mps_draw_data_sysall.mark_min_point_sysn_l[i]->SetMarkerColor(kRed);
                    mps_draw_data_sysall.mark_min_point_sysn_l[i]->Draw();
                    if(false)
                    {
                        mps_draw_data_sysall.line_min_point_sysn_l[i]->Draw();
                    }
                }

                if((mps_draw_data_sysall.min_point_sysn_h[i][0] != 0.0) &&
                   (mps_draw_data_sysall.min_point_sysn_h[i][1] != 0.0))
                {
                    std::cout << "draw! (h) "
                              << mps_draw_data_sysall.min_point_sysn_h[i][0] << " "
                              << mps_draw_data_sysall.min_point_sysn_h[i][1]
                              << std::endl;
//                    mps_draw_data_sysall.mark_min_point_sysn_h[i]->SetMarkerColor(kRed);
                    mps_draw_data_sysall.mark_min_point_sysn_h[i]->Draw();
                    if(false)
                    {
                        mps_draw_data_sysall.line_min_point_sysn_h[i]->Draw();
                    }
                }
                std::cout << std::endl;
            }
            else
            {
                std::cout << "MIN_POINT: i=" << i << " disabled" << std::endl;
                std::cout << std::endl;
            }
        }

/*
        if(ll_walk_save.size() > 0)
        {
            std::vector<TLine*> linesteps;
            for(std::size_t ix_walk = 0; ix_walk < ll_walk_save.size() - 1; ++ ix_walk)
            {
                std::pair<double, double> p1 = ll_walk_save.at(ix_walk);
                std::pair<double, double> p2 = ll_walk_save.at(ix_walk + 1);
                Double_t x1 = p1.first;
                Double_t x2 = p2.first;
                Double_t y1 = p1.second;
                Double_t y2 = p2.second;
                //std::cout << "ix_walk=" << ix_walk << " " << x1 << " " << y1 << std::endl;
                TLine *linestep = new TLine(x1, y1, x2, y2);
                linestep->SetLineColorAlpha(kRed, 0.1);
                linestep->SetLineWidth(2);
                linestep->Draw();
                linesteps.push_back(linestep);
            }
        }
*/

        TString texts[N_SYSTEMATICS] = 
        {
            "ignore 0.1 MeV",
            "Energy Calibration 1.2 %",
            "Efficiency",
            "ignore enrichment",
            "Energy Calibration 3 keV",
            "Foil Thickness virtual",
            "Energy Loss dE/dX virtual",
            "ignore Brem virtual",
            "Foil Thickness nominal",
            "Energy Loss dE/dX nominal",
            "ignore Brem nominal"
        };
        TLegend *leg = new TLegend(0.03, 0.72, 0.31, 0.96);
        for(int i = 0; i < N_SYSTEMATICS; ++ i)
        {
            if(ENABLE_MIN_POINT_SYSn[i] == true)
            {
                TString text;
                text.Form("%i", i);
                leg->AddEntry(mps_draw_data_sysall.mark_min_point_sysn_h[i],
                              texts[i], "P");
            }
        }
        leg->SetMargin(0.12);
        leg->SetBorderSize(5);
        leg->SetShadowColor(kGray + 2);
        leg->SetTextFont(63);
        leg->SetTextSize(15);
        leg->Draw("BR");


    }
}







void newloglikfitter_mps_draw_helper
(
    TCanvas *&c_mps,
    mpsdrawdata &mps_draw_data_sysall,
    mpsdrawdata &mps_draw_data_sysnone
)
{

    std::cout << __func__ << std::endl;

    ///////////////////////////////////////////////////////////////////
    // c_mps before/after
    ///////////////////////////////////////////////////////////////////

    if(0 || (MODE_PARALLEL == 0))
    {
        TString c_mps_name_base = "c_mps";
        TString c_mps_name = c_mps_name_base + "_" + "after";

        /*TCanvas **/c_mps = new TCanvas(c_mps_name, c_mps_name);
        c_mps->SetTicks(2, 2);
        //c_mps->SetRightMargin(0.15);
        c_mps->SetBottomMargin(0.15);
        //c_mps->SetLogz();
        TVirtualPad *padret = c_mps->cd();
        if(padret == nullptr)
        {
            std::cout << "PAD FAIL" << std::endl;
            std::cin.get();
        }
        //c_mps->GetPad()->cd();
        //c_mps_v.push_back(c_mps);
        //c_mps = nullptr;
        //c_mps->cd();


        TH2D *h_mps = mps_draw_data_sysnone.h_mps;
        
        /*
        TH2D *h_mps_both = new TH2D("h_mps_both", "h_mps_both",
                                    2 * h_mps->GetNbinsX(),
                                    h_mps->GetXaxis()->GetBinLowEdge(1),
                                    h_mps->GetXaxis()->GetBinUpEdge(h_mps->GetNbinsX()),
                                    h_mps->GetNbinsY(),
                                    h_mps->GetYaxis()->GetBinLowEdge(1),
                                    h_mps->GetYaxis()->GetBinUpEdge(h_mps->GetNbinsY()));
        h_mps_both->SetStats(0);
        h_mps_both->SetTitle(0);
        h_mps_both->SetContour(1000);
        c_mps->SetLogz();
        h_mps_both->SetTitle("");
        h_mps_both->SetStats(0);
        h_mps_both->GetZaxis()->SetLabelOffset(0.005);
        h_mps_both->GetXaxis()->SetLabelSize(17.0);
        h_mps_both->GetXaxis()->SetLabelFont(43);
        h_mps_both->GetYaxis()->SetLabelSize(17.0);
        h_mps_both->GetYaxis()->SetLabelFont(43);
        h_mps_both->GetZaxis()->SetLabelSize(17.0);
        h_mps_both->GetZaxis()->SetLabelFont(43);
        h_mps_both->GetXaxis()->SetTitleSize(18.0);
        h_mps_both->GetXaxis()->SetTitleFont(43);
        h_mps_both->GetYaxis()->SetTitleSize(18.0);
        h_mps_both->GetYaxis()->SetTitleFont(43);
        h_mps_both->GetYaxis()->SetTitle("^{150}Nd Amplitude Scale Factor");
        h_mps_both->GetXaxis()->SetTitle("#xi^{2#nu#beta#beta}_{31}");
        h_mps_both->GetXaxis()->SetTitleOffset(1.5);
        h_mps_both->GetYaxis()->SetTitleOffset(1.2);
        h_mps_both->GetXaxis()->SetLabelOffset(0.01);
        h_mps_both->GetYaxis()->SetLabelOffset(0.01);

        for(Int_t j = 1; j <= h_mps->GetNbinsY(); ++ j)
        {
            for(Int_t i = 1; i <= h_mps->GetNbinsX(); ++ i)
            {
                TH2D *h_mps_sysall = mps_draw_data_sysall.h_mps;
                Double_t content_sysnone = h_mps->GetBinContent(i, j);
                Double_t content_sysall = h_mps_sysall->GetBinContent(i, j);
                h_mps_both->SetBinContent(2 * (i - 1) + 0 + 1, j, content_sysnone);
                h_mps_both->SetBinContent(2 * (i - 1) + 1 + 1, j, content_sysall);
            }
        }
        h_mps_both->Draw("colz");
        */

        double min = mps_draw_data_sysnone.min;
        double min_x = mps_draw_data_sysnone.min_x;
        double min_y = mps_draw_data_sysnone.min_y;
        double min_point[2] = {mps_draw_data_sysnone.min_point[0],
                               mps_draw_data_sysnone.min_point[1]};
        double param_1_min = mps_draw_data_sysnone.param_1_min;
        double param_1_max = mps_draw_data_sysnone.param_1_max;
        double param_2_min = mps_draw_data_sysnone.param_2_min;
        double param_2_max = mps_draw_data_sysnone.param_2_max;

        double min_fval = mps_draw_data_sysnone.min_fval;
        std::cout << "min_fval=" << min_fval << std::endl;

        h_mps->SetTitle("");
        h_mps->SetStats(0);
        h_mps->GetZaxis()->SetLabelOffset(0.005);
        h_mps->GetXaxis()->SetLabelSize(17.0);
        h_mps->GetXaxis()->SetLabelFont(43);
        h_mps->GetYaxis()->SetLabelSize(17.0);
        h_mps->GetYaxis()->SetLabelFont(43);
        h_mps->GetZaxis()->SetLabelSize(17.0);
        h_mps->GetZaxis()->SetLabelFont(43);
        h_mps->GetXaxis()->SetTitleSize(18.0);
        h_mps->GetXaxis()->SetTitleFont(43);
        h_mps->GetYaxis()->SetTitleSize(18.0);
        h_mps->GetYaxis()->SetTitleFont(43);
        h_mps->GetYaxis()->SetTitle("^{150}Nd Amplitude Scale Factor                    ");
        h_mps->GetXaxis()->SetTitle("#xi^{2#nu#beta#beta}_{31}");
        h_mps->GetXaxis()->SetTitleOffset(1.5);
        h_mps->GetYaxis()->SetTitleOffset(1.2);
        h_mps->GetXaxis()->SetLabelOffset(0.01);
        h_mps->GetYaxis()->SetLabelOffset(0.01);
        TH2D *h_mps_contour = (TH2D*)h_mps->Clone("h_mps_1_0_clone");
        //h_mps->Draw("colz");
        h_mps->Draw("AXIS");


        std::cout << "min=" << min << " min_x=" << min_x << " min_y=" << min_y << std::endl;
        //double clevels[3] = {min + 1.0, min + 2.0, min + 3.0};
        double clevels[3] = {min_fval + 2.30, min_fval + 4.61, min_fval + 9.21};
        //double clevels[3] = {2.30, 4.61, 9.21}; // true minimum is 0.0 for HSD
        //h_mps_contour->SetLineColor(kGreen); //kRed
        //h_mps_contour->SetLineColor(kCyan); //kRed
        h_mps_contour->SetLineColor(kGray + 2); //kRed
        h_mps_contour->SetLineWidth(2);
//        h_mps_contour->SetLineStyle(0); // in drawmps.C but not used
        h_mps_contour->SetContour(3, clevels);

        c_mps->Update();
        TPaletteAxis *palette = (TPaletteAxis*)h_mps->GetListOfFunctions()->FindObject("palette");
        //TPaletteAxis *palette = (TPaletteAxis*)h_mps_both->GetListOfFunctions()->FindObject("palette");
/*        palette->SetX1NDC(0.88 + 0.03);
        palette->SetX2NDC(0.92 + 0.03);
        palette->SetY1NDC(0.15);
        palette->SetY2NDC(0.9);
        palette->Draw();*/
        gPad->Modified();
        gPad->Update();
        c_mps->Modified();

        TLine *lineHSD = new TLine(0.0, param_2_min, 0.0, param_2_max);
        TLine *lineSSD = new TLine(0.296, param_2_min, 0.296, param_2_max);
        TLine *lineY = new TLine(param_1_min, 1.0, param_1_max, 1.0);
        //TLine *lineXc = new TLine(param_1_min, min_y, param_1_max, min_y);
        //TLine *lineYc = new TLine(min_x, param_2_min, min_x, param_2_max);
        TLine *lineXc = new TLine(param_1_min, min_point[1], param_1_max, min_point[1]);
        TLine *lineYc = new TLine(min_point[0], param_2_min, min_point[0], param_2_max);
        //lineHSD->SetLineColor(kWhite);
        //lineSSD->SetLineColor(kWhite);
        //lineY->SetLineColor(kWhite);
        lineHSD->SetLineColorAlpha(kBlack, 1.0);
        lineSSD->SetLineColorAlpha(kBlack, 1.0);
        lineSSD->SetLineStyle(1);

        lineY->SetLineColorAlpha(kBlack, 1.1);
        lineXc->SetLineColorAlpha(kMagenta, 1.0);
        lineYc->SetLineColorAlpha(kMagenta, 1.0);
        Int_t min_ix = h_mps->GetXaxis()->FindBin(min_x);
        Int_t min_iy = h_mps->GetXaxis()->FindBin(min_y);
        Int_t ix_0 = h_mps->GetXaxis()->FindBin(0.0);
        Int_t iy_1 = h_mps->GetXaxis()->FindBin(1.0);
        if(mps_draw_data_sysnone.mode_fake_data_flag == false)
        {
            lineXc->Draw();
            lineYc->Draw();
        }
        if(std::abs(0.0 - min_point[0]) > 1.0e-2)
        {
            lineHSD->Draw();
        }
        if(std::abs(0.296 - min_point[0]) > 1.0e-2)
        {
            lineSSD->Draw();
        }
        if(std::abs(1.0 - min_point[1]) > 1.0e-2)
        {
            lineY->Draw();
        }

        // draw sysnone contour
        h_mps_contour->Draw("cont3same");

        // draw sysall contour
        TH2D *h_mps_sysall = mps_draw_data_sysall.h_mps;
        TH2D *h_mps_contour_sysall = (TH2D*)h_mps_sysall->Clone("h_mps_1_0_clone_sysall");
        //h_mps_contour_sysall->SetLineColor(kGreen);
        h_mps_contour_sysall->SetLineColor(kGray);
        h_mps_contour_sysall->SetLineWidth(2);
        h_mps_contour_sysall->SetContour(3, clevels);
        h_mps_contour_sysall->Draw("cont3same");


        ///////////////////////////////////////////////////////////////////////
        // PLOT MARKERS FOR SYSTEMATICS BEST FIT POINTS
        ///////////////////////////////////////////////////////////////////////

        for(int i = 0; i < N_SYSTEMATICS; ++ i)
        {
            if(ENABLE_MIN_POINT_SYSn[i] == true)
            {
                std::cout << "MIN_POINT: i=" << i << " enabled" << std::endl;
                if((mps_draw_data_sysall.min_point_sysn_l[i][0] != 0.0) &&
                   (mps_draw_data_sysall.min_point_sysn_l[i][1] != 0.0))
                {
                    std::cout << "draw! (l) "
                              << mps_draw_data_sysall.min_point_sysn_l[i][0] << " "
                              << mps_draw_data_sysall.min_point_sysn_l[i][1]
                              << std::endl;
//                    mps_draw_data_sysall.mark_min_point_sysn_l[i]->SetMarkerColor(kRed);
                    mps_draw_data_sysall.mark_min_point_sysn_l[i]->Draw();
                    if(false)
                    {
                        mps_draw_data_sysall.line_min_point_sysn_l[i]->Draw();
                    }
                }

                if((mps_draw_data_sysall.min_point_sysn_h[i][0] != 0.0) &&
                   (mps_draw_data_sysall.min_point_sysn_h[i][1] != 0.0))
                {
                    std::cout << "draw! (h) "
                              << mps_draw_data_sysall.min_point_sysn_h[i][0] << " "
                              << mps_draw_data_sysall.min_point_sysn_h[i][1]
                              << std::endl;
//                    mps_draw_data_sysall.mark_min_point_sysn_h[i]->SetMarkerColor(kRed);
                    mps_draw_data_sysall.mark_min_point_sysn_h[i]->Draw();
                    if(false)
                    {
                        mps_draw_data_sysall.line_min_point_sysn_h[i]->Draw();
                    }
                }
                std::cout << std::endl;
            }
            else
            {
                std::cout << "MIN_POINT: i=" << i << " disabled" << std::endl;
                std::cout << std::endl;
            }
        }

/*
        if(ll_walk_save.size() > 0)
        {
            std::vector<TLine*> linesteps;
            for(std::size_t ix_walk = 0; ix_walk < ll_walk_save.size() - 1; ++ ix_walk)
            {
                std::pair<double, double> p1 = ll_walk_save.at(ix_walk);
                std::pair<double, double> p2 = ll_walk_save.at(ix_walk + 1);
                Double_t x1 = p1.first;
                Double_t x2 = p2.first;
                Double_t y1 = p1.second;
                Double_t y2 = p2.second;
                //std::cout << "ix_walk=" << ix_walk << " " << x1 << " " << y1 << std::endl;
                TLine *linestep = new TLine(x1, y1, x2, y2);
                linestep->SetLineColorAlpha(kRed, 0.1);
                linestep->SetLineWidth(2);
                linestep->Draw();
                linesteps.push_back(linestep);
            }
        }
*/

        TString texts[N_SYSTEMATICS] = 
        {
            "ignore 0.1 MeV",
            "Energy Calibration 1.2 %",
            "Efficiency",
            "ignore enrichment",
            "Energy Calibration 3 keV",
            "Foil Thickness virtual",
            "Energy Loss dE/dX virtual",
            "ignore Brem virtual",
            "Foil Thickness nominal",
            "Energy Loss dE/dX nominal",
            "ignore Brem nominal"
        };
        TLegend *leg = new TLegend(0.03, 0.72, 0.31, 0.96);
        for(int i = 0; i < N_SYSTEMATICS; ++ i)
        {
            if(ENABLE_MIN_POINT_SYSn[i] == true)
            {
                TString text;
                text.Form("%i", i);
                leg->AddEntry(mps_draw_data_sysall.mark_min_point_sysn_h[i],
                              texts[i], "P");
            }
        }
        leg->SetMargin(0.12);
        leg->SetBorderSize(5);
        leg->SetShadowColor(kGray + 2);
        leg->SetTextFont(63);
        leg->SetTextSize(15);
        leg->Draw("BR");


    }
}


#endif // NEWLOGLIKFITTER_DRAW_MPS_HELPER_H
