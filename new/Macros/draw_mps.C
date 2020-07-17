TCanvas *c;

void save();

void draw_mps()
{

    TFile *f = new TFile("h_mps_1_0_2020-07-14_highlow.root");
    TH2F* h = (TH2F*)f->Get("h_mps_1_0");
    c = new TCanvas("c_mps_1_0");
    c->SetTicks(2, 2);
    c->SetRightMargin(0.15);
    c->SetBottomMargin(0.15);
    c->SetLogz();
    gStyle->SetPalette(kLightTemperature);
    h->SetTitle("");
    Double_t MIN_Z = 1.0e-02;
    Double_t MAX_Z = 1.0e+03;
    Double_t min = h->GetBinContent(1, 1);
    Double_t min_x, min_y;
    Int_t min_ix, min_iy;
    for(Int_t ix = 1; ix <= h->GetNbinsX(); ++ ix)
    {
        for(Int_t iy = 1; iy <= h->GetNbinsY(); ++ iy)
        {
            Double_t content = h->GetBinContent(ix, iy);
            if(content < MIN_Z)
            {
                h->SetBinContent(ix, iy, MIN_Z);
            }
            if(content < min)
            {
                min = content;
                min_x = h->GetXaxis()->GetBinCenter(ix);
                min_y = h->GetYaxis()->GetBinCenter(iy);
                min_ix = ix;
                min_iy = iy;
            }
        }
    }
//    h->GetZaxis()->SetRangeUser(MIN_Z, MAX_Z);
    h->GetZaxis()->SetLabelOffset(0.005);
    h->GetXaxis()->SetLabelSize(17.0);
    h->GetXaxis()->SetLabelFont(63);
    h->GetYaxis()->SetLabelSize(17.0);
    h->GetYaxis()->SetLabelFont(63);
    h->GetZaxis()->SetLabelSize(17.0);
    h->GetZaxis()->SetLabelFont(63);
    h->GetXaxis()->SetTitleSize(18.0);
    h->GetXaxis()->SetTitleFont(43);
    h->GetYaxis()->SetTitleSize(18.0);
    h->GetYaxis()->SetTitleFont(43);
    h->GetYaxis()->SetTitle("^{150}Nd Amplitude Scale Factor");
    h->GetXaxis()->SetTitle("#xi^{2#nu#beta#beta}_{31}");
    h->GetXaxis()->SetTitleOffset(1.5);
    h->GetXaxis()->SetLabelOffset(0.01);
    h->GetYaxis()->SetLabelOffset(0.01);
    //std::cout << h->GetXaxis()->GetLabelOffset() << std::endl;
    //h->DrawClone("colz");
    TH2F* h2 = (TH2F*)h->Clone();
    h->Draw("colz");


    //Double_t min = 0.0;
    // 1 sigma, 90 %, 99 %
    Double_t clevels[3] = {min + 2.30, min + 4.61, min + 9.21};
    h2->SetLineColor(kBlack);
    h2->SetContour(3, clevels);
    //h2->DrawClone("cont2same");


//    TPaletteAxis *palette = new TPaletteAxis(10.0, 0.2, 20.0, 15.0, h);
//    h->GetListOfFunctions()->Add(palette, "br");


    //gPad->Update();
    c->Update();
    TPaletteAxis *palette = (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette");
    //((TPave*)palette)->SetX1NDC(0.7);
    //((TPave*)palette)->SetX2NDC(0.8);
    palette->SetX1NDC(0.88 + 0.02);
    palette->SetX2NDC(0.92 + 0.02);
    palette->SetY1NDC(0.15);
    palette->SetY2NDC(0.9);
    palette->Draw();
    gPad->Modified();
    gPad->Update();
    c->Modified();

    //c->Modified();
    //c->Update();


    Double_t param_2_min = h->GetYaxis()->GetBinLowEdge(h->GetYaxis()->GetFirst()); //GetMinimum();
    Double_t param_2_max = h->GetYaxis()->GetBinUpEdge(h->GetYaxis()->GetLast()); //GetMaximum();
    Double_t param_1_min = h->GetXaxis()->GetBinLowEdge(h->GetXaxis()->GetFirst()); //GetMaximum();
    Double_t param_1_max = h->GetXaxis()->GetBinUpEdge(h->GetXaxis()->GetLast()); //GetMaximum();
    TLine *lineHSD = new TLine(0.0, param_2_min, 0.0, param_2_max);
    TLine *lineSSD = new TLine(0.296, param_2_min, 0.296, param_2_max);
    TLine *lineY = new TLine(param_1_min, 1.0, param_1_max, 1.0);
    TLine *lineXc = new TLine(param_1_min, min_y, param_1_max, min_y);
    TLine *lineYc = new TLine(min_x, param_2_min, min_x, param_2_max);
    //lineHSD->SetLineColor(kWhite);
    //lineSSD->SetLineColor(kWhite);
    //lineY->SetLineColor(kWhite);
    lineHSD->SetLineColorAlpha(kWhite, 0.5);
    lineSSD->SetLineColorAlpha(kWhite, 0.5);
    lineY->SetLineColorAlpha(kWhite, 0.5);
    lineXc->SetLineColorAlpha(kBlack, 0.5);
    lineYc->SetLineColorAlpha(kBlack, 0.5);
    lineHSD->Draw();
    lineSSD->Draw();
    lineY->Draw();
    Int_t ix_0 = h->GetXaxis()->FindBin(0.0);
    Int_t iy_1 = h->GetYaxis()->FindBin(1.0);
    if(min_ix != ix_0 && min_iy != iy_1)
    {
        lineXc->Draw();
        lineYc->Draw();
    }
    //TMarker *bestfitpoint = new TMarker(min_x, min_y, 106);
    //bestfitpoint->SetMarkerColorAlpha(kBlack, 0.5);
    //bestfitpoint->SetMarkerSize(2.0);
    //bestfitpoint->Draw();
    h2->Draw("cont2same");

}

void save()
{

    c->SaveAs("draw_mps_c_2020-07-13_highlow.pdf");

}
