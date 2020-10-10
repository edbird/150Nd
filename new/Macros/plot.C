
std::size_t LD_REWEIGHT_DATA_2(
    TH2D *&h_out,
    const std::string &filename,
    const std::string& h_name,
    const std::string& name//,
    //const double low,
    //const double high
    )
{

    

    //TTree *t_G0 = new TTree();
    //Long64_t count_G0 = t_G0->ReadFile("../../data/150Nd-data/dG150Nd/G0/dG0.dat", "electronEnergy1:electronEnergy2:electronWeight");
    //std::cout << "count_G0=" << count_G0 << std::endl;

    //Long64_t count_G0 = 0;
    std::vector<std::vector<double>> data;
    std::size_t read_count = 0;
    std::ifstream if_Gn(filename);
    double electronEnergy1, electronEnergy2, electronWeight;
    while(if_Gn.good())
    {
        if_Gn >> electronEnergy1 >> electronEnergy2 >> electronWeight;
        //++ count_G0;
        ++ read_count;

       
        //if(read_count < 10)
        //{
        //    std::cout << "electronEnergy1=" << electronEnergy1
        //          << " electronEnergy2=" << electronEnergy2
        //          << " electronWeight=" << electronWeight << std::endl;
        //}
        //std::cin.get();


        data.emplace_back(std::vector<double>());
        data.back().emplace_back(electronEnergy1);
        data.back().emplace_back(electronEnergy2);
        data.back().emplace_back(electronWeight);
    }
    std::cout << "read_count=" << read_count << std::endl;

    //TTree *t_G2 = new TTree();
    //Long64_t count_G2 = t_G2->ReadFile("../../data/150Nd-data/dG150Nd/G2/dG2.dat", "electronEnergy1:electronEnergy2:electronWeight");
    //std::cout << "count_G2=" << count_G2 << std::endl;

    //output = out;
    //outputG0 = t_G0;
    //outputG2 = t_G2;
    //countG0 = count_G0;
    //countG2 = count_G2;

    // solve quadratic
    //const double A = 0.5;
    //const double B = 0.5;
    //double C = -(double)count;
    //double dimension_xy = (-B + std::sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
    Long64_t count = read_count;
    double dimension_xy = (int)(-0.5 + std::sqrt(0.25 - 2.0 * (-(double)count)));
    if(dimension_xy != 3371)
    {
        std::cout << "error! dimension_xy != 3371, dimension_xy=" << dimension_xy << std::endl;
    }
    //const double dimension_xy{3371};

    std::cout << "dimension_xy=" << dimension_xy << std::endl;


    ////////////////////////////////////////////////////////////////////////////
    // CONVERT INPUT DATA TO HISTOGRAM FORMAT
    // Note: Added 2018-04-23 (After INTERMEDIATE DATA below)
    // Note: These histograms do NOT have the phase space variable included
    ////////////////////////////////////////////////////////////////////////////
                                      
    //TH2D *h_nEqNull = new TH2D("h_nEqNull", "nEqNull", dimension_xy, 0.0, bb_Q, dimension_xy, 0.0, bb_Q);
    //TH2D *h_nEqTwo = new TH2D("h_nEqTwo", "nEqTwo", dimension_xy, 0.0, bb_Q, dimension_xy, 0.0, bb_Q);


    //const Int_t dimension_xy{3371};
    //const Int_t dimension_xy{1001};
    // don't plot raw data
    //TH2D *h_nEqNull = new TH2D("h_nEqNull", "", dimension_xy, 0.0, bb_Q, dimension_xy, 0.0, bb_Q);
    //TH2D *h_nEqTwo = new TH2D("h_nEqTwo", "", dimension_xy, 0.0, bb_Q, dimension_xy, 0.0, bb_Q);
    //TH2D *h_nEqNull = new TH2D("h_nEqNull", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //TH2D *h_nEqTwo = new TH2D("h_nEqTwo", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //TH2D *h_ratio = new TH2D("h_ratio", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //h_nEqNull = new TH2D("h_nEqNull", "nEqNull", dimension_xy, 0.0, bb_Q, dimension_xy, 0.0, bb_Q);
    //h_nEqTwo = new TH2D("h_nEqTwo", "nEqTwo", dimension_xy, 0.0, bb_Q, dimension_xy, 0.0, bb_Q);
    //TH2D *h_ = new TH2D(hname.c_str(), name.c_str(), dimension_xy, low, high, dimension_xy, low, high);
    const Double_t keV = 1.0e-03;
    const Double_t E_min = 0.38e-03;
    const Double_t limit_lower  = E_min - 0.5 * keV;
    const Double_t limit_higher = E_min - 0.5 * keV + dimension_xy * keV;
    TH2D *h_ = new TH2D(h_name.c_str(), h_name.c_str(),
                        dimension_xy, limit_lower, limit_higher,
                        dimension_xy, limit_lower, limit_higher);
    //h_nEqTwo = new TH2D("h_nEqTwo", "nEqTwo", dimension_xy, 0.0, bb_Q, dimension_xy, 0.0, bb_Q);
    h_->SetStats(0);
    //h_nEqTwo->SetStats(0);
    //h_ratio->SetStats(0);

    // required because Nd-150 data is triangular
    std::size_t rd_ix{0};
    for(std::size_t i{0}; i < dimension_xy; ++ i)
    {
        for(std::size_t j{0}; j < dimension_xy; ++ j)
        {
            if(i < dimension_xy - j)
            {
                //h_nEqNull->SetBinContent(i, j, data_nEqNull.at(rd_ix)[2]);
                //h_nEqTwo->SetBinContent(i, j, data_nEqTwo.at(rd_ix)[2]);
                h_->SetBinContent(i, j, data.at(rd_ix)[2]);
                //if(i < 5 && j < 10)
                //{
                //    std::cout << "fill: " << i << " " << j << " " << data.at(rd_ix)[2] << std::endl;
                //}
            
                ++ rd_ix;
            }
            else
            {
                //h_nEqNull->SetBinContent(i, j, 0.0);
                //h_nEqTwo->SetBinContent(i, j, 0.0);
                h_->SetBinContent(i, j, 0.0);
            }
            //h_nEqNull->SetBinContent(i, j, data_nEqNull.at(i * dimension_xy + j)[2]);
            //h_nEqTwo->SetBinContent(i, j, data_nEqTwo.at(i * dimension_xy + j)[2]);
            //if(i < dimension_xy - j)
            //if(i + j < dimension_xy - 1)
            //{
                // TODO: move above lines to inside this if
                //h_ratio->SetBinContent(i, j, ratio.at(i * dimension_xy + j)[2]);
            //}
        }
    }
    std::cout << "rd_ix_max=" << rd_ix << std::endl;
    std::cout << "Finished constructing input data histogram" << std::endl;
    
    //h_nEqNull_return = h_nEqNull;
    //h_nEqTwo_return = h_nEqTwo;
    h_out = h_;

    //TCanvas *c_ = new TCanvas("c_", "c_");
    //h_->Draw("colz");
    //c_->SaveAs("hout.png");
    //std::cout << "histogram done" << std::endl;
    //std::cin.get();

    return read_count;

}





void plot()
{

    TH2D *h_nEqNull;
    TH2D *h_nEqTwo;
    TH2D *h_nEqTwoTwo;
    TH2D *h_nEqFour;
    double psiN0;
    double psiN2;
    double psiN22;
    double psiN4;

    
    // load data
    std::cout << "attempting to load spectral data from file" << std::endl;

    //bb_Q = 3.368;
    //bb_Q = 3.37138;
    // ramdisk, should be faster?
    std::size_t count_G0 =
        LD_REWEIGHT_DATA_2(h_nEqNull,
                           //"/home/ecb/150Nd/150Nd-data/dG150Nd/G0/dG0.dat",
                           "/home/ecb/100Mo-150Nd/gA_theoretical_files/psf-nuclei/150Nd/0-N0/nEqNull.dat",
                           "h_nEqNull",
                           "nEqNull"//,
                           //0.0, bb_Q
                           );

    std::size_t count_G2 =
        LD_REWEIGHT_DATA_2(h_nEqTwo,
                           //"/home/ecb/150Nd/150Nd-data/dG150Nd/G2/dG2.dat",
                           "/home/ecb/100Mo-150Nd/gA_theoretical_files/psf-nuclei/150Nd/1-N2/nEqTwo.dat",
                           "h_nEqTwo",
                           "nEqTwo"//,  0.0, bb_Q
                           );

    std::size_t count_G22 =
        LD_REWEIGHT_DATA_2(h_nEqTwoTwo,
                           //"/home/ecb/150Nd/150Nd-data/dG150Nd/G2/dG2.dat",
                           "/home/ecb/100Mo-150Nd/gA_theoretical_files/psf-nuclei/150Nd/2-N22/nEqTwoTwo.dat",
                           "h_nEqTwoTwo",
                           "nEqTwoTwo"//,  0.0, bb_Q
                           );

    std::size_t count_G4 =
        LD_REWEIGHT_DATA_2(h_nEqFour,
                           //"/home/ecb/150Nd/150Nd-data/dG150Nd/G2/dG2.dat",
                           "/home/ecb/100Mo-150Nd/gA_theoretical_files/psf-nuclei/150Nd/3-N4/nEqFour.dat",
                           "h_nEqFour",
                           "nEqFour"//,  0.0, bb_Q
                           );

    // phase space integrals
    const Double_t G0_ps_integral_MeV = 0.759721E-45; //0.420438E-45;
    const Double_t G0_ps_integral_yrinv = 0.364244E-16; //0.201577E-16;
    const Double_t G2_ps_integral_MeV = 0.429791E-45; //0.744684E-45;
    const Double_t G2_ps_integral_yrinv = 0.206061E-16; //0.357034E-16;
    const Double_t G22_ps_integral_MeV = 0.704956E-46;
    const Double_t G22_ps_integral_yrinv = 0.337987E-17;
    const Double_t G4_ps_integral_MeV = 0.304316E-45;
    const Double_t G4_ps_integral_yrinv = 0.145903E-16;

    //xi_31_baseline = 0.0;

    double count = 0;
    if((count_G0 == count_G2) &&
       (count_G0 == count_G22) &&
       (count_G0 == count_G4))
    {
        count = count_G0;
    }
    else
    {
        std::cout << "error: "
                  << "count_G0=" << count_G0 << ", "
                  << "count_G2=" << count_G2 << ", "
                  << "count_G22=" << count_G22 << ", "
                  << "count_G4=" << count_G4 << std::endl;
    }

    psiN0 = G0_ps_integral_MeV;
    psiN2 = G2_ps_integral_MeV; // TODO: check this is the correct option
    psiN22 = G22_ps_integral_MeV;
    psiN4 = G4_ps_integral_MeV;
    
    std::cout << "histogram format constructed" << std::endl;


    const Int_t TITLE_FONT = 43;
    const Int_t TITLE_SIZE = 25; // 20
    const Int_t LABEL_FONT = 43;
    const Int_t LABEL_SIZE = 25; // 20

    const Double_t X_TITLE_OFFSET = 1.3; // 1.5
    const Double_t Y_TITLE_OFFSET = 1.0; // 1.2
    const Double_t Z_TITLE_OFFSET = 1.2; // 1.3

    const Double_t LATEX_X = 0.77;
    const Double_t LATEX_Y = 0.85;

    TCanvas *c_nEqNull;
    c_nEqNull = new TCanvas("c_nEqNull", "c_nEqNull"); //, 4000, 3000);
    TH2D *h_nEqNull_clone = (TH2D*)h_nEqNull->Clone();
    h_nEqNull_clone->SetTitle("");
    h_nEqNull_clone->GetXaxis()->SetTitle("E_{e_{1}} Electron Energy [MeV]");
    h_nEqNull_clone->GetXaxis()->SetTitleFont(TITLE_FONT);
    h_nEqNull_clone->GetXaxis()->SetTitleSize(TITLE_SIZE);
    h_nEqNull_clone->GetXaxis()->SetTitleOffset(X_TITLE_OFFSET);
    h_nEqNull_clone->GetXaxis()->SetLabelFont(LABEL_FONT);
    h_nEqNull_clone->GetXaxis()->SetLabelSize(LABEL_SIZE);
    h_nEqNull_clone->GetYaxis()->SetTitle("E_{e_{2}} Electron Energy [MeV]");
    h_nEqNull_clone->GetYaxis()->SetTitleFont(TITLE_FONT);
    h_nEqNull_clone->GetYaxis()->SetTitleSize(TITLE_SIZE);
    h_nEqNull_clone->GetYaxis()->SetTitleOffset(Y_TITLE_OFFSET);
    h_nEqNull_clone->GetYaxis()->SetLabelFont(LABEL_FONT);
    h_nEqNull_clone->GetYaxis()->SetLabelSize(LABEL_SIZE);
    h_nEqNull_clone->GetZaxis()->SetTitle("Decay Rate");
    h_nEqNull_clone->GetZaxis()->SetTitleFont(TITLE_FONT);
    h_nEqNull_clone->GetZaxis()->SetTitleSize(TITLE_SIZE);
    h_nEqNull_clone->GetZaxis()->SetTitleOffset(Z_TITLE_OFFSET);
    h_nEqNull_clone->GetZaxis()->SetLabelFont(LABEL_FONT);
    h_nEqNull_clone->GetZaxis()->SetLabelSize(LABEL_SIZE);
    c_nEqNull->SetRightMargin(0.18);
    c_nEqNull->SetLeftMargin(0.12);
    c_nEqNull->SetBottomMargin(0.15);
    h_nEqNull_clone->SetContour(1000);
    h_nEqNull_clone->Draw("colz");

    TLatex latexlabel;
    latexlabel.SetNDC();
    latexlabel.SetTextAlign(33);
    latexlabel.SetTextFont(43);
    latexlabel.SetTextSize(30);

    latexlabel.DrawLatex(LATEX_X, LATEX_Y, "#frac{d^{2}G_{0}}{dT_{e_{1}}dT_{e_{2}}}");

    c_nEqNull->SaveAs("c_nEqNull.png");
    c_nEqNull->SaveAs("c_nEqNull.pdf");
    c_nEqNull->SaveAs("c_nEqNull.C");
    //delete c_nEqNull;

    TCanvas *c_nEqTwo;
    c_nEqTwo = new TCanvas("c_nEqTwo", "c_nEqTwo"); //, 4000, 3000);
    TH2D *h_nEqTwo_clone = (TH2D*)h_nEqTwo->Clone();
    h_nEqTwo_clone->SetTitle("");
    h_nEqTwo_clone->GetXaxis()->SetTitle("E_{e_{1}} Electron Energy [MeV]");
    h_nEqTwo_clone->GetXaxis()->SetTitleFont(TITLE_FONT);
    h_nEqTwo_clone->GetXaxis()->SetTitleSize(TITLE_SIZE);
    h_nEqTwo_clone->GetXaxis()->SetTitleOffset(X_TITLE_OFFSET);
    h_nEqTwo_clone->GetXaxis()->SetLabelFont(LABEL_FONT);
    h_nEqTwo_clone->GetXaxis()->SetLabelSize(LABEL_SIZE);
    h_nEqTwo_clone->GetYaxis()->SetTitle("E_{e_{2}} Electron Energy [MeV]");
    h_nEqTwo_clone->GetYaxis()->SetTitleFont(TITLE_FONT);
    h_nEqTwo_clone->GetYaxis()->SetTitleSize(TITLE_SIZE);
    h_nEqTwo_clone->GetYaxis()->SetTitleOffset(Y_TITLE_OFFSET);
    h_nEqTwo_clone->GetYaxis()->SetLabelFont(LABEL_FONT);
    h_nEqTwo_clone->GetYaxis()->SetLabelSize(LABEL_SIZE);
    h_nEqTwo_clone->GetZaxis()->SetTitle("Decay Rate");
    h_nEqTwo_clone->GetZaxis()->SetTitleFont(TITLE_FONT);
    h_nEqTwo_clone->GetZaxis()->SetTitleSize(TITLE_SIZE);
    h_nEqTwo_clone->GetZaxis()->SetTitleOffset(Z_TITLE_OFFSET);
    h_nEqTwo_clone->GetZaxis()->SetLabelFont(LABEL_FONT);
    h_nEqTwo_clone->GetZaxis()->SetLabelSize(LABEL_SIZE);
    c_nEqTwo->SetRightMargin(0.18);
    c_nEqTwo->SetLeftMargin(0.12);
    c_nEqTwo->SetBottomMargin(0.15);
    //h_nEqTwo_clone->SetNumberContours(256);
    h_nEqTwo_clone->SetContour(1000);
    h_nEqTwo_clone->Draw("colz");

    latexlabel.DrawLatex(LATEX_X, LATEX_Y, "#frac{d^{2}G_{2}}{dT_{e_{1}}dT_{e_{2}}}");

    c_nEqTwo->SaveAs("c_nEqTwo.png");
    c_nEqTwo->SaveAs("c_nEqTwo.pdf");
    c_nEqTwo->SaveAs("c_nEqTwo.C");
    //delete c_nEqTwo;

    TCanvas *c_nEqTwoTwo;
    c_nEqTwoTwo = new TCanvas("c_nEqTwoTwo", "c_nEqTwoTwo"); //, 4000, 3000);
    TH2D *h_nEqTwoTwo_clone = (TH2D*)h_nEqTwoTwo->Clone();
    h_nEqTwoTwo_clone->SetTitle("");
    h_nEqTwoTwo_clone->GetXaxis()->SetTitle("E_{e_{1}} Electron Energy [MeV]");
    h_nEqTwoTwo_clone->GetXaxis()->SetTitleFont(TITLE_FONT);
    h_nEqTwoTwo_clone->GetXaxis()->SetTitleSize(TITLE_SIZE);
    h_nEqTwoTwo_clone->GetXaxis()->SetTitleOffset(X_TITLE_OFFSET);
    h_nEqTwoTwo_clone->GetXaxis()->SetLabelFont(LABEL_FONT);
    h_nEqTwoTwo_clone->GetXaxis()->SetLabelSize(LABEL_SIZE);
    h_nEqTwoTwo_clone->GetYaxis()->SetTitle("E_{e_{2}} Electron Energy [MeV]");
    h_nEqTwoTwo_clone->GetYaxis()->SetTitleFont(TITLE_FONT);
    h_nEqTwoTwo_clone->GetYaxis()->SetTitleSize(TITLE_SIZE);
    h_nEqTwoTwo_clone->GetYaxis()->SetTitleOffset(Y_TITLE_OFFSET);
    h_nEqTwoTwo_clone->GetYaxis()->SetLabelFont(LABEL_FONT);
    h_nEqTwoTwo_clone->GetYaxis()->SetLabelSize(LABEL_SIZE);
    h_nEqTwoTwo_clone->GetZaxis()->SetTitle("Decay Rate");
    h_nEqTwoTwo_clone->GetZaxis()->SetTitleFont(TITLE_FONT);
    h_nEqTwoTwo_clone->GetZaxis()->SetTitleSize(TITLE_SIZE);
    h_nEqTwoTwo_clone->GetZaxis()->SetTitleOffset(Z_TITLE_OFFSET);
    h_nEqTwoTwo_clone->GetZaxis()->SetLabelFont(LABEL_FONT);
    h_nEqTwoTwo_clone->GetZaxis()->SetLabelSize(LABEL_SIZE);
    c_nEqTwoTwo->SetRightMargin(0.18);
    c_nEqTwoTwo->SetLeftMargin(0.12);
    c_nEqTwoTwo->SetBottomMargin(0.15);
    //h_nEqTwoTwo_clone->SetNumberContours(256);
    h_nEqTwoTwo_clone->SetContour(1000);
    h_nEqTwoTwo_clone->Draw("colz");

    latexlabel.DrawLatex(LATEX_X, LATEX_Y, "#frac{d^{2}G_{22}}{dT_{e_{1}}dT_{e_{2}}}");

    c_nEqTwoTwo->SaveAs("c_nEqTwoTwo.png");
    c_nEqTwoTwo->SaveAs("c_nEqTwoTwo.pdf");
    c_nEqTwoTwo->SaveAs("c_nEqTwoTwo.C");
    //delete c_nEqTwoTwo;
    
    TCanvas *c_nEqFour;
    c_nEqFour = new TCanvas("c_nEqFour", "c_nEqFour"); //, 4000, 3000);
    TH2D *h_nEqFour_clone = (TH2D*)h_nEqFour->Clone();
    h_nEqFour_clone->SetTitle("");
    h_nEqFour_clone->GetXaxis()->SetTitle("E_{e_{1}} Electron Energy [MeV]");
    h_nEqFour_clone->GetXaxis()->SetTitleFont(TITLE_FONT);
    h_nEqFour_clone->GetXaxis()->SetTitleSize(TITLE_SIZE);
    h_nEqFour_clone->GetXaxis()->SetTitleOffset(X_TITLE_OFFSET);
    h_nEqFour_clone->GetXaxis()->SetLabelFont(LABEL_FONT);
    h_nEqFour_clone->GetXaxis()->SetLabelSize(LABEL_SIZE);
    h_nEqFour_clone->GetYaxis()->SetTitle("E_{e_{2}} Electron Energy [MeV]");
    h_nEqFour_clone->GetYaxis()->SetTitleFont(TITLE_FONT);
    h_nEqFour_clone->GetYaxis()->SetTitleSize(TITLE_SIZE);
    h_nEqFour_clone->GetYaxis()->SetTitleOffset(Y_TITLE_OFFSET);
    h_nEqFour_clone->GetYaxis()->SetLabelFont(LABEL_FONT);
    h_nEqFour_clone->GetYaxis()->SetLabelSize(LABEL_SIZE);
    h_nEqFour_clone->GetZaxis()->SetTitle("Decay Rate");
    h_nEqFour_clone->GetZaxis()->SetTitleFont(TITLE_FONT);
    h_nEqFour_clone->GetZaxis()->SetTitleSize(TITLE_SIZE);
    h_nEqFour_clone->GetZaxis()->SetTitleOffset(Z_TITLE_OFFSET);
    h_nEqFour_clone->GetZaxis()->SetLabelFont(LABEL_FONT);
    h_nEqFour_clone->GetZaxis()->SetLabelSize(LABEL_SIZE);
    c_nEqFour->SetRightMargin(0.18);
    c_nEqFour->SetLeftMargin(0.12);
    c_nEqFour->SetBottomMargin(0.15);
    //h_nEqFour_clone->SetNumberContours(256);
    h_nEqFour_clone->SetContour(1000);
    h_nEqFour_clone->Draw("colz");

    latexlabel.DrawLatex(LATEX_X, LATEX_Y, "#frac{d^{2}G_{4}}{dT_{e_{1}}dT_{e_{2}}}");

    c_nEqFour->SaveAs("c_nEqFour.png");
    c_nEqFour->SaveAs("c_nEqFour.pdf");
    c_nEqFour->SaveAs("c_nEqFour.C");
    //delete c_nEqFour;


}
