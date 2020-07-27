

#include <string>
#include <iostream>
#include <fstream>

void scriptfinalize()
{
    
    const int n_param_xy = 301;
    int n_param_1 = n_param_xy;
    int n_param_2 = n_param_xy;
    int n_param_max = n_param_1 * n_param_2;
    //int c_param = 0;

    double param_1_min;
    double param_1_max;

    // param 1 is gA
    // custom range
    param_1_min = 0.1;
    param_1_max = 1.7;
    // fake data values
    //if(mode_fake_data)
    //{
    //    param_1_min = -0.4;
    //    param_1_max = 0.6;
    //}

    double param_2_min;
    double param_2_max;
    
    // param 2 is 150Nd amplitude
    // custom range
    param_2_min = 0.8; //1.1; //0.0; //0.0;
    param_2_max = 2.6; //2.6; //1.8; //2.0; //2.0; //4.0;
    // fake data values
    //if(mode_fake_data)
    //{
    //    param_2_min = 0.2;
    //    param_2_max = 1.8;
    //}


    TString h_mps_name_base;
    h_mps_name_base = "h_mps";
    TString h_mps_name = h_mps_name_base + "_" + to_string(1) + "_" + to_string(0);
    TH2D *h_mps = new TH2D(h_mps_name, h_mps_name,
                           n_param_1, param_1_min, param_1_max,
                           n_param_2, param_2_min, param_2_max);

    h_mps->SetContour(1000);

    TString param_1_name_str = TString("g_{A} Axial Vector #xi_{31}");
    TString param_2_name_str = TString("{150}Nd 2#nu#beta#beta");

    h_mps->GetXaxis()->SetTitle(param_1_name_str);
    h_mps->GetYaxis()->SetTitle(param_2_name_str);

    double min = std::numeric_limits<double>::infinity();
    double min_x = -1.0; //-0.085;
    double min_y = -1.0; //0.87;
    
    double min_before = std::numeric_limits<double>::infinity();
    double min_x_before = -1.0; //-0.085;
    double min_y_before = -1.0; //0.87;

    double fval_min = 0.0;
    fval_min = std::numeric_limits<double>::infinity();


    for(int i = 0; i < 301; ++ i)
    {
        std::string name = std::string("script") + std::to_string(i) + ".txt";
        //std::ofstream ofs(name);
        std::cout << "NUMBER=" << std::to_string(i) << std::endl;
        std::string output_name = std::string("mps_resultsmatrix") + "_"
                                + "after" + "_"
                                + "JID"
                                + to_string(i) + ".txt";
        //ofs << "START_INDEX=" << std::to_string(i) << std::endl;
        //ofs << "STOP_INDEX=" << std::to_string(i + 1) << std::endl;
        //ofs << "RUNNING=false" << std::endl;
        

        std::ifstream ifs(output_name);
        if(!ifs.is_open())
        {
            std::cout << "Error failed to open file " << output_name << std::endl;
            continue;
        }


        
        double min_stripe = std::numeric_limits<double>::infinity();
        double min_stripe_y = 0.0;

        while(!ifs.eof() && ifs.good())
        {

            int n_1, n_2;
            double t_param_1, t_param_2, fval;
            std::vector<double> params;
            std::vector<double> param_errs;

            std::string line;
            std::getline(ifs, line);
            //std::stringstream ss;
            //ss << s;
            std::size_t token_index = 0;
            for(;;)
            {
                std::string token;
                std::string remaining_string;
                std::size_t pos = line.find(" ");
                if(pos == std::string::npos)
                {
                    break;
                }
                token = line.substr(0, pos);
                remaining_string = line.substr(pos + 1);
                //std::cout << "token=" << token << "." << std::endl;
                //std::cout << "remaining=" << remaining_string << "." << std::endl;
                line = remaining_string;

                if(token_index == 0)
                {
                    n_1 = std::stoi(token);
                    std::cout << "n_1=" << n_1;
                }
                else if(token_index == 1)
                {
                    n_2 = std::stoi(token);
                    std::cout << " n_2=" << n_2;
                }
                else if(token_index == 2)
                {
                    t_param_1 = std::stod(token);
                    std::cout << " t_param_1=" << t_param_1;
                }
                else if(token_index == 3)
                {
                    t_param_2 = std::stod(token);
                    std::cout << " t_param_2=" << t_param_2;
                }
                else if(token_index == 4)
                {
                    fval = std::stod(token);
                    std::cout << " fval=" << fval;
                }
                else
                {
                    std::size_t param_ix = token_index - 5;
                    if(param_ix % 2 == 0)
                    {
                        params.push_back(std::stod(token));
                    }
                    else
                    {
                        param_errs.push_back(std::stod(token));
                    }
                }

                ++ token_index;
            }
            std::cout << " params.size()=" << params.size()
                      << " param_errs.size()=" << param_errs.size() << std::endl;

            std::cout << "fval=" << fval << std::endl;
            h_mps->SetBinContent(n_1, n_2, fval);

            if(fval < min)
            {
                min = fval;
                min_x = t_param_1;
                min_y = t_param_2;
            }

            // TODO: don't think this currently works need to reset on
            // each stripe
            if(fval < min_stripe)
            {
                min_stripe = fval;
                min_stripe_y = t_param_2;
            }

        }
        std::cout << "EOF" << std::endl;
        // while

        ifs.close();

    }
    // for i 0 .. 301

    
    ///////////////////////////////////////////////////////////////////
    // c_mps
    ///////////////////////////////////////////////////////////////////

    int param_1_ix_external = 1;
    int param_2_ix_external = 0;

    TString param_1_ix_str_external;
    param_1_ix_str_external.Form("%i", param_1_ix_external);
    TString param_2_ix_str_external;
    param_2_ix_str_external.Form("%i", param_2_ix_external);

    TString c_mps_name_base = "c_mps";
    TString c_mps_name = c_mps_name_base + "_" + param_1_ix_str_external + "_" + param_2_ix_str_external;

    TCanvas *c_mps = new TCanvas(c_mps_name, c_mps_name);
    c_mps->SetTicks(2, 2);
    c_mps->SetRightMargin(0.15);
    c_mps->SetBottomMargin(0.15);
    c_mps->SetLogz();
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
    h_mps->SetStats(0);
    h_mps->SetTitle("");
    h_mps->GetZaxis()->SetLabelOffset(0.005);
    h_mps->GetXaxis()->SetLabelSize(17.0);
    h_mps->GetXaxis()->SetLabelFont(63);
    h_mps->GetYaxis()->SetLabelSize(17.0);
    h_mps->GetYaxis()->SetLabelFont(63);
    h_mps->GetZaxis()->SetLabelSize(17.0);
    h_mps->GetZaxis()->SetLabelFont(63);
    h_mps->GetXaxis()->SetTitleSize(18.0);
    h_mps->GetXaxis()->SetTitleFont(43);
    h_mps->GetYaxis()->SetTitleSize(18.0);
    h_mps->GetYaxis()->SetTitleFont(43);
    h_mps->GetYaxis()->SetTitle("^{150}Nd Amplitude Scale Factor");
    h_mps->GetXaxis()->SetTitle("#xi^{2#nu#beta#beta}_{31}");
    h_mps->GetXaxis()->SetTitleOffset(1.5);
    h_mps->GetXaxis()->SetLabelOffset(0.01);
    h_mps->GetYaxis()->SetLabelOffset(0.01);
    TH2D *h_mps_contour = (TH2D*)h_mps->Clone("h_mps_1_0_clone");
    h_mps->Draw("colz");


    std::cout << "min=" << min << " min_x=" << min_x << " min_y=" << min_y << std::endl;
    //double clevels[3] = {min + 1.0, min + 2.0, min + 3.0};
    double clevels[3] = {min + 2.30, min + 4.61, min + 9.21};
    //double clevels[3] = {2.30, 4.61, 9.21}; // true minimum is 0.0 for HSD
    h_mps_contour->SetLineColor(kBlack);
    h_mps_contour->SetContour(3, clevels);

    c_mps->Update();
    TPaletteAxis *palette = (TPaletteAxis*)h_mps->GetListOfFunctions()->FindObject("palette");
    //((TPave*)palette)->SetX1NDC(0.7);
    //((TPave*)palette)->SetX2NDC(0.8);
    palette->SetX1NDC(0.88 + 0.02);
    palette->SetX2NDC(0.92 + 0.02);
    palette->SetY1NDC(0.15);
    palette->SetY2NDC(0.9);
    palette->Draw();
    gPad->Modified();
    gPad->Update();
    c_mps->Modified();
    

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
    Int_t min_ix = h_mps->GetXaxis()->FindBin(min_x);
    Int_t min_iy = h_mps->GetXaxis()->FindBin(min_y);
    Int_t ix_0 = h_mps->GetXaxis()->FindBin(0.0);
    Int_t iy_1 = h_mps->GetXaxis()->FindBin(1.0);
    if(min_ix != ix_0 && min_iy != iy_1)
    {
        lineXc->Draw();
        lineYc->Draw();
    }
    //TMarker *bestfitpoint = new TMarker(min_x, min_y, 106);
    //bestfitpoint->SetMarkerColorAlpha(kBlack, 0.5);
    //bestfitpoint->SetMarkerSize(2.0);
    //bestfitpoint->Draw();

    /*
    std::vector<TLine*> linesteps;
    for(std::size_t ix_walk = 0; ix_walk < ll_walk_save.size() - 1; ++ ix_walk)
    {
        std::pair<double, double> p1 = ll_walk_save.at(ix_walk);
        std::pair<double, double> p2 = ll_walk_save.at(ix_walk + 1);
        Double_t x1 = p1.first;
        Double_t x2 = p2.first;
        Double_t y1 = p1.second;
        Double_t y2 = p2.second;
        std::cout << "ix_walk=" << ix_walk << " " << x1 << " " << y1 << std::endl;
        TLine *linestep = new TLine(x1, y1, x2, y2);
        linestep->SetLineColorAlpha(kRed, 0.1);
        linestep->SetLineWidth(2);
        linestep->Draw();
        linesteps.push_back(linestep);
    }
    */

    h_mps_contour->Draw("cont2same");
    //TString c_fname_png = c_mps_name + datetimestamp_TString + ".png";
    //TString c_fname_pdf = c_mps_name + datetimestamp_TString + ".pdf";
    TString c_fname = c_mps_name + "_scriptfinalize";
    TString c_fname_png = c_fname + ".png";
    TString c_fname_pdf = c_fname + ".pdf";
    c_mps->SaveAs(c_fname_png);
    c_mps->SaveAs(c_fname_pdf);
    //h_mps = nullptr;



}
