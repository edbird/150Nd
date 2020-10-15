







void drawmps()
{

    bool g_mode_fake_data = true;
    const bool V_ENABLE_SYSALL = true;

    const int number_job_id = 0;
    const std::string &output_name = "noparallel";
    const int start_index = 0;
    const int stop_index = 51;


    ///////////////////////////////////////////////////////////////////////////
    // testing - my phase space
    // 2020-08-19: new version, cuts out a lot of the stuff from the code
    // I origionally copied from which didn't make sense in this config
    // plot phase space for Nd150 and Mo100 parameters
    ///////////////////////////////////////////////////////////////////////////


    TString c_mps_name_base = "c_mps_after_drawmps";
    TString c_mps_name = c_mps_name_base;
    
    std::cout << "rendering: " << c_mps_name << std::endl;

    // TODO: NOTE: have to change value of "stop_index" as well
    const int n_param_xy = 31;//301; // 1001
    int n_param_1 = n_param_xy; //300;
    int n_param_2 = n_param_xy; //300;
    int n_param_max = n_param_1 * n_param_2;
    int c_param = 0;

    //double param_1 = AdjustActs[param_1_ix];

    double param_1_min;
    double param_1_max;

    // param 1 is gA
    // custom range
    param_1_min = 0.1; //-0.4; //-0.5; //1.0; //-0.5;
    param_1_max = 1.7; //0.6; //1.6; //0.5; //2.5; //5.0; //2.5;
    // after changing the psiN0, psiN2 values...
    param_1_min = -0.3;
    param_1_max = 1.7; // adjusted for SYS1 energy offset was 1.5
    //param_1_min = -0.4;
    //param_1_max = 1.6; TODO
    // fake data values
    if(g_mode_fake_data == true)
    {
//        param_1_min = -0.4;
//        param_1_max = 0.6;
        param_1_min = -0.7;
        param_1_max = 0.7;
    }
    
    // with systematics
    if(g_mode_fake_data == false)
    {
        param_1_min = -0.7;
        param_1_max = 2.1;
    }
    else if(g_mode_fake_data == true)
    {
        param_1_min = -0.5;
        param_1_max = 0.7;
    }

    // hack to get HSD
    //param_1_min = -0.1;
    //param_1_max = +0.1;

    //double param_2 = AdjustActs[param_2_ix];

    double param_2_min;
    double param_2_max;
    
    // param 2 is 150Nd amplitude
    // custom range
    param_2_min = 0.8; //1.1; //0.0; //0.0;
    param_2_max = 2.6; //2.6; //1.8; //2.0; //2.0; //4.0;
    // after changing the psiN0, psiN2 values...
    param_2_min = 0.95;
    param_2_max = 1.3;
    //param_2_min = 0.0;
    //param_2_max = 3.0;  //TODO
    // fake data values
    if(g_mode_fake_data == true)
    {
//        param_2_min = 0.2;
//        param_2_max = 1.8
        param_2_max = 1.15;
        param_2_min = 0.85;
    }
    
    // with systematics
    if(g_mode_fake_data == false)
    {
        param_2_min = 0.75;
        param_2_max = 1.5;
    }
    else if(g_mode_fake_data == true)
    {
        param_2_min = 0.85;
        param_2_max = 1.15;
    }

    // hack to get HSD
    //param_1_min = 0.9999;
    //param_1_max = 1.0001;


    TString h_mps_name_base;
    if(g_mode_fake_data == true)
    {
        h_mps_name_base = "h_mps_fake_data";
    }
    if(g_mode_fake_data == false)
    {
        h_mps_name_base = "h_mps";
    }
    TString h_mps_name = h_mps_name_base;

    //std::cout << h_mps_name << " param_1=" << param_1 << " sigma_1=" << sigma_1
    //                        << " param_1_min=" << param_1_min << " param_1_max=" << param_1_max
    //                        << " param_2=" << param_2 << " sigma_2=" << sigma_2
    //                        << " param_2_min=" << param_2_min << " param_2_max=" << param_2_max
    //                        << std::endl;

    TH2D *h_mps = new TH2D(h_mps_name, h_mps_name,
                           n_param_1, param_1_min, param_1_max,
                           n_param_2, param_2_min, param_2_max); 
    //h_mps_v.push_back(h_mps);
    //h_mps = nullptr;

    //h_mps->GetZaxis()->SetRangeUser(0.0, 1.0e+04);
    h_mps->SetContour(1000);
    
    //TString param_1_name_str = TString(paramNameMap[param_1_ix_external]);
    //TString param_2_name_str = TString(paramNameMap[param_2_ix_external]);

    //h_mps->GetXaxis()->SetTitle(param_1_name_str);
    //h_mps->GetYaxis()->SetTitle(param_2_name_str);
    h_mps->GetXaxis()->SetTitle("^{150}Nd Amplitude Scale Factor");
    h_mps->GetYaxis()->SetTitle("#xi_{31}^{2#nu#beta#beta}");

    // reset params array
    // now code moved to new function, simply use new variables (local)
    //std::vector<double> params = theParameterState.Params();
    //std::vector<double> param_errs = theParameterState.Errors();



    double min = std::numeric_limits<double>::infinity();
    double min_x = -1.0; //-0.085;
    double min_y = -1.0; //0.87;
    
    double min_before = std::numeric_limits<double>::infinity();
    double min_x_before = -1.0; //-0.085;
    double min_y_before = -1.0; //0.87;



    TString h_mps_name_base_before;
    if(g_mode_fake_data == true)
    {
        h_mps_name_base_before = "h_mps_fake_data_before";
    }
    if(g_mode_fake_data == false)
    {
        h_mps_name_base_before = "h_mps_before";
    }
    TString h_mps_name_before = h_mps_name_base_before;

    TH2D *h_mps_before = new TH2D(h_mps_name_before, h_mps_name_before,
                           n_param_1, param_1_min, param_1_max,
                           n_param_2, param_2_min, param_2_max); 

    h_mps_before->SetContour(1000);
    
    //TString param_1_name_str = TString(paramNameMap[param_1_ix_external]);
    //TString param_2_name_str = TString(paramNameMap[param_2_ix_external]);

    //h_mps_before->GetXaxis()->SetTitle(param_1_name_str);
    //h_mps_before->GetYaxis()->SetTitle(param_2_name_str);
    h_mps_before->GetXaxis()->SetTitle("^{150}Nd Amplitude Scale Factor");
    h_mps_before->GetYaxis()->SetTitle("#xi_{31}^{2#nu#beta#beta}");








    std::string output_name_append;
    if(V_ENABLE_SYSALL == false)
    {
        output_name_append += "_STAT";
    }
    else if(V_ENABLE_SYSALL == true)
    {
        output_name_append += "_STATSYS";
    }
    if(g_mode_fake_data == false)
    {
        output_name_append += "_data";
    }
    else if(g_mode_fake_data == true)
    {
        output_name_append += "_fake";
    }

    std::string ofs_resultsmatrix_before_fname =
        output_name + output_name_append + "_before" + "_"
        + "JID" + std::to_string(number_job_id)
        + ".txt";

    //"mps_resultsmatrix_after"
    std::string ofs_resultsmatrix_after_fname =
        output_name + output_name_append + "_after" + "_"
        + "JID" + std::to_string(number_job_id)
        + ".txt";

    std::ifstream ofs_resultsmatrix_before(ofs_resultsmatrix_before_fname);
    std::ifstream ofs_resultsmatrix_after(ofs_resultsmatrix_after_fname);
    
    if(!ofs_resultsmatrix_before.is_open())
    {
        std::cout << "Error: could not open " << ofs_resultsmatrix_before_fname << std::endl;
        return -1;
    }
    if(!ofs_resultsmatrix_after.is_open())
    {
        std::cout << "Error: could not open " << ofs_resultsmatrix_after_fname << std::endl;
        return -1;
    }

    std::cout << "*****************************************************" << std::endl;
    std::cout << "*****************************************************" << std::endl;
    std::cout << "reading data from " << ofs_resultsmatrix_before_fname << std::endl;
    std::cout << "reading data from " << ofs_resultsmatrix_after_fname << std::endl;
    std::cout << "*****************************************************" << std::endl;
    std::cout << "*****************************************************" << std::endl;


//    double min = 0.0;
//    min = std::numeric_limits<double>::infinity();

    double min_stripe = std::numeric_limits<double>::infinity();
    double min_stripe_y = 0.0;


    ///////////////////////////////////////////////////////////////////////////
    // read "before"
    ///////////////////////////////////////////////////////////////////////////

    std::size_t line_count = 1;
    int n_1_last = -1;
    int n_2_last = -1;
    double t_param_1, t_param_2;
    while(!ofs_resultsmatrix_before.eof())
    {
    //std::cin.get();
        std::stringstream ss;
        std::string s;
        std::getline(ofs_resultsmatrix_before, s);
        ++ line_count;

        ss << s;

        int n_1, n_2;
        /*double t_param_1, t_param_2;*/
        double fval_before;
        std::vector<double> params_before;
        std::vector<double> param_errs_before;

        ss >> n_1 >> n_2;
        ss >> t_param_1 >> t_param_2;
        ss >> fval_before;
        //std::cout << "n_1=" << n_1 << " n_2=" << n_2 << std::endl;
        //std::cout << "t_param_1=" << t_param_1 << " t_param_2=" << t_param_2 << std::endl;
        //std::cout << "fval_before=" << fval_before << std::endl;
        for(;;)
        {
            try
            {
                if(ss.peek() == std::char_traits<wchar_t>::eof())
                {
                    break;
                }
                double tmp1, tmp2;
                ss >> tmp1 >> tmp2;
                params_before.push_back(tmp1);
                param_errs_before.push_back(tmp2);
            }
            catch(...)
            {
                break;
            }
        }
        //std::cout << "line: " << line_count << " params_before.size()=" << params_before.size() << " param_errs_before.size()=" << param_errs_before.size() << std::endl;

        // This detects n_1 changing (loop of outer for loop)
        if(n_1 != n_1_last)
        {
            if(n_1_last != -1)
            {
                std::cout << "min_stripe=" << min_stripe << " min_stripe_x=" << t_param_1 << " min_stripe_y=" << min_stripe_y << std::endl;

                min_stripe = std::numeric_limits<double>::infinity();
                min_stripe_y = 0.0;
            }

            n_1_last = n_1;
        }

        if(n_2 != n_2_last)
        {

            n_2_last = n_2;
        }

        int bin_ix = h_mps_before->GetNbinsX() - n_1;
        //std::cout << "bin_ix=" << bin_ix << std::endl;
        //std::cout << "t_param_1=" << t_param_1 << std::endl;
        //std::cout << "bin center: " << h_mps_before->GetXaxis()->GetBinCenter(bin_ix) << std::endl;

        int bin_iy = h_mps_before->GetNbinsY() - n_2;
        //std::cout << "bin_iy=" << bin_iy << std::endl;
        //std::cout << "t_param_2=" << t_param_2 << std::endl;
        //std::cout << "bin center: " << h_mps_before->GetYaxis()->GetBinCenter(bin_iy) << std::endl;


        if(fval_before < min_stripe)
        {
            min_stripe = fval_before;
            min_stripe_y = t_param_2;
        }

        if(fval_before < min_before)
        {
            min_before = fval_before;
            min_x_before = t_param_1;
            min_y_before = t_param_2;
        }

        std::cout << "ix=" << bin_ix << " iy=" << bin_iy << " fval=" << fval_before << std::endl;
        h_mps_before->SetBinContent(bin_ix, bin_iy, fval_before);

        ++ line_count;
    }
    std::cout << "min_stripe=" << min_stripe << " min_stripe_x=" << t_param_1 << " min_stripe_y=" << min_stripe_y << std::endl;
    std::cout << "read: " << ofs_resultsmatrix_before_fname << " -> done" << std::endl;



    ///////////////////////////////////////////////////////////////////////////
    // read "after"
    ///////////////////////////////////////////////////////////////////////////

    line_count = 1;
    n_1_last = -1;
    n_2_last = -1;
    while(!ofs_resultsmatrix_after.eof())
    {
        
        std::stringstream ss;
        std::string s;
        std::getline(ofs_resultsmatrix_after, s);
        ++ line_count;

        ss << s;

        int n_1, n_2;
        /*double t_param_1, t_param_2;*/
        double fval_after;
        std::vector<double> params_after;
        std::vector<double> param_errs_after;

        ss >> n_1 >> n_2;
        ss >> t_param_1 >> t_param_2;
        ss >> fval_after;
        for(;;)
        {
            try
            {
                if(ss.peek() == std::char_traits<wchar_t>::eof())
                {
                    break;
                }
                double tmp1, tmp2;
                ss >> tmp1 >> tmp2;
                params_after.push_back(tmp1);
                param_errs_after.push_back(tmp2);
            }
            catch(...)
            {
                break;
            }
        }
        std::cout << "line: " << line_count << " params_after.size()=" << params_after.size() << " param_errs_after.size()=" << param_errs_after.size() << std::endl;

        // This detects n_1 changing (loop of outer for loop)
        if(n_1 != n_1_last)
        {
            std::cout << "min_stripe=" << min_stripe << " min_stripe_x=" << t_param_1 << " min_stripe_y=" << min_stripe_y << std::endl;

            min_stripe = std::numeric_limits<double>::infinity();
            min_stripe_y = 0.0;
            n_1_last = n_1;
        }

        if(n_2 != n_2_last)
        {

            n_2_last = n_2;
        }

        int bin_ix = h_mps->GetNbinsX() - n_1;
        std::cout << "bin_ix=" << bin_ix << std::endl;
        std::cout << "t_param_1=" << t_param_1 << std::endl;
        std::cout << "bin center: " << h_mps->GetXaxis()->GetBinCenter(bin_ix) << std::endl;

        int bin_iy = h_mps->GetNbinsY() - n_2;
        std::cout << "bin_iy=" << bin_iy << std::endl;
        std::cout << "t_param_2=" << t_param_2 << std::endl;
        std::cout << "bin center: " << h_mps->GetYaxis()->GetBinCenter(bin_iy) << std::endl;


        if(fval_after < min_stripe)
        {
            min_stripe = fval_after;
            min_stripe_y = t_param_2;
        }
        
        if(fval_after < min)
        {
            min = fval_after;
            min_x = t_param_1;
            min_y = t_param_2;
        }
       
        h_mps->SetBinContent(bin_ix, bin_iy, fval_after);

        ++ line_count;
    }
    std::cout << "read: " << ofs_resultsmatrix_after_fname << " -> done" << std::endl;











    ///////////////////////////////////////////////////////////////////
    // c_mps
    ///////////////////////////////////////////////////////////////////
    //if(0 || (MODE_PARALLEL == 0))
    {
        TCanvas *c_mps = new TCanvas(c_mps_name, c_mps_name);
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
        h_mps->SetTitle("");
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
        h_mps->GetYaxis()->SetTitle("^{150}Nd Amplitude Scale Factor");
        h_mps->GetXaxis()->SetTitle("#xi^{2#nu#beta#beta}_{31}");
        h_mps->GetXaxis()->SetTitleOffset(1.5);
        h_mps->GetYaxis()->SetTitleOffset(1.2);
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

        // TODO:
        /*
        if((min_point_sys1_l[0] != 0.0) &&
           (min_point_sys1_l[1] != 0.0))
        {
        std::cout << "DRAW MARK SYS1L" << std::endl;
            std::cout << "SYS1L: " << min_point_sys1_l[0] << " " << min_point_sys1_l[1] << std::endl;
            TMarker *mark_min_point_sys1_l = new TMarker(min_point_sys1_l[0],
                                                        min_point_sys1_l[1],
                                                        106);
            mark_min_point_sys1_l->SetMarkerColorAlpha(kWhite, 0.0);
            mark_min_point_sys1_l->SetMarkerSize(2.0);
            mark_min_point_sys1_l->Draw();
        }

        if((min_point_sys1_h[0] != 0.0) &&
           (min_point_sys1_h[1] != 0.0))
        {
        std::cout << "DRAW MARK SYS1H" << std::endl;
            std::cout << "SYS1H: " << min_point_sys1_h[0] << " " << min_point_sys1_h[1] << std::endl;
            TMarker *mark_min_point_sys1_h = new TMarker(min_point_sys1_h[0],
                                                        min_point_sys1_h[1],
                                                        106);
            mark_min_point_sys1_h->SetMarkerColorAlpha(kWhite, 0.0);
            mark_min_point_sys1_h->SetMarkerSize(2.0);
            mark_min_point_sys1_h->Draw();
        }
        */

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

        h_mps_contour->Draw("cont2same");
        //TString c_fname_png = c_mps_name + datetimestamp_TString + ".png";
        //TString c_fname_pdf = c_mps_name + datetimestamp_TString + ".pdf";
        TString c_fname = c_mps_name + "_"
                       + "JID" + std::to_string(number_job_id);// + "_"
                       //+ datetimestamp_TString;
        TString c_fname_png = c_fname + ".png";
        TString c_fname_pdf = c_fname + ".pdf";
        std::cout << "*****************************************************" << std::endl;
        std::cout << "c_fname=" << c_fname << std::endl;
        std::cout << "is the filename legal?" << std::endl;
        std::cout << "*****************************************************" << std::endl;
        c_mps->SaveAs(c_fname_png);
        c_mps->SaveAs(c_fname_pdf);
        //h_mps = nullptr;
    }


    ///////////////////////////////////////////////////////////////////
    // c_mps_before
    ///////////////////////////////////////////////////////////////////
    //if(0 || (MODE_PARALLEL == 0))
    {
        TString c_mps_name_base_before = "c_mps_before_drawmps";
        TString c_mps_name_before = c_mps_name_base_before;

        TCanvas *c_mps_before = new TCanvas(c_mps_name_before, c_mps_name_before);
        c_mps_before->SetTicks(2, 2);
        c_mps_before->SetRightMargin(0.15);
        c_mps_before->SetBottomMargin(0.15);
        //c_mps_before->SetLogz();
        //c_mps_before->GetPad()->cd();
        //c_mps_before_v.push_back(c_mps_before);
        //c_mps_before = nullptr;
        //c_mps_before->cd();
        h_mps_before->SetTitle("");
        h_mps_before->GetZaxis()->SetLabelOffset(0.005);
        h_mps_before->GetXaxis()->SetLabelSize(17.0);
        h_mps_before->GetXaxis()->SetLabelFont(43);
        h_mps_before->GetYaxis()->SetLabelSize(17.0);
        h_mps_before->GetYaxis()->SetLabelFont(43);
        h_mps_before->GetZaxis()->SetLabelSize(17.0);
        h_mps_before->GetZaxis()->SetLabelFont(43);
        h_mps_before->GetXaxis()->SetTitleSize(18.0);
        h_mps_before->GetXaxis()->SetTitleFont(43);
        h_mps_before->GetYaxis()->SetTitleSize(18.0);
        h_mps_before->GetYaxis()->SetTitleFont(43);
        h_mps_before->GetYaxis()->SetTitle("^{150}Nd Amplitude Scale Factor");
        h_mps_before->GetXaxis()->SetTitle("#xi^{2#nu#beta#beta}_{31}");
        h_mps_before->GetXaxis()->SetTitleOffset(1.5);
        h_mps_before->GetYaxis()->SetTitleOffset(1.2);
        h_mps_before->GetXaxis()->SetLabelOffset(0.01);
        h_mps_before->GetYaxis()->SetLabelOffset(0.01);
        TH2D *h_mps_contour_before = (TH2D*)h_mps_before->Clone("h_mps_before_1_0_clone");
        h_mps_before->Draw("colz");


        std::cout << "min_before=" << min_before << " min_x_before=" << min_x_before << " min_y_before=" << min_y_before << std::endl;
        //double clevels[3] = {min + 1.0, min + 2.0, min + 3.0};
        double clevels_before[3] = {min_before + 2.30, min_before + 4.61, min_before + 9.21};
        //double clevels[3] = {2.30, 4.61, 9.21}; // true minimum is 0.0 for HSD
        h_mps_contour_before->SetLineColor(kBlack);
        h_mps_contour_before->SetContour(3, clevels_before);

        c_mps_before->Update();
        TPaletteAxis *palette_before = (TPaletteAxis*)h_mps_before->GetListOfFunctions()->FindObject("palette");
        //((TPave*)palette)->SetX1NDC(0.7);
        //((TPave*)palette)->SetX2NDC(0.8);
        palette_before->SetX1NDC(0.88 + 0.03);
        palette_before->SetX2NDC(0.92 + 0.03);
        palette_before->SetY1NDC(0.15);
        palette_before->SetY2NDC(0.9);
        palette_before->Draw();
        gPad->Modified();
        gPad->Update();
        c_mps_before->Modified();
        

        TLine *lineHSD_before = new TLine(0.0, param_2_min, 0.0, param_2_max);
        TLine *lineSSD_before = new TLine(0.296, param_2_min, 0.296, param_2_max);
        TLine *lineY_before = new TLine(param_1_min, 1.0, param_1_max, 1.0);
        TLine *lineXc_before = new TLine(param_1_min, min_y_before, param_1_max, min_y_before);
        TLine *lineYc_before = new TLine(min_x_before, param_2_min, min_x_before, param_2_max);
        //lineHSD->SetLineColor(kWhite);
        //lineSSD->SetLineColor(kWhite);
        //lineY->SetLineColor(kWhite);
        lineHSD_before->SetLineColorAlpha(kWhite, 0.5);
        lineSSD_before->SetLineColorAlpha(kWhite, 0.5);
        lineY_before->SetLineColorAlpha(kWhite, 0.5);
        lineXc_before->SetLineColorAlpha(kBlack, 0.5);
        lineYc_before->SetLineColorAlpha(kBlack, 0.5);
        lineHSD_before->Draw();
        lineSSD_before->Draw();
        lineY_before->Draw();
        Int_t min_ix_before = h_mps_before->GetXaxis()->FindBin(min_x_before);
        Int_t min_iy_before = h_mps_before->GetXaxis()->FindBin(min_y_before);
        Int_t ix_0_before = h_mps_before->GetXaxis()->FindBin(0.0);
        Int_t iy_1_before = h_mps_before->GetXaxis()->FindBin(1.0);
        if(min_ix_before != ix_0_before && min_iy_before != iy_1_before)
        {
            lineXc_before->Draw();
            lineYc_before->Draw();
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

        h_mps_contour_before->Draw("cont2same");
        //TString c_fname_before_png = c_mps_name_before + datetimestamp_TString + ".png";
        //TString c_fname_before_pdf = c_mps_name_before + datetimestamp_TString + ".pdf";
        TString c_fname_before = c_mps_name_before + "_"
                               + "JID" + std::to_string(number_job_id);// + "_"
                               //+ datetimestamp_TString;
        TString c_fname_before_png = c_fname_before + ".png";
        TString c_fname_before_pdf = c_fname_before + ".pdf";
        std::cout << "*****************************************************" << std::endl;
        std::cout << "c_fname_beofre=" << c_fname_before << std::endl;
        std::cout << "is the filename legal?" << std::endl;
        std::cout << "*****************************************************" << std::endl;
        c_mps_before->SaveAs(c_fname_before_png);
        c_mps_before->SaveAs(c_fname_before_pdf);
        //h_mps_before = nullptr;
    }


    ///////////////////////////////////////////////////////////////////////////
    // draw the minimum and draw the point (0,1)
    ///////////////////////////////////////////////////////////////////////////

    #if 0
    if(1)
    {
        params[param_1_ix] = 0.0;
        params[param_2_ix] = 1.0;

        /*
        std::cout << "the parameters before drawing HSD are" << std::endl;
        for(int param_ix = 0; param_ix < params.size(); ++ param_ix)
        {
            std::cout << params[param_ix] << std::endl;
            if(param_ix < 2) continue;
            params[param_ix] = 1.0;
        }
        */

        // TODO: note, if multiple channels enabled then fval is wrong
        // as it should be split between the different channels
        // not printed as the same for both/more than 1

        double fval;
        //logLikelihood(n_params, nullptr, fval, params, 0);
        fval = theFCN.operator()(params);
        std::cout << "fval(" << params[param_1_ix] << "," << params[param_2_ix] << ")=" << fval << std::endl;

        TString savename;
        savename.Form("%s_%d_%d_HSD.png", h_mps_name.Data(), 1, 0);
        //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", true);
        draw(params, param_errs, fval, std::string(savename), ".", g_mode_fake_data, 1);
    }

    if(1)
    {
        //draw_channel(1, params, -1.0, "NOSAVE");
        //std::cin.get();
        //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", true);


        params[param_1_ix] = ll_walk_save.back().first;
        params[param_2_ix] = ll_walk_save.back().second;

        double fval;
        //logLikelihood(n_params, nullptr, fval, params, 0);
        fval = theFCN.operator()(params);
        std::cout << "fval(" << params[param_1_ix] << "," << params[param_2_ix] << ")=" << fval << std::endl;

        TString savename;
        savename.Form("%s_%d_%d_minuit_1_minimum.png", h_mps_name.Data(), 1, 0);
        //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", true);
        draw(params, param_errs, fval, std::string(savename), ".", g_mode_fake_data, 1);
    }

    if(1)
    {
        //draw_channel(1, params, -1.0, "NOSAVE");
        //std::cin.get();
        //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", true);


        params[param_1_ix] = min_x;
        params[param_2_ix] = min_y;

        double fval;
        //logLikelihood(n_params, nullptr, fval, params, 0);
        fval = theFCN.operator()(params);
        std::cout << "fval(" << params[param_1_ix] << "," << params[param_2_ix] << ")=" << fval << std::endl;

        TString savename;
        savename.Form("%s_%d_%d_mps_measured_minimum.png", h_mps_name.Data(), 1, 0);
        //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", true);
        draw(params, param_errs, fval, std::string(savename), ".", g_mode_fake_data, 1);
    }

    if(0)
    {
        #if 0
        params[param_1_ix] = 0.005941;
        params[param_2_ix] = 1.017822;

        //logLikelihood(n_params, nullptr, fval, params, 0);
        fval = theFCN.operator()(params);
        std::cout << "fval(" << params[param_1_ix] << "," << params[param_2_ix] << ")=" << fval << std::endl;

        //TH1D *junk1, *junk2, *junk3, *junk4;
        //TString savename;
        savename.Form("%s_%d_%d_predicted_minimum.png", h_mps_name.Data(), 1, 0);
        //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", true);
        draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", g_mode_fake_data);




        params[param_1_ix] = 0.296;
        params[param_2_ix] = 1.5;

        //logLikelihood(n_params, nullptr, fval, params, 0);
        fval = theFCN.operator()(params);
        std::cout << "fval(" << params[param_1_ix] << "," << params[param_2_ix] << ")=" << fval << std::endl;

        //TH1D *junk1, *junk2, *junk3, *junk4;
        //TString savename;
        savename.Form("%s_%d_%d_expected_scaled_SSD_minimum.png", h_mps_name.Data(), 1, 0);
        //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", true);
        draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", g_mode_fake_data);
        #endif
    }
    #endif
}
