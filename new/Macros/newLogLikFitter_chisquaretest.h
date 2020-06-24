#ifndef NEWLOGLIKFITTER_CHISQUARETEST_H
#define NEWLOGLIKFITTER_CHISQUARETEST_H





void do_test_xi_31_test1(double *const AdjustActs, double* const AdjustActs_Err)
{

        //double xi_31_values[] = {-0.55, 0.0, 0.246, 0.286, 0.295, 0.297, 0.306, 0.346, 1.0, 1.2, 1.5, 10.0, 100.0};
        double xi_31_values[] = {-0.5, 0.0, 0.246, 0.286, 0.295, 0.297, 0.306, 0.346, 1.0, 1.2, 1.5};

        // draw some gA values as output
        //const int i_max = 1;
        int i_max = sizeof(xi_31_values) / sizeof(xi_31_values[0]);
        for(int i = 0; i < i_max; ++ i)
        {

            // -1.5 does not appear to be consistent - see negative values
            /*
            Double_t xi_31_default = 0.296; //+ 1.0;
            Double_t xi_31_half_range = 2.5;

            Double_t xi_31_offset = 0.0;
            Double_t xi_31_min = xi_31_default - xi_31_half_range + xi_31_offset;
            Double_t xi_31_max = xi_31_default + xi_31_half_range + xi_31_offset;
            */
            //Double_t xi_31_min = -0.565;
            //Double_t xi_31_max = -0.550;
            //Double_t xi_31_min = 100.0;
            //Double_t xi_31_max = 1000.0;
            //Double_t xi_31_min = 0.296 - 0.1;
            //Double_t xi_31_max = 0.296 + 0.1;

            //Double_t xi_31_value = ((double)i / (double)i_max) * (xi_31_max - xi_31_min) + xi_31_min;
            Double_t xi_31_value = xi_31_values[i];
            AdjustActs[1] = xi_31_value / xi_31_baseline;

            TString xi_31_str;
            xi_31_str.Form("%f", xi_31_value);

            // TODO, put in custom directory with text file containing params
            TH1F *hHighEnergy_allMC = nullptr;
            TH1F *hLowEnergy_allMC = nullptr;
            TH1F *hHighEnergy_data = nullptr;
            TH1F *hLowEnergy_data = nullptr;

            draw(AdjustActs, AdjustActs_Err, std::string("hTotalE_") + std::string(xi_31_str) + std::string(".png"),
                    hHighEnergy_allMC, hLowEnergy_allMC, hHighEnergy_data, hLowEnergy_data);
            
            draw_2D(AdjustActs, AdjustActs_Err, std::string("hHighLowEnergy_") + std::string(xi_31_str) + std::string(".png"),
                    hHighEnergy_allMC, hLowEnergy_allMC, hHighEnergy_data, hLowEnergy_data);

            //draw_outputdiff(AdjustActs, 0.296, std::string("houtputdiff_") + std::string(xi_31_str) + std::string(".png"), -3);
            draw_outputdiff(AdjustActs, 0.0, std::string("houtputdiff_") + std::string(xi_31_str) + std::string(".png"), -3);
        }

        std::cout << "done, check output folder for figures" << std::endl;
}



void newloglikfitter_gA_chisquaretest(
    TMinuit *minuit,
    const double* const AdjustActs,
    const double* const AdjustActs_Err)
{

    ///////////////////////////////////////////////////////////////////////////
    // testing
    
    // run chisquare tests

    std::cout << "running chi-square tests (gA): " << "variable: g_A parameter (1)" << std::endl;

    int n_tests = 1; //20;
    // 100 Mo
    int axial_vector_parameter_0_index = paramNumberToMinuitParamNumberMap.at(1);
    std::cout << "the internal index for parameter 1 is " << axial_vector_parameter_0_index << std::endl;
    //std::cin.get();
    // These are in units of minuit internal parameter units
    // To convert to external parameter units, multiply by the value of the
    // external input parameter initial activity
    // Caution: For cases where the fitted parameter minimum is not at 1.0
    // the errors must be treated as upper and lower bound separatly by adding
    // them (internal param & error) to the central value fit parameter
    // external_param_error_lowerbound = (internal_param_CV - internal_param_error) * external_param_init_value
    // similar for upperbound, then subtract and / 2.0
    double test_central_value = AdjustActs[axial_vector_parameter_0_index];
    double test_range = 0.0; //1.0; //10.0 * AdjustActs_Err[axial_vector_parameter_0_index];
    // this range should hit delta sigma = 1.0 at 66 % of the width, but it
    // doesn't.
    double test_start = test_central_value - 0.5 * test_range;
    double test_end   = test_central_value + 0.5 * test_range;
    double *test_values = new double[n_tests];
    double test_step = test_range / (double)n_tests;
    std::cout << "test_central_value=" << test_central_value << "\n"
              << "test_range=" << test_range << "\n"
              << "test_start=" << test_start << "\n"
              << "test_end=" << test_end
              << std::endl;
    int n_params = minuit->GetNumPars();
    double *params = new double[n_params];
    double *param_errs = new double[n_params];
    for(int jx = 0; jx < n_params; ++ jx)
    {
        minuit->GetParameter(jx, params[jx], param_errs[jx]);
    }
    std::ofstream ofstream_testvalue("testvalue_gA.txt", std::ios::out | std::ios::app);
    timestamp(ofstream_testvalue);
    for(int ix = 0; ix < n_tests; ++ ix)
    {
        test_values[ix] = test_start + test_step * ix;

        // get chisquare value for test
        double fval = 0.;

        // set parameter for 100Mo
        double test_value = test_values[ix];
        params[axial_vector_parameter_0_index] = test_value;
        
        std::cout << "test: ix=" << ix << ", " << "test_value=" << test_value << std::endl; //  ", "; << "fval=" << fval << std::endl;

        // TODO: reenable
        logLikelihood(n_params, nullptr, fval, params, 0);
        std::cout << "fval=" << fval << std::endl;

        // save canvas to file
        std::string saveas_filename("testvalue_gA_");
        saveas_filename += std::to_string(ix) + ".png";
        //draw(AdjustActs, AdjustActs_Err, saveas_filename);

        TH1F *hHighEnergy_allMC = nullptr;
        TH1F *hLowEnergy_allMC = nullptr;
        TH1F *hHighEnergy_data = nullptr;
        TH1F *hLowEnergy_data = nullptr;
        draw(params, param_errs, saveas_filename, hHighEnergy_allMC, hLowEnergy_allMC, hHighEnergy_data, hLowEnergy_data);

        ofstream_testvalue << "value," << test_value << ",chisquare," << fval << std::endl;

        //void logLikelihood(Int_t & /*nPar*/, Double_t* /*grad*/, Double_t &fval, Double_t *p, Int_t /*iflag */)

    }
    ofstream_testvalue.close();
    delete [] test_values;
    //delete [] params;
    //delete [] param_errs;

}




void newloglikfitter_100Mo_chisquaretest(
    TMinuit *minuit,
    const double* const AdjustActs,
    const double* const AdjustActs_Err)
{
    ///////////////////////////////////////////////////////////////////////////
    // testing
    
    // run chisquare tests

    std::cout << "running chi-square tests: " << "variable: 100Mo (10)" << std::endl;

    int n_tests = 250;
    // 100 Mo
    int mo100_99_rot_2n2b_m14_index = paramNumberToMinuitParamNumberMap.at(10);
    std::cout << "the internal index for parameter 10 is " << mo100_99_rot_2n2b_m14_index << std::endl;
    // These are in units of minuit internal parameter units
    // To convert to external parameter units, multiply by the value of the
    // external input parameter initial activity
    // Caution: For cases where the fitted parameter minimum is not at 1.0
    // the errors must be treated as upper and lower bound separatly by adding
    // them (internal param & error) to the central value fit parameter
    // external_param_error_lowerbound = (internal_param_CV - internal_param_error) * external_param_init_value
    // similar for upperbound, then subtract and / 2.0
    double test_central_value = AdjustActs[mo100_99_rot_2n2b_m14_index];
    double test_range = 10.0 * AdjustActs_Err[mo100_99_rot_2n2b_m14_index];
    // this range should hit delta sigma = 1.0 at 66 % of the width, but it
    // doesn't.
    double test_start = test_central_value - 0.5 * test_range;
    double test_end =   test_central_value + 0.5 * test_range;
    double *test_values = new double[n_tests];
    double test_step = test_range / (double)n_tests;
    std::cout << "test_central_value=" << test_central_value << "\n"
              << "test_range=" << test_range << "\n"
              << "test_start=" << test_start << "\n"
              << "test_end=" << test_end
              << std::endl;
    int n_params = minuit->GetNumPars();
    double *params = new double[n_params];
    double *param_errs = new double[n_params];
    for(int jx = 0; jx < n_params; ++ jx)
    {
        minuit->GetParameter(jx, params[jx], param_errs[jx]);
    }
    std::ofstream ofstream_testvalue("testvalue.txt");
    for(int ix = 0; ix < n_tests; ++ ix)
    {

        test_values[ix] = test_start + test_step * ix;

        // get chisquare value for test
        double fval = 0.;

        // set parameter for 100Mo
        double test_value = test_values[ix];
        params[mo100_99_rot_2n2b_m14_index] = test_value;

        std::cout << "test: ix=" << ix << ", " << "test_value=" << test_value << std::endl;

        logLikelihood(n_params, nullptr, fval, params, 0);

        ofstream_testvalue << "value=," << test_value << ",chisquare=," << fval << std::endl;

        //void logLikelihood(Int_t & /*nPar*/, Double_t* /*grad*/, Double_t &fval, Double_t *p, Int_t /*iflag */)

    }
    ofstream_testvalue.close();
    delete [] test_values;
    //delete [] params;
    //delete [] param_errs;
    

    ///////////////////////////////////////////////////////////////////////////
}



void newloglikfitter_testmyphasespace(
    TMinuit *minuit,
    const double* const AdjustActs,
    const double* const AdjustActs_Err)
{

    ///////////////////////////////////////////////////////////////////////////
    // testing - my phase space
    // plot phase space for Nd150 and Mo100 parameters
    ///////////////////////////////////////////////////////////////////////////


    std::vector<TCanvas*> c_mps_v;
    std::vector<TH2D*> h_mps_v;

    std::cout << "numberEnabledParams=" << numberEnabledParams << std::endl;
    std::cout << "numberFreeParams=" << free_params.size() << std::endl;

    // loop over all combinations of parameters
    //for(int param_1_ix = 0; param_1_ix < numberEnabledParams; ++ param_1_ix)
    for(int free_params_index_1 = 0; free_params_index_1 < free_params.size(); ++ free_params_index_1)
    {
        for(int free_params_index_2 = 0; free_params_index_2 < free_params_index_1; ++ free_params_index_2)
        {
            int free_param_1 = free_params.at(free_params_index_1);
            int free_param_2 = free_params.at(free_params_index_2);

            int param_1_ix = paramNumberToMinuitParamNumberMap.at(free_param_1);
            int param_2_ix = paramNumberToMinuitParamNumberMap.at(free_param_2);

            int param_1_ix_external = minuitParamNumberToParamNumberMap.at(param_1_ix);
            int param_2_ix_external = minuitParamNumberToParamNumberMap.at(param_2_ix);


            std::cout << "free_params_index_1=" << free_params_index_1
                      << " free_params_index_2=" << free_params_index_2
                      << " free_param_1=" << free_param_1
                      << " free_param_2=" << free_param_2
                      << " param_1_ix=" << param_1_ix
                      << " param_2_ix=" << param_2_ix
                      << " param_1_ix_external=" << param_1_ix_external
                      << " param_2_ix_external=" << param_2_ix_external
                      << std::endl;
                      //std::cin.get();



            
            TString param_1_ix_str_external;
            param_1_ix_str_external.Form("%i", param_1_ix_external);
            TString param_2_ix_str_external;
            param_2_ix_str_external.Form("%i", param_2_ix_external);

            TString c_mps_name_base = "c_mps";
            TString c_mps_name = c_mps_name_base + "_" + param_1_ix_str_external + "_" + param_2_ix_str_external;
            
            std::cout << "rendering: " << c_mps_name << std::endl;

            TCanvas *c_mps = new TCanvas(c_mps_name, c_mps_name);
            c_mps_v.push_back(c_mps);
            //c_mps = nullptr;

            int n_param_1 = 200;
            int n_param_2 = 1000;
            int n_param_max = n_param_1 * n_param_2;
            int c_param = 0;

            double param_1 = AdjustActs[param_1_ix];
            double sigma_1 = AdjustActs_Err[param_1_ix];
            double width_1 = 5.0;
            double param_1_min = param_1 + width_1 * sigma_1 * (-0.5); //(-n_param_1 / 2);
            double param_1_max = param_1 + width_1 * sigma_1 * 0.5; //(n_param_1 - n_param_1 / 2);

            // param 1 is gA
            // custom range
            param_1_min = -1.0;
            param_1_max = 3.0;
            sigma_1 = 1.0;
            width_1 = param_1_max - param_1_min;
            width_1 *= 2.0;
    

            double param_2 = AdjustActs[param_2_ix];
            double sigma_2 = AdjustActs_Err[param_2_ix];
            double width_2 = 5.0;
            double param_2_min = param_2 + width_2 * sigma_2 * (-0.5); //(-n_param_2 / 2);
            double param_2_max = param_2 + width_2 * sigma_2 * 0.5; //(n_param_2 - n_param_2 / 2);
            
            // param 2 is 150Nd amplitude
            // custom range
            param_2_min = 0.0;
            param_2_max = 4.0;
            sigma_2 = 1.0;
            width_2 = param_1_max - param_1_min;
            width_2 *= 2.0;

            // -1.0, 1.0
            // 0.0 .. 2.5

            // x: 0.42: 0.36 .. 0.56
            // y: 1.6: 1.5 .. 1.7

            // x: -0.6 .. -0.4
            // y: 0 .. 0.2

            TString h_mps_name_base = "h_mps";
            TString h_mps_name = h_mps_name_base + "_" + param_1_ix_str_external + "_" + param_2_ix_str_external;

            std::cout << h_mps_name << " param_1=" << param_1 << " sigma_1=" << sigma_1
                                    << " param_1_min=" << param_1_min << " param_1_max=" << param_1_max
                                    << " param_2=" << param_2 << " sigma_2=" << sigma_2
                                    << " param_2_min=" << param_2_min << " param_2_max=" << param_2_max
                                    << std::endl;

            TH2D *h_mps = new TH2D(h_mps_name, h_mps_name,
                                   n_param_1, param_1_min, param_1_max,
                                   n_param_2, param_2_min, param_2_max); 
            h_mps_v.push_back(h_mps);
            //h_mps = nullptr;

            //h_mps->GetZaxis()->SetRangeUser(0.0, 1.0e+04);
            h_mps->SetContour(1000);
            
            TString param_1_name_str = TString(paramNameMap[param_1_ix_external].c_str());
            TString param_2_name_str = TString(paramNameMap[param_2_ix_external]);

            h_mps->GetXaxis()->SetTitle(param_1_name_str);
            h_mps->GetYaxis()->SetTitle(param_2_name_str);

            // reset params array
            // now code moved to new function, simply use new variables (local)
            int n_params = minuit->GetNumPars();
            double *params = new double[n_params];
            double *param_errs = new double[n_params];
            for(int jx = 0; jx < n_params; ++ jx)
            {
                minuit->GetParameter(jx, params[jx], param_errs[jx]);
            }

            // get minimum
            double fval_min = 0.0;
            logLikelihood(n_params, nullptr, fval_min, params, 0);

            double min = std::numeric_limits<double>::infinity();

            // modify parameters
            //for(int n_1 = 0; n_1 <= n_param_1; ++ n_1)
            for(int n_1 = 0; n_1 < n_param_1; ++ n_1)
            {
                //for(int n_2 = 0; n_2 <= n_param_2; ++ n_2)
                for(int n_2 = 0; n_2 < n_param_2; ++ n_2)
                {
                    // TODO: try using GetBinCenter() and looping over bins
                    // in combination with Fill method

                    double fval = 0.;

                    double a_1 = (double)n_1 / (double)n_param_1 - 0.5;
                    double a_2 = (double)n_2 / (double)n_param_2 - 0.5;

                    double t_param_1 = param_1 + width_1 * sigma_1 * a_1;
                    double t_param_2 = param_2 + width_2 * sigma_2 * a_2;
                    t_param_1 = param_1_min + a_1 * (param_1_max - param_1_min);
                    t_param_2 = param_2_min + a_2 * (param_2_max - param_2_min);


                    //t_param_1 = h_mps->GetXaxis()->GetBinCenter(1 + n_1);
                    t_param_1 = h_mps->GetXaxis()->GetBinCenter(h_mps->GetNbinsX() - n_1);
                    //t_param_2 = h_mps->GetYaxis()->GetBinCenter(1 + n_2);
                    t_param_2 = h_mps->GetYaxis()->GetBinCenter(h_mps->GetNbinsY() - n_2);

                    std::cout << "t_param_1=" << t_param_1 << " t_param_2=" << t_param_2 << std::endl;

                    params[param_1_ix] = t_param_1;
                    params[param_2_ix] = t_param_2;

                    logLikelihood(n_params, nullptr, fval, params, 0);

                    if(fval < min)
                        min = fval;

                    /*
                    if(m == 50)
                    {
                        std::cout << "n=" << n << " a_nd150=" << a_nd150 << " p_nd150=" << p_nd150 << " fval=" << fval << std::endl;
                    }
                    */

                    //std::cout << n << " " << m << " " << fval << std::endl;
                    //std::cin.get();

                    //h_mps->Fill(n, m, fval);
                    //h_mps->SetBinContent(n_1 + 1, n_2 + 1, fval - fval_min);
                    double step_1 = width_1 * sigma_1 * (double)1 / (double)n_param_1;
                    double step_2 = width_2 * sigma_2 * (double)1 / (double)n_param_2;
                    //h_mps->Fill(t_param_1 + step_1 / 2.0, t_param_2 + step_2 / 2.0, fval - fval_min);
                    h_mps->Fill(t_param_1, t_param_2, fval);
                    // TODO: fval_min does not appear to always be the minimum

                    /*
                    if(fval - fval_min <= 0.0)
                    {
                        std::cout << "dbg1: " << n_1 << " " << n_2 << " " << h_mps->GetBinContent(n_1, n_2) << std::endl;
                    }
                    if(n_1 == n_param_1 / 2)
                    {
                        if(n_2 == n_param_2 / 2)
                        {
                            std::cout << "dbg2: " << n_param_1 / 2 << " " << n_param_2 / 2 << " " << h_mps->GetBinContent(n_1, n_2) << std::endl;
                        }
                    }
                    */

                    ++ c_param;
                    std::cout << c_param << " / " << n_param_max << std::endl;
                }
            }

            h_mps->Draw("colz");
            std::cout << "min=" << min << std::endl;

            h_mps = nullptr;
            c_mps->SaveAs("mps.png");
            
        }
    }


#if 0
    // reset params array
    // now code moved to new function, simply use new variables (local)
    int n_params = minuit->GetNumPars();
    double *params = new double[n_params];
    double *param_errs = new double[n_params];
    for(int jx = 0; jx < n_params; ++ jx)
    {
        minuit->GetParameter(jx, params[jx], param_errs[jx]);
    }

    TCanvas *c_mps = new TCanvas("c_mps", "c_mps");

    int nd150_rot_2n2b_m4_index = paramNumberToMinuitParamNumberMap.at(0);
    int mo100_99_rot_2n2b_m14_index = paramNumberToMinuitParamNumberMap.at(10);
    // TODO: this will fail if index changes, I think it is now wrong

    int n_nd150 = 100;
    int n_mo100 = 100;

    double p_nd150_cv = AdjustActs[nd150_rot_2n2b_m4_index];
    double p_nd150_sigma = AdjustActs_Err[nd150_rot_2n2b_m4_index];
    double w_nd150 = 5.0;
    double p_min_nd150 = p_nd150_cv + w_nd150 * p_nd150_sigma * (-n_nd150 / 2);
    double p_max_nd150 = p_nd150_cv + w_nd150 * p_nd150_sigma * (n_nd150 - n_nd150 / 2);
    
    double p_mo100_cv = AdjustActs[mo100_99_rot_2n2b_m14_index];
    double p_mo100_sigma = AdjustActs_Err[mo100_99_rot_2n2b_m14_index];
    double w_mo100 = 5.0;
    double p_min_mo100 = p_mo100_cv + w_mo100 * p_mo100_sigma * (double)(-n_mo100 / 2);
    double p_max_mo100 = p_mo100_cv + w_mo100 * p_mo100_sigma * (double)(n_mo100 - n_mo100 / 2);
    
    std::cout << p_nd150_cv << ", " << p_nd150_sigma << ", " << p_min_nd150 << ", " << p_max_nd150 << std::endl;
    std::cout << p_mo100_cv << ", " << p_mo100_sigma << ", " << p_min_mo100 << ", " << p_max_mo100 << std::endl;

    TH2F *h_mps = new TH2F("h_mps", "h_mps", n_nd150, p_min_nd150, p_max_nd150, n_mo100, p_min_mo100, p_max_mo100); 
    //h_mps->GetZaxis()->SetRangeUser(0.0, 1.0e+04);
    h_mps->SetContour(256);
    h_mps->GetXaxis()->SetTitle("^{150}Nd");
    h_mps->GetYaxis()->SetTitle("^{100}Mo");

    // get minimum
    double fval_min = 0.0;
    logLikelihood(n_params, nullptr, fval_min, params, 0);

    for(int n = 0; n < n_nd150; ++ n)
    {
        for(int m = 0; m < n_mo100; ++ m)
        {
            double fval = 0.;

            double a_nd150 = (double)n / (double)n_nd150 - 0.5;
            double a_mo100 = (double)m / (double)n_mo100 - 0.5;

            double p_nd150 = p_nd150_cv + w_nd150 * p_nd150_sigma * a_nd150;
            double p_mo100 = p_mo100_cv + w_mo100 * p_mo100_sigma * a_mo100;

            params[nd150_rot_2n2b_m4_index] = p_nd150;
            params[mo100_99_rot_2n2b_m14_index] = p_mo100;

            logLikelihood(n_params, nullptr, fval, params, 0);

            /*
            if(m == 50)
            {
                std::cout << "n=" << n << " a_nd150=" << a_nd150 << " p_nd150=" << p_nd150 << " fval=" << fval << std::endl;
            }
            */

            //std::cout << n << " " << m << " " << fval << std::endl;
            //std::cin.get();

            //h_mps->Fill(n, m, fval);
            h_mps->SetBinContent(n + 1, m + 1, fval - fval_min);
        }
    }

    h_mps->Draw("colz");
#endif

    
    //delete [] params;
    //delete [] param_errs;


}

#endif // NEWLOGLIKFITTER_CHISQUARETEST_H
