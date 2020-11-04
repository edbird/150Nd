#ifndef NEWLOGLIKFITTER_CHISQUARETEST_H
#define NEWLOGLIKFITTER_CHISQUARETEST_H



#if 0
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
            TH1D *hHighEnergy_allMC = nullptr;
            TH1D *hLowEnergy_allMC = nullptr;
            TH1D *hHighEnergy_data = nullptr;
            TH1D *hLowEnergy_data = nullptr;

            draw(AdjustActs, AdjustActs_Err, -1.0,
                    hHighEnergy_allMC, hLowEnergy_allMC, hHighEnergy_data, hLowEnergy_data,
                    std::string("hTotalE_") + std::string(xi_31_str) + std::string(".png"));
            
            draw_2D(AdjustActs, AdjustActs_Err, std::string("hHighLowEnergy_") + std::string(xi_31_str) + std::string(".png"),
                    hHighEnergy_allMC, hLowEnergy_allMC, hHighEnergy_data, hLowEnergy_data);

            //draw_outputdiff(AdjustActs, 0.296, std::string("houtputdiff_") + std::string(xi_31_str) + std::string(".png"), -3);
            draw_outputdiff(AdjustActs, 0.0, std::string("houtputdiff_") + std::string(xi_31_str) + std::string(".png"), -3);
        }

        std::cout << "done, check output folder for figures" << std::endl;
}
#endif


#if 0
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

        TH1D *hHighEnergy_allMC = nullptr;
        TH1D *hLowEnergy_allMC = nullptr;
        TH1D *hHighEnergy_data = nullptr;
        TH1D *hLowEnergy_data = nullptr;
        draw(params, param_errs, fval, hHighEnergy_allMC, hLowEnergy_allMC, hHighEnergy_data, hLowEnergy_data, saveas_filename, ".", false);

        ofstream_testvalue << "value," << test_value << ",chisquare," << fval << std::endl;

        //void logLikelihood(Int_t & /*nPar*/, Double_t* /*grad*/, Double_t &fval, Double_t *p, Int_t /*iflag */)

    }
    ofstream_testvalue.close();
    delete [] test_values;
    //delete [] params;
    //delete [] param_errs;

}
#endif


#if 0
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
#endif




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


























void newloglikfitter_testmyphasespace_newversion(
    //ROOT::Minuit2::MnUserParameterState &theParameterState,
    //MinimizeFCNAxialVector &theFCN//,
//    ROOT::Minuit2::FunctionMinimum &FCN_min
    const int number_job_id,
    const std::string &output_name,
    const int start_index,
    const int stop_index
    )
{


    ///////////////////////////////////////////////////////////////////////////
    // testing - my phase space
    // 2020-08-19: new version, cuts out a lot of the stuff from the code
    // I originally copied from which didn't make sense in this config
    // plot phase space for Nd150 and Mo100 parameters
    ///////////////////////////////////////////////////////////////////////////


    TString c_mps_name_base = "c_mps_after";
    TString c_mps_name = c_mps_name_base;
    
    std::cout << "rendering: " << c_mps_name << std::endl;

    // TODO: NOTE: have to change value of "stop_index" as well
    const int n_param_xy = 31; //51;//301; // 1001
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



    #if 0
    if(0) 
    {
        ///////////////////////////////////////////////////////////////////////////
        // draw the minimum and draw the point (0,1)
        ///////////////////////////////////////////////////////////////////////////

        std::vector<double> params_copy = params;
        std::vector<double> param_errs_copy = param_errs;

        params_copy[param_1_ix] = 0.0;
        params_copy[param_2_ix] = 1.0;

        double fval;
        fval = theFCN.operator()(params_copy);
        //logLikelihood(n_params, nullptr, fval, params_copy, 0);
        std::cout << "fval(" << params_copy[param_1_ix] << "," << params_copy[param_2_ix] << ")=" << fval << std::endl;

        TString savename;
        savename.Form("%s_%d_%d_expected_minimum.png", h_mps_name.Data(), 1, 0);
        draw(params_copy, param_errs_copy, fval, std::string(savename), ".", g_mode_fake_data);


        //draw_channel(1, params_copy, -1.0, "NOSAVE");
        //std::cin.get();
        //draw(params_copy, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", true);


        params_copy[param_1_ix] = 1.651043;//-0.085; //min_x;
        params_copy[param_2_ix] = 0.521986;//0.87; //min_y;

        fval = theFCN.operator()(params_copy);
        //logLikelihood(n_params_copy, nullptr, fval, params, 0);
        std::cout << "fval(" << params_copy[param_1_ix] << "," << params_copy[param_2_ix] << ")=" << fval << std::endl;

        //TString savename;
        savename.Form("%s_%d_%d_minuit_measured_minimum.png", h_mps_name.Data(), 1, 0);
        draw(params_copy, param_errs_copy, fval, std::string(savename), ".", g_mode_fake_data);

        std::cin.get();
    }
    #endif





// TODO: working here. params is not initialized seneibly.
// check index order where MARK is




// 2020-09-26: does not appear to be used?
/*
    // get minimum
    double fval_min = 0.0;
    // TODO: remove all calls to logLikelihood
    //logLikelihood(n_params, nullptr, fval_min, params, 0);
    //fval_min = theFCN.operator()(params);
    fval_min = std::numeric_limits<double>::infinity();
*/

    std::string os_fname = std::string("mps_")
                         + std::string("JID") + std::to_string(number_job_id)
                         + std::string(".log");
    std::ofstream os(os_fname, std::ios::out | std::ios::app);

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

    std::ofstream ofs_resultsmatrix_before(ofs_resultsmatrix_before_fname);
    std::ofstream ofs_resultsmatrix_after(ofs_resultsmatrix_after_fname);

    std::cout << "*****************************************************" << std::endl;
    std::cout << "*****************************************************" << std::endl;
    std::cout << "saving data to " << ofs_resultsmatrix_before_fname << std::endl;
    std::cout << "saving data to " << ofs_resultsmatrix_after_fname << std::endl;
    std::cout << "*****************************************************" << std::endl;
    std::cout << "*****************************************************" << std::endl;


    if(1)
    {
        std::cout << "START_INDEX=" << start_index << " STOP_INDEX=" << stop_index << std::endl;
        std::cout << "n_param_1=" << n_param_1 << " n_param_2=" << n_param_2 << std::endl;

        // modify parameters
        //for(int n_1 = 0; n_1 <= n_param_1; ++ n_1)
        for(int n_1 = 0; n_1 < n_param_1; ++ n_1)
        {

            std::cout << "n_1=" << n_1 << std::endl;

            if((start_index <= n_1) && (n_1 < stop_index))
            {
                // do nothing
            }
            else
            {
                // managed by another job, skip
                continue;
            }

            // TODO: preserve the minuit parameters as we go down in
            // Y axis
            // to speed up the fit, start params using previous params
            // obtained
            // does this introduce a bias?

            int bin_ix = h_mps->GetNbinsX() - n_1;

            double t_param_1 = 0.0;
            t_param_1 = h_mps->GetXaxis()->GetBinCenter(bin_ix);
            //std::cout << "test: t_param_1=" << t_param_1 << std::endl;

            double min_stripe = std::numeric_limits<double>::infinity();
            double min_stripe_y = 0.0;

            //for(int n_2 = 0; n_2 <= n_param_2; ++ n_2)
            for(int n_2 = 0; n_2 < n_param_2; ++ n_2)
            {

                std::chrono::system_clock::time_point start_time = std::chrono::high_resolution_clock::now();

                // TODO: try using GetBinCenter() and looping over bins
                // in combination with Fill method

                //double fval = 0.;

                //double a_1 = (double)n_1 / (double)n_param_1 - 0.5;
                //double a_2 = (double)n_2 / (double)n_param_2 - 0.5;

                //double t_param_1 = param_1 + width_1 * sigma_1 * a_1;
                //double t_param_2 = param_2 + width_2 * sigma_2 * a_2;
                //t_param_1 = param_1_min + a_1 * (param_1_max - param_1_min);
                //t_param_2 = param_2_min + a_2 * (param_2_max - param_2_min);

                //double t_param_1 = 0.0;
                double t_param_2 = 0.0;

                //t_param_1 = h_mps->GetXaxis()->GetBinCenter(1 + n_1);
                //t_param_1 = h_mps->GetXaxis()->GetBinCenter(h_mps->GetNbinsX() - n_1);
                //t_param_2 = h_mps->GetYaxis()->GetBinCenter(1 + n_2);
                int bin_iy = h_mps->GetNbinsY() - n_2;
                t_param_2 = h_mps->GetYaxis()->GetBinCenter(bin_iy);

                //std::cout << "t_param_1=" << t_param_1 << " t_param_2=" << t_param_2 << std::endl;

                //params[param_1_ix] = t_param_1; // TODO: this params object is not used in the _init function
                //params[param_2_ix] = t_param_2;

                //TH1D *junk1, *junk2, *junk3, *junk4;
                //TString savename;
                //savename.Form("%s_%d_%d.png", h_mps_name.Data(), n_1, n_2);
                //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", true);
                //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", g_mode_fake_data);
                //draw_channel(1, params, -1.0, "NOSAVE");

                //std::cin.get();

                //logLikelihood(n_params, nullptr, fval, params, 0);
                
                // do this for test with no additional minimization
                ROOT::Minuit2::MnUserParameterState theParameterStateBefore;
                ROOT::Minuit2::VariableMetricMinimizer theMinimizer;
                MinimizeFCNAxialVector theFCN;

                const double Nd150_A_value = t_param_2;
                const double Nd150_A_error = 0.1;
                const double xi_31_value = t_param_1;
                const double xi_31_error = 0.1;

                fitBackgrounds_phasespace_init(
                    theParameterStateBefore,
                    theMinimizer,
                    Nd150_A_value,
                    Nd150_A_error,
                    xi_31_value,
                    xi_31_error);
                
                // get initial parameters
                std::vector<double> params_before = theParameterStateBefore.Params();
                std::vector<double> param_errs_before = theParameterStateBefore.Errors();
                double fval_before = theFCN.operator()(params_before);
                //int ndf = theFCN.ndf - theParameterStateBefore.VariableParameters();
                int nch = theFCN.nch;
                //int nfp = g_pg.get_number_free_params();
                int nfp = gNumberFreeParams;
                int ndf = nch - nfp;
 
            /*
                std::cout << "fval=" << fval_before << std::endl;
                for(int i = 0; i < params_before.size(); ++ i)
                {
                    std::cout << "i=" << i << " params_before[" << i << "]=" << params_before[i] << std::endl;
                }
            */  
            /*
                std::cout << "V_ENABLE_SYSALL=" << V_ENABLE_SYSALL << std::endl;
                std::cout << "V_ENABLE_SYS1=" << V_ENABLE_SYS1 << std::endl;
                std::cout << "V_ENABLE_SYS2=" << V_ENABLE_SYS2 << std::endl;
                std::cout << "V_ENABLE_SYS3=" << V_ENABLE_SYS3 << std::endl;
                std::cout << "V_ENABLE_SYS4=" << V_ENABLE_SYS4 << std::endl;
            */

                // do minuit2 fit
                ROOT::Minuit2::FunctionMinimum FCN_min =
                fitBackgrounds_phasespace_exec(
                    theParameterStateBefore,
                    theMinimizer,
                    theFCN);

                // get best fit parameters
                ROOT::Minuit2::MnUserParameterState theParameterStateAfter = FCN_min.UserParameters();
                std::vector<double> params_after = theParameterStateAfter.Params(); // TODO this overwrites params
                std::vector<double> param_errs_after = theParameterStateAfter.Errors();
                double fval_after = theFCN.operator()(params_after);
                //ndf = theFCN.ndf - theParameterStateAfter.VariableParameters();
                nch = theFCN.nch;
                nfp = gNumberFreeParams;
                ndf = nch - nfp;

            /*
                std::cout << "fval=" << fval_after << std::endl;
                for(int i = 0; i < params_after.size(); ++ i)
                {
                    std::cout << "i=" << i << " params_after[" << i << "]=" << params_after[i] << std::endl;
                }
                std::cin.get();
            */

                //print_adjustacts(os, params_after, param_errs_after);

                //std::cout << "fval(" << params[param_1_ix] << "," << params[param_2_ix] << ")=" << fval << std::endl;
               //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", true);
                //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", g_mode_fake_data);

                os << "fval_before=" << fval_before << std::endl;
                //fval = FCN_min.Fval();
                //double fval_after = theFCN.operator()(params_after); // TODO: this produces a different result to above?
                os << "fval_after=" << fval_after << std::endl;

                if((n_2 == 0) || (n_2 == n_param_2 - 1))
                {
                    std::cout << "t_param_1=" << t_param_1
                              << " t_param_2=" << t_param_2
                              << " fval_before=" << fval_before
                              << " fval_after=" << fval_after
                              << std::endl;
                }

                /*
                for(int pix = 0; pix < params_after.size(); ++ pix)
                {
                    std::cout << pix << " " << params_before[pix] << " " << params_after[pix] << std::endl;
                }
                std::cin.get();
                */

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

                if(fval_before < min_before)
                {
                    min_before = fval_before;
                    min_x_before = t_param_1;
                    min_y_before = t_param_2;
                }

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
                //double step_1 = width_1 * sigma_1 * (double)1 / (double)n_param_1;
                //double step_2 = width_2 * sigma_2 * (double)1 / (double)n_param_2;
                //h_mps->Fill(t_param_1 + step_1 / 2.0, t_param_2 + step_2 / 2.0, fval - fval_min);
                //h_mps->Fill(t_param_1, t_param_2, fval);
                //h_mps_before->Fill(t_param_1, t_param_2, fval_before);
                h_mps->SetBinContent(bin_ix, bin_iy, fval_after);
                h_mps_before->SetBinContent(bin_ix, bin_iy, fval_before);
                // TODO: fval_min does not appear to always be the minimum


                ofs_resultsmatrix_before << n_1 << " "
                                         << n_2 << " "
                                         << t_param_1 << " "
                                         << t_param_2 << " "
                                         << fval_before << " ";
                for(int pix = 0; pix < params_before.size(); ++ pix)
                {
                    ofs_resultsmatrix_before << params_before.at(pix) << " "
                                             << param_errs_before.at(pix) << " ";
                }
                ofs_resultsmatrix_before << std::endl;

                ofs_resultsmatrix_after << n_1 << " "
                                        << n_2 << " "
                                        << t_param_1 << " "
                                        << t_param_2 << " "
                                        << fval_after << " ";
                for(int pix = 0; pix < params_after.size(); ++ pix)
                {
                    ofs_resultsmatrix_after << params_after.at(pix) << " "
                                            << param_errs_after.at(pix) << " ";
                }
                ofs_resultsmatrix_after << std::endl;


                std::string mps_output_name_before = "mps_output_singleenergy_before_"
                                            + std::to_string(bin_ix)
                                            + "_"
                                            + std::to_string(bin_iy)
                                            + ".png";
                /*
                draw(params_before, param_errs_before, fval_before,
                     number_job_id,
                     mps_output_name_before,
                     "mps_output_before",
                     false, 1);
                */

                std::string mps_output_name_after = "mps_output_singleenergy_after_"
                                            + std::to_string(bin_ix)
                                            + "_"
                                            + std::to_string(bin_iy)
                                            + ".png";
                /*
                draw(params_after, param_errs_after, fval_after,
                     number_job_id,
                     mps_output_name_after,
                     "mps_output_after",
                     false, 1);
                */

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
                //std::cout << c_param << " / " << n_param_max << std::endl;

                std::chrono::system_clock::time_point end_time = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> runtime_sec = end_time - start_time;
                if((n_2 == 0) || (n_2 == n_param_2 - 1))
                {
                    std::cout << "Exec time: " << runtime_sec.count() << " s (n_2=" << n_2 << ")" << std::endl;
                }
            }
            std::cout << c_param << " / " << n_param_max << std::endl;
            std::cout << "min_stripe=" << min_stripe << " min_stripe_x=" << t_param_1 << " min_stripe_y=" << min_stripe_y << std::endl;
        }

    //print_adjustacts(os, params, param_errs);

    } 
    os.close();

    // TODO: this is not working, but I probably do not want a datetimestamp
    //TString datetimestamp_TString = TString(g_datetimestamp_string);
    /*
    TString f_name = TString("h_mps_1_0_")
                   + "JID" + std::to_string(number_job_id) + "_"
                   + datetimestamp_TString
                   + "_singleenergy.root";
    std::cout << "*****************************************************" << std::endl;
    std::cout << "f_name=" << f_name << std::endl;
    std::cout << "is the filename legal?" << std::endl;
    std::cout << "*****************************************************" << std::endl;
    TFile *f = new TFile(f_name, "recreate");
    h_mps->Write();
    h_mps_before->Write();
    f->Close();
    */


    ///////////////////////////////////////////////////////////////////
    // c_mps
    ///////////////////////////////////////////////////////////////////
    if(0 || (MODE_PARALLEL == 0))
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
//        h_mps_contour->SetLineStyle(0); // in drawmps.C but not used
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
        //TLine *lineXc = new TLine(param_1_min, min_y, param_1_max, min_y);
        //TLine *lineYc = new TLine(min_x, param_2_min, min_x, param_2_max);
        TLine *lineXc = new TLine(param_1_min, min_point[1], param_1_max, min_point[1]);
        TLine *lineYc = new TLine(min_point[0], param_2_min, min_point[0], param_2_max);
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

        std::cout << "********************************************************************************" << std::endl;
        std::cout << "                          PLOTTING MARKERS" << std::endl;
        std::cout << "********************************************************************************" << std::endl;


        // SYS n
        // TODO
        /*
        TMarker *mark_min_point_sys1_l = nullptr;
        TMarker *mark_min_point_sys1_h = nullptr;
        TLine *line_min_point_sys1_l = nullptr;
        TLine *line_min_point_sys1_h = nullptr;
        min_point_marker_helper(
            mark_min_point_sys1_l,
            mark_min_point_sys1_h,
            line_min_point_sys1_l,
            line_min_point_sys1_h,
            min_point,
            min_point_fake_data,
            min_point_sys1_l,
            min_point_sys1_h,
            ENABLE_MIN_POINT_SYS1,
            106, kRed, kRed
            );*/

        /*
        // SYS 1
        TMarker *mark_min_point_sys1_l = nullptr;
        TMarker *mark_min_point_sys1_h = nullptr;
        TLine *line_min_point_sys1_l = nullptr;
        TLine *line_min_point_sys1_h = nullptr;
        min_point_marker_helper(
            mark_min_point_sys1_l,
            mark_min_point_sys1_h,
            line_min_point_sys1_l,
            line_min_point_sys1_h,
            min_point,
            min_point_fake_data,
            min_point_sys1_l,
            min_point_sys1_h,
            ENABLE_MIN_POINT_SYS1,
            106, kRed, kRed
            );

        // SYS 2
        TMarker *mark_min_point_sys2_l = nullptr;
        TMarker *mark_min_point_sys2_h = nullptr;
        TLine *line_min_point_sys2_l = nullptr;
        TLine *line_min_point_sys2_h = nullptr;
        min_point_marker_helper(
            mark_min_point_sys2_l,
            mark_min_point_sys2_h,
            line_min_point_sys2_l,
            line_min_point_sys2_h,
            min_point,
            min_point_fake_data,
            min_point_sys2_l,
            min_point_sys2_h,
            ENABLE_MIN_POINT_SYS2,
            22, kOrange, kOrange
            );

        // SYS 3
        TMarker *mark_min_point_sys3_l = nullptr;
        TMarker *mark_min_point_sys3_h = nullptr;
        TLine *line_min_point_sys3_l = nullptr;
        TLine *line_min_point_sys3_h = nullptr;
        min_point_marker_helper(
            mark_min_point_sys3_l,
            mark_min_point_sys3_h,
            line_min_point_sys3_l,
            line_min_point_sys3_h,
            min_point,
            min_point_fake_data,
            min_point_sys3_l,
            min_point_sys3_h,
            ENABLE_MIN_POINT_SYS3,
            23, kYellow, kYellow
            );

        // SYS 4
        TMarker *mark_min_point_sys4_l = nullptr;
        TMarker *mark_min_point_sys4_h = nullptr;
        TLine *line_min_point_sys4_l = nullptr;
        TLine *line_min_point_sys4_h = nullptr;
        min_point_marker_helper(
            mark_min_point_sys4_l,
            mark_min_point_sys4_h,
            line_min_point_sys4_l,
            line_min_point_sys4_h,
            min_point,
            min_point_fake_data,
            min_point_sys4_l,
            min_point_sys4_h,
            ENABLE_MIN_POINT_SYS4,
            45, kGreen, kGreen
            );

        // SYS 5
        TMarker *mark_min_point_sys5_l = nullptr;
        TMarker *mark_min_point_sys5_h = nullptr;
        TLine *line_min_point_sys5_l = nullptr;
        TLine *line_min_point_sys5_h = nullptr;
        min_point_marker_helper(
            mark_min_point_sys5_l,
            mark_min_point_sys5_h,
            line_min_point_sys5_l,
            line_min_point_sys5_h,
            min_point,
            min_point_fake_data,
            min_point_sys5_l,
            min_point_sys5_h,
            ENABLE_MIN_POINT_SYS5,
            43, kMagenta, kMagenta
            );

        if(ENABLE_MIN_POINT_SYS1 == true)
        {
            if((min_point_sys1_l[0] != 0.0) &&
               (min_point_sys1_l[1] != 0.0))
            {
                mark_min_point_sys1_l->Draw();
                line_min_point_sys1_l->Draw();
            }

            if((min_point_sys1_h[0] != 0.0) &&
               (min_point_sys1_h[1] != 0.0))
            {
                mark_min_point_sys1_h->Draw();
                line_min_point_sys1_h->Draw();
            }
        }

        if(ENABLE_MIN_POINT_SYS2 == true)
        {
            if((min_point_sys2_l[0] != 0.0) &&
               (min_point_sys2_l[1] != 0.0))
            {
                mark_min_point_sys2_l->Draw();
                line_min_point_sys2_l->Draw();
            }

            if((min_point_sys2_h[0] != 0.0) &&
               (min_point_sys2_h[1] != 0.0))
            {
                mark_min_point_sys2_h->Draw();
                line_min_point_sys2_h->Draw();
            }
        }

        if(ENABLE_MIN_POINT_SYS3 == true)
        {
            if((min_point_sys3_l[0] != 0.0) &&
               (min_point_sys3_l[1] != 0.0))
            {
                mark_min_point_sys3_l->Draw();
                line_min_point_sys3_l->Draw();
            }

            if((min_point_sys3_h[0] != 0.0) &&
               (min_point_sys3_h[1] != 0.0))
            {
                mark_min_point_sys3_h->Draw();
                line_min_point_sys3_h->Draw();
            }
        }

        if(ENABLE_MIN_POINT_SYS4 == true)
        {
            if((min_point_sys4_l[0] != 0.0) &&
               (min_point_sys4_l[1] != 0.0))
            {
                mark_min_point_sys4_l->Draw();
                line_min_point_sys4_l->Draw();
            }

            if((min_point_sys4_h[0] != 0.0) &&
               (min_point_sys4_h[1] != 0.0))
            {
                mark_min_point_sys4_h->Draw();
                line_min_point_sys4_h->Draw();
            }
        }

        if(ENABLE_MIN_POINT_SYS5 == true)
        {
            if((min_point_sys5_l[0] != 0.0) &&
               (min_point_sys5_l[1] != 0.0))
            {
                mark_min_point_sys5_l->Draw();
                line_min_point_sys5_l->Draw();
            }

            if((min_point_sys5_h[0] != 0.0) &&
               (min_point_sys5_h[1] != 0.0))
            {
                mark_min_point_sys5_h->Draw();
                line_min_point_sys5_h->Draw();
            }
        }
        */

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
    if(0 || (MODE_PARALLEL == 0))
    {
        TString c_mps_name_base_before = "c_mps_before";
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

    TH2D *h_mps = new TH2D("h_mps", "h_mps", n_nd150, p_min_nd150, p_max_nd150, n_mo100, p_min_mo100, p_max_mo100); 
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

    






#if 0
void newloglikfitter_testmyphasespace(
    //ROOT::Minuit2::MnUserParameterState &theParameterState,
    //MinimizeFCNAxialVector &theFCN//,
//    ROOT::Minuit2::FunctionMinimum &FCN_min
    const int number_job_id,
    const std::string &output_name,
    const int start_index,
    const int stop_index
    )
{

    /*
    double* AdjustActs = new double[theParameterState.Params().size()];
    double* AdjustActs_Err = new double[theParameterState.Params().size()];
    for(int i = 0; i < theParameterState.Params().size(); ++ i)
    {
        AdjustActs[i] = theParameterState.Params().at(i);
        AdjustActs_Err[i] = theParameterState.Errors().at(i);
    }
    */

    ///////////////////////////////////////////////////////////////////////////
    // testing - my phase space
    // plot phase space for Nd150 and Mo100 parameters
    ///////////////////////////////////////////////////////////////////////////


    //std::vector<TCanvas*> c_mps_v;
    //std::vector<TH2D*> h_mps_v;

    std::cout << "numberEnabledParams=" << numberEnabledParams << std::endl;
    std::cout << "numberFreeParams=" << free_params.size() << std::endl;

    // loop over all combinations of parameters
    //for(int param_1_ix = 0; param_1_ix < numberEnabledParams; ++ param_1_ix)
    for(int free_params_index_1 = 0; free_params_index_1 < free_params.size(); ++ free_params_index_1)
    {
        for(int free_params_index_2 = 0; free_params_index_2 < free_params_index_1; ++ free_params_index_2)
        {
            // TODO: this is essential
            // for code to behave nice with JOBID
            // (in parallel)
            if(free_params_index_1 != 1)
                continue;
            if(free_params_index_2 != 0)
                continue;

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

            const int n_param_xy = 151;//301; // 1001
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
            param_1_max = 1.5;
            //param_1_min = -0.4;
            //param_1_max = 1.6; TODO
            // fake data values
            if(g_mode_fake_data)
            {
                param_1_min = -0.4;
                param_1_max = 0.6;
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
            if(g_mode_fake_data)
            {
                param_2_min = 0.2;
                param_2_max = 1.8;
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
            TString h_mps_name = h_mps_name_base + "_" + param_1_ix_str_external + "_" + param_2_ix_str_external;

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
            
            TString param_1_name_str = TString(paramNameMap[param_1_ix_external]);
            TString param_2_name_str = TString(paramNameMap[param_2_ix_external]);

            h_mps->GetXaxis()->SetTitle(param_1_name_str);
            h_mps->GetYaxis()->SetTitle(param_2_name_str);

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
            TString h_mps_name_before = h_mps_name_base_before + "_" + param_1_ix_str_external + "_" + param_2_ix_str_external;

            TH2D *h_mps_before = new TH2D(h_mps_name_before, h_mps_name_before,
                                   n_param_1, param_1_min, param_1_max,
                                   n_param_2, param_2_min, param_2_max); 

            h_mps_before->SetContour(1000);
            
            //TString param_1_name_str = TString(paramNameMap[param_1_ix_external]);
            //TString param_2_name_str = TString(paramNameMap[param_2_ix_external]);

            h_mps_before->GetXaxis()->SetTitle(param_1_name_str);
            h_mps_before->GetYaxis()->SetTitle(param_2_name_str);



            #if 0
            if(0) 
            {
                ///////////////////////////////////////////////////////////////////////////
                // draw the minimum and draw the point (0,1)
                ///////////////////////////////////////////////////////////////////////////

                std::vector<double> params_copy = params;
                std::vector<double> param_errs_copy = param_errs;

                params_copy[param_1_ix] = 0.0;
                params_copy[param_2_ix] = 1.0;

                double fval;
                fval = theFCN.operator()(params_copy);
                //logLikelihood(n_params, nullptr, fval, params_copy, 0);
                std::cout << "fval(" << params_copy[param_1_ix] << "," << params_copy[param_2_ix] << ")=" << fval << std::endl;

                TString savename;
                savename.Form("%s_%d_%d_expected_minimum.png", h_mps_name.Data(), 1, 0);
                draw(params_copy, param_errs_copy, fval, std::string(savename), ".", g_mode_fake_data);


                //draw_channel(1, params_copy, -1.0, "NOSAVE");
                //std::cin.get();
                //draw(params_copy, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", true);


                params_copy[param_1_ix] = 1.651043;//-0.085; //min_x;
                params_copy[param_2_ix] = 0.521986;//0.87; //min_y;

                fval = theFCN.operator()(params_copy);
                //logLikelihood(n_params_copy, nullptr, fval, params, 0);
                std::cout << "fval(" << params_copy[param_1_ix] << "," << params_copy[param_2_ix] << ")=" << fval << std::endl;

                //TString savename;
                savename.Form("%s_%d_%d_minuit_measured_minimum.png", h_mps_name.Data(), 1, 0);
                draw(params_copy, param_errs_copy, fval, std::string(savename), ".", g_mode_fake_data);

                std::cin.get();
            }
            #endif





// TODO: working here. params is not initialized seneibly.
// check index order where MARK is





            // get minimum
            double fval_min = 0.0;
            // TODO: remove all calls to logLikelihood
            //logLikelihood(n_params, nullptr, fval_min, params, 0);
            //fval_min = theFCN.operator()(params);
            fval_min = std::numeric_limits<double>::infinity();

            std::string os_fname = std::string("mps_")
                                 + std::string("JID") + std::to_string(number_job_id)
                                 + std::string(".log");
            std::ofstream os(os_fname, std::ios::out | std::ios::app);


            std::string ofs_resultsmatrix_before_fname =
                output_name + "_before" + "_"
                + "JID" + std::to_string(number_job_id)
                + ".txt";

            //"mps_resultsmatrix_after"
            std::string ofs_resultsmatrix_after_fname =
                output_name + "_after" + "_"
                + "JID" + std::to_string(number_job_id)
                + ".txt";

            std::ofstream ofs_resultsmatrix_before(ofs_resultsmatrix_before_fname);
            std::ofstream ofs_resultsmatrix_after(ofs_resultsmatrix_after_fname);

            std::cout << "*****************************************************" << std::endl;
            std::cout << "*****************************************************" << std::endl;
            std::cout << "saving data to " << ofs_resultsmatrix_before_fname << std::endl;
            std::cout << "saving data to " << ofs_resultsmatrix_after_fname << std::endl;
            std::cout << "*****************************************************" << std::endl;
            std::cout << "*****************************************************" << std::endl;


            if(1)
            {
                std::cout << "START_INDEX=" << start_index << " STOP_INDEX=" << stop_index << std::endl;
                std::cout << "n_param_1=" << n_param_1 << " n_param_2=" << n_param_2 << std::endl;

                // modify parameters
                //for(int n_1 = 0; n_1 <= n_param_1; ++ n_1)
                for(int n_1 = 0; n_1 < n_param_1; ++ n_1)
                {
                    if((start_index <= n_1) && (n_1 < stop_index))
                    {
                        // do nothing
                    }
                    else
                    {
                        // managed by another job, skip
                        continue;
                    }

                    // TODO: preserve the minuit parameters as we go down in
                    // Y axis
                    // to speed up the fit, start params using previous params
                    // obtained
                    // does this introduce a bias?

                    int bin_ix = h_mps->GetNbinsX() - n_1;

                    double t_param_1 = 0.0;
                    t_param_1 = h_mps->GetXaxis()->GetBinCenter(bin_ix);
                    //std::cout << "test: t_param_1=" << t_param_1 << std::endl;

                    double min_stripe = std::numeric_limits<double>::infinity();
                    double min_stripe_y = 0.0;

                    //for(int n_2 = 0; n_2 <= n_param_2; ++ n_2)
                    for(int n_2 = 0; n_2 < n_param_2; ++ n_2)
                    {
                        // TODO: try using GetBinCenter() and looping over bins
                        // in combination with Fill method

                        //double fval = 0.;

                        //double a_1 = (double)n_1 / (double)n_param_1 - 0.5;
                        //double a_2 = (double)n_2 / (double)n_param_2 - 0.5;

                        //double t_param_1 = param_1 + width_1 * sigma_1 * a_1;
                        //double t_param_2 = param_2 + width_2 * sigma_2 * a_2;
                        //t_param_1 = param_1_min + a_1 * (param_1_max - param_1_min);
                        //t_param_2 = param_2_min + a_2 * (param_2_max - param_2_min);

                        //double t_param_1 = 0.0;
                        double t_param_2 = 0.0;

                        //t_param_1 = h_mps->GetXaxis()->GetBinCenter(1 + n_1);
                        //t_param_1 = h_mps->GetXaxis()->GetBinCenter(h_mps->GetNbinsX() - n_1);
                        //t_param_2 = h_mps->GetYaxis()->GetBinCenter(1 + n_2);
                        int bin_iy = h_mps->GetNbinsY() - n_2;
                        t_param_2 = h_mps->GetYaxis()->GetBinCenter(bin_iy);

                        //std::cout << "t_param_1=" << t_param_1 << " t_param_2=" << t_param_2 << std::endl;

                        //params[param_1_ix] = t_param_1; // TODO: this params object is not used in the _init function
                        //params[param_2_ix] = t_param_2;

                        //TH1D *junk1, *junk2, *junk3, *junk4;
                        //TString savename;
                        //savename.Form("%s_%d_%d.png", h_mps_name.Data(), n_1, n_2);
                        //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", true);
                        //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", g_mode_fake_data);
                        //draw_channel(1, params, -1.0, "NOSAVE");

                        //std::cin.get();

                        //logLikelihood(n_params, nullptr, fval, params, 0);
                        
                        // do this for test with no additional minimization
                        ROOT::Minuit2::MnUserParameterState theParameterStateBefore;
                        ROOT::Minuit2::VariableMetricMinimizer theMinimizer;
                        MinimizeFCNAxialVector theFCN;

                        const double Nd150_A = t_param_2;
                        const double xi_31 = t_param_1;

                        fitBackgrounds_phasespace_init(
                            theParameterStateBefore,
                            theMinimizer,
                            //AdjustActs, //  TODO: the _init function uses AdjustActs not params
                            //AdjustActs_Err,
                            Nd150_A,
                            xi_31);
                        
                        // get initial parameters
                        std::vector<double> params_before = theParameterStateBefore.Params();
                        std::vector<double> param_errs_before = theParameterStateBefore.Errors();
                        double fval_before = theFCN.operator()(params_before);

                        // do minuit2 fit
                        ROOT::Minuit2::FunctionMinimum FCN_min =
                        fitBackgrounds_phasespace_exec(
                            theParameterStateBefore,
                            theMinimizer,
                            theFCN);

                        // get best fit parameters
                        ROOT::Minuit2::MnUserParameterState theParameterStateAfter = FCN_min.UserParameters();
                        std::vector<double> params_after = theParameterStateAfter.Params(); // TODO this overwrites params
                        std::vector<double> param_errs_after = theParameterStateAfter.Errors();
                        print_adjustacts(os, params_after, param_errs_after);

                        //std::cout << "fval(" << params[param_1_ix] << "," << params[param_2_ix] << ")=" << fval << std::endl;
                       //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", true);
                        //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", g_mode_fake_data);

                        os << "fval_before=" << fval_before << std::endl;
                        //fval = FCN_min.Fval();
                        double fval_after = theFCN.operator()(params_after); // TODO: this produces a different result to above?
                        os << "fval_after=" << fval_after << std::endl;

                        if((n_2 == 0) || (n_2 == n_param_2 - 1))
                        {
                            std::cout << "t_param_1=" << t_param_1
                                      << " t_param_2=" << t_param_2
                                      << " fval_before=" << fval_before
                                      << " fval_after=" << fval_after
                                      << std::endl;
                        }

                        /*
                        for(int pix = 0; pix < params_after.size(); ++ pix)
                        {
                            std::cout << pix << " " << params_before[pix] << " " << params_after[pix] << std::endl;
                        }
                        std::cin.get();
                        */

                        if(fval_after < min)
                        {
                            min = fval_after;
                            min_x = t_param_1;
                            min_y = t_param_2;
                        }

                        if(fval_after < min_stripe)
                        {
                            min_stripe = fval_after;
                            min_stripe_y = t_param_2;
                        }

                        if(fval_before < min_before)
                        {
                            min_before = fval_before;
                            min_x_before = t_param_1;
                            min_y_before = t_param_2;
                        }

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
                        //double step_1 = width_1 * sigma_1 * (double)1 / (double)n_param_1;
                        //double step_2 = width_2 * sigma_2 * (double)1 / (double)n_param_2;
                        //h_mps->Fill(t_param_1 + step_1 / 2.0, t_param_2 + step_2 / 2.0, fval - fval_min);
                        //h_mps->Fill(t_param_1, t_param_2, fval);
                        //h_mps_before->Fill(t_param_1, t_param_2, fval_before);
                        h_mps->SetBinContent(bin_ix, bin_iy, fval_after);
                        h_mps_before->SetBinContent(bin_ix, bin_iy, fval_before);
                        // TODO: fval_min does not appear to always be the minimum


                        ofs_resultsmatrix_before << n_1 << " "
                                                 << n_2 << " "
                                                 << t_param_1 << " "
                                                 << t_param_2 << " "
                                                 << fval_before << " ";
                        for(int pix = 0; pix < params_before.size(); ++ pix)
                        {
                            ofs_resultsmatrix_before << params_before.at(pix) << " "
                                                     << param_errs_before.at(pix) << " ";
                        }
                        ofs_resultsmatrix_before << std::endl;

                        ofs_resultsmatrix_after << n_1 << " "
                                                << n_2 << " "
                                                << t_param_1 << " "
                                                << t_param_2 << " "
                                                << fval_after << " ";
                        for(int pix = 0; pix < params_after.size(); ++ pix)
                        {
                            ofs_resultsmatrix_after << params_after.at(pix) << " "
                                                    << param_errs_after.at(pix) << " ";
                        }
                        ofs_resultsmatrix_after << std::endl;


                        std::string mps_output_name_before = "mps_output_singleenergy_before_"
                                                    + std::to_string(bin_ix)
                                                    + "_"
                                                    + std::to_string(bin_iy)
                                                    + ".png";
                        /*
                        draw(params_before, param_errs_before, fval_before,
                             number_job_id,
                             mps_output_name_before,
                             "mps_output_before",
                             false, 1);
                        */

                        std::string mps_output_name_after = "mps_output_singleenergy_after_"
                                                    + std::to_string(bin_ix)
                                                    + "_"
                                                    + std::to_string(bin_iy)
                                                    + ".png";
                        /*
                        draw(params_after, param_errs_after, fval_after,
                             number_job_id,
                             mps_output_name_after,
                             "mps_output_after",
                             false, 1);
                        */

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
                        //std::cout << c_param << " / " << n_param_max << std::endl;
                    }
                    std::cout << c_param << " / " << n_param_max << std::endl;
                    std::cout << "min_stripe=" << min_stripe << " min_stripe_x=" << t_param_1 << " min_stripe_y=" << min_stripe_y << std::endl;
                }

            //print_adjustacts(os, params, param_errs);

            } 
            os.close();

            // TODO: this is not working, but I probably do not want a datetimestamp
            TString datetimestamp_TString = TString(g_datetimestamp_string);
            /*
            TString f_name = TString("h_mps_1_0_")
                           + "JID" + std::to_string(number_job_id) + "_"
                           + datetimestamp_TString
                           + "_singleenergy.root";
            std::cout << "*****************************************************" << std::endl;
            std::cout << "f_name=" << f_name << std::endl;
            std::cout << "is the filename legal?" << std::endl;
            std::cout << "*****************************************************" << std::endl;
            TFile *f = new TFile(f_name, "recreate");
            h_mps->Write();
            h_mps_before->Write();
            f->Close();
            */


            ///////////////////////////////////////////////////////////////////
            // c_mps
            ///////////////////////////////////////////////////////////////////
            if(0 || (MODE_PARALLEL == 0))
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

                h_mps_contour->Draw("cont2same");
                //TString c_fname_png = c_mps_name + datetimestamp_TString + ".png";
                //TString c_fname_pdf = c_mps_name + datetimestamp_TString + ".pdf";
                TString c_fname = c_mps_name + "_"
                               + "JID" + std::to_string(number_job_id) + "_"
                               + datetimestamp_TString;
                TString c_fname_png = c_fname + ".png";
                TString c_fname_pdf = c_fname + ".pdf";
                std::cout << "*****************************************************" << std::endl;
                std::cout << "c_fname=" << c_fname << std::endl;
                std::cout << "is the filename legal?" << std::endl;
                std::cout << "*****************************************************" << std::endl;

                TString dir;
                if(V_ENABLE_SYS1 == true)
                {
                    dir += "./SYS1/";
                }
                if(V_ENABLE_SYS2 == true)
                {
                    dir += "./SYS2/";
                }
                if(V_ENABLE_SYS3 == true)
                {
                    dir += "./SYS3/";
                }
                if(V_ENABLE_SYS4 == true)
                {
                    dir += "./SYS4/";
                }
                if(V_ENABLE_SYS5 == true)
                {
                    dir += "./SYS5/";
                }
                
                c_mps->SaveAs(dir + c_fname_png);
                c_mps->SaveAs(dir + c_fname_pdf);
                //h_mps = nullptr;
            }


            ///////////////////////////////////////////////////////////////////
            // c_mps_before
            ///////////////////////////////////////////////////////////////////
            if(0 || (MODE_PARALLEL == 0))
            {
                TString c_mps_name_base_before = "c_mps_before";
                TString c_mps_name_before = c_mps_name_base_before + "_" + param_1_ix_str_external + "_" + param_2_ix_str_external;

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
                                       + "JID" + std::to_string(number_job_id) + "_"
                                       + datetimestamp_TString;
                TString c_fname_before_png = c_fname_before + ".png";
                TString c_fname_before_pdf = c_fname_before + ".pdf";
                std::cout << "*****************************************************" << std::endl;
                std::cout << "c_fname_before=" << c_fname_before << std::endl;
                std::cout << "is the filename legal?" << std::endl;
                std::cout << "*****************************************************" << std::endl;

                TString dir;
                if(V_ENABLE_SYS1 == true)
                {
                    dir += "./SYS1/";
                }
                if(V_ENABLE_SYS2 == true)
                {
                    dir += "./SYS2/";
                }
                if(V_ENABLE_SYS3 == true)
                {
                    dir += "./SYS3/";
                }
                if(V_ENABLE_SYS4 == true)
                {
                    dir += "./SYS4/";
                }
                if(V_ENABLE_SYS5 == true)
                {
                    dir += "./SYS5/";
                }
                
                c_mps_before->SaveAs(dir + c_fname_before_png);
                c_mps_before->SaveAs(dir + c_fname_before_pdf);
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

    TH2D *h_mps = new TH2D("h_mps", "h_mps", n_nd150, p_min_nd150, p_max_nd150, n_mo100, p_min_mo100, p_max_mo100); 
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
#endif

#endif // NEWLOGLIKFITTER_CHISQUARETEST_H
