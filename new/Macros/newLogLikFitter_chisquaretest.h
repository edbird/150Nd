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


struct ThreadData
{
    int threadid;

    double min;
    int n_1_start;
    int n_1_stop;

    int n_param_1;
    int n_param_2;
    int n_param_max_thread;

    double param_1;
    double param_2;

    double width_1;
    double width_2;

    double sigma_1;
    double sigma_2;

    double param_1_max;
    double param_1_min;
    double param_2_max;
    double param_2_min;

    int param_1_ix;
    int param_2_ix;

    TH2D *h_mps;

    int n_params;
    double *params;
};


void threadhandle(void *p)
{

    ThreadData *threaddata = (ThreadData*)p;
    //ThreadData &threaddata = *threaddata_p;
    int threadid = threaddata->threadid;
    double min = threaddata->min;
    int n_1_start = threaddata->n_1_start;
    int n_1_stop = threaddata->n_1_stop;
    int n_param_1 = threaddata->n_param_1;
    int n_param_2 = threaddata->n_param_2;
    int n_param_max_thread = threaddata->n_param_max_thread;
    double param_1 = threaddata->param_1;
    double param_2 = threaddata->param_2;
    double width_1 = threaddata->width_1;
    double width_2 = threaddata->width_2;
    double sigma_1 = threaddata->sigma_1;
    double sigma_2 = threaddata->sigma_2;
    double param_1_max = threaddata->param_1_max;
    double param_1_min = threaddata->param_1_min;
    double param_2_max = threaddata->param_2_max;
    double param_2_min = threaddata->param_2_min;
    int param_1_ix = threaddata->param_1_ix;
    int param_2_ix = threaddata->param_2_ix;
    TH2D *h_mps = threaddata->h_mps;
    int n_params = threaddata->n_params;
    double *params = threaddata->params;

    int c_param = 0;


    TThread::Printf("MSG from THREADID=%d: n_1_start=%d n_1_stop=%d\n", threadid, n_1_start, n_1_stop);

    // modify parameters
    //for(int n_1 = 0; n_1 <= n_param_1; ++ n_1)
    //for(int n_1 = 0; n_1 < n_param_1; ++ n_1)
    for(int n_1 = n_1_start; n_1 < n_1_stop; ++ n_1)
    {
        //for(int n_2 = 0; n_2 <= n_param_2; ++ n_2)
        for(int n_2 = 0; n_2 < n_param_2; ++ n_2)
        {
            // TODO: try using GetBinCenter() and looping over bins
            // in combination with Fill method

            double fval = 0.;

            double t_param_1;
            double t_param_2;

            //TThread::Lock();
            //t_param_1 = h_mps->GetXaxis()->GetBinCenter(1 + n_1);
            t_param_1 = h_mps->GetXaxis()->GetBinCenter(h_mps->GetNbinsX() - n_1);
            //t_param_2 = h_mps->GetYaxis()->GetBinCenter(1 + n_2);
            t_param_2 = h_mps->GetYaxis()->GetBinCenter(h_mps->GetNbinsY() - n_2);
            //TThread::UnLock();

//            std::cout << "t_param_1=" << t_param_1 << " t_param_2=" << t_param_2 << std::endl;

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
            //TThread::Lock();
            h_mps->Fill(t_param_1, t_param_2, fval);
            //TThread::UnLock();
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
            if(c_param % 1 == 0)
            {
                TThread::Printf("MSG from THREADID=%d: %d / %d\n", threadid, c_param, n_param_max_thread);
                //std::cout << "MSG from THREADID=" << threadid << ": " <<  c_param << " / " << n_param_max_thread << std::endl;
            }
	    }
    }

    threaddata->min = min;


    TThread::Printf("MSG from THREADID=%d finished\n", threadid);
}

void newloglikfitter_testmyphasespace(
    ROOT::Minuit2::MnUserParameterState &theParameterState,
    MinimizeFCNAxialVector &theFCN,
    ROOT::Minuit2::FunctionMinimum &FCN_min
    )
{

    double* AdjustActs = new double[theParameterState.Params().size()];
    double* AdjustActs_Err = new double[theParameterState.Params().size()];
    for(int i = 0; i < theParameterState.Params().size(); ++ i)
    {
        AdjustActs[i] = theParameterState.Params().at(i);
        AdjustActs_Err[i] = theParameterState.Errors().at(i);
    }

    ///////////////////////////////////////////////////////////////////////////
    // testing - my phase space
    // plot phase space for Nd150 and Mo100 parameters
    ///////////////////////////////////////////////////////////////////////////


    bool mode_fake_data = false; //true;

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

            const int n_param_xy = 301; // 1001
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
            //param_1_min = -0.4;
            //param_1_max = 1.6; TODO
    

            //double param_2 = AdjustActs[param_2_ix];

            double param_2_min;
            double param_2_max;
            
            // param 2 is 150Nd amplitude
            // custom range
            param_2_min = 0.8; //1.1; //0.0; //0.0;
            param_2_max = 2.6; //2.6; //1.8; //2.0; //2.0; //4.0;
            //param_2_min = 0.0;
            //param_2_max = 3.0;  //TODO


            TString h_mps_name_base = "h_mps";
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
            std::vector<double> params = theParameterState.Params();
            std::vector<double> param_errs = theParameterState.Errors();



            double min = std::numeric_limits<double>::infinity();
            double min_x = -1.0; //-0.085;
            double min_y = -1.0; //0.87;





            if(1) 
            {
                ///////////////////////////////////////////////////////////////////////////
                // draw the minimum and draw the point (0,1)
                ///////////////////////////////////////////////////////////////////////////


                params[param_1_ix] = 0.0;
                params[param_2_ix] = 1.0;

                double fval;
                fval = theFCN.operator()(params);
                //logLikelihood(n_params, nullptr, fval, params, 0);
                std::cout << "fval(" << params[param_1_ix] << "," << params[param_2_ix] << ")=" << fval << std::endl;

                TH1D *junk1, *junk2, *junk3, *junk4;
                TString savename;
                savename.Form("%s_%d_%d_expected_minimum.png", h_mps_name.Data(), 1, 0);
                draw(params, param_errs, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", mode_fake_data);


                //draw_channel(1, params, -1.0, "NOSAVE");
                //std::cin.get();
                //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", true);


                params[param_1_ix] = 1.651043;//-0.085; //min_x;
                params[param_2_ix] = 0.521986;//0.87; //min_y;

                fval = theFCN.operator()(params);
                //logLikelihood(n_params, nullptr, fval, params, 0);
                std::cout << "fval(" << params[param_1_ix] << "," << params[param_2_ix] << ")=" << fval << std::endl;

                //TH1D *junk1, *junk2, *junk3, *junk4;
                //TString savename;
                savename.Form("%s_%d_%d_minuit_measured_minimum.png", h_mps_name.Data(), 1, 0);
                draw(params, param_errs, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", mode_fake_data);

                std::cin.get();
            }











            // get minimum
            double fval_min = 0.0;
            // TODO: remove all calls to logLikelihood
            //logLikelihood(n_params, nullptr, fval_min, params, 0);
            fval_min = theFCN.operator()(params);

//<<<<<<< HEAD
//            double min = std::numeric_limits<double>::infinity();
//            const int NUM_THREADS = 2;
//            TThread *thread[NUM_THREADS];
//            ThreadData td[NUM_THREADS];
//            // TODO: can I make this multithreaded
//            for(int i = 0; i < NUM_THREADS; ++ i)
//=======
            if(1)
//>>>>>>> f680bdeea5541e237bbccd1b9f86a5ed017405cf
            {

                td[i].threadid = i;
                td[i].min = min;
                //td[i].n_1_start = n_param_1 / NUM_THREADS * i;
                //td[i].n_1_stop = td[i].n_1_start + n_param_1 / NUM_THREADS - 1;
                td[i].n_1_start = 0 + n_param_1 / NUM_THREADS * i;
                td[i].n_1_stop = n_param_1 / NUM_THREADS * (i + 1);
                td[i].n_param_1 = n_param_1;
                td[i].n_param_2 = n_param_2;
                td[i].n_param_max_thread = n_param_max / NUM_THREADS;
                td[i].param_1 = param_1;
                td[i].param_2 = param_2;
                td[i].width_1 = width_1;
                td[i].width_2 = width_2;
                td[i].sigma_1 = sigma_1;
                td[i].sigma_2 = sigma_2;
                td[i].param_1_max = param_1_max;
                td[i].param_1_min = param_1_min;
                td[i].param_1_ix = param_1_ix;
                td[i].param_2_ix = param_2_ix;
                td[i].h_mps = h_mps;
                td[i].n_params = n_params;
                td[i].params = new double[n_params];
                for(int jx = 0; jx < n_params; ++ jx)
                {
                    td[i].params[jx] = params[jx];
                }

                //std::cout << "Exec thread " << i << std::endl;
                TThread::Printf("Exec thread %d\n", i);
                TString threadname;
                threadname.Form("thread_%i", i);
                thread[i] = new TThread(threadname, threadhandle, (void*)&(td[i]));
                thread[i]->Run();
            }

            TThread::Ps();
            //std::cin.get();

            //std::cout << "joining" << std::endl;
            TThread::Printf("joining\n");
            for(int i = 0; i < NUM_THREADS; ++ i)
            {
            //if(thread[i]->GetState() == kFinishedState)
            //{
                thread[i]->Join();
                //std::cout << "thread i=" << i << " joined" << std::endl;
                TThread::Printf("thread i=%d joined\n", i);
            //}
            // would be nice to do some other work here but not necessary
            }
            //std::cout << "all threads joined" << std::endl;
            TThread::Printf("all threads joined\n");
            TThread::Ps();

            for(int i = 0; i < NUM_THREADS; ++ i)
            {
                //delete thread[i];
            }


#if 0
            // modify parameters
            //for(int n_1 = 0; n_1 <= n_param_1; ++ n_1)
            for(int n_1 = 0; n_1 < n_param_1; ++ n_1)
            {

                double t_param_1 = 0.0;
                t_param_1 = h_mps->GetXaxis()->GetBinCenter(h_mps->GetNbinsX() - n_1);
                std::cout << "test: t_param_1=" << t_param_1 << std::endl;

                double min_stripe = std::numeric_limits<double>::infinity();
                double min_stripe_y = 0.0;

                //for(int n_2 = 0; n_2 <= n_param_2; ++ n_2)
                for(int n_2 = 0; n_2 < n_param_2; ++ n_2)
                {
                    // TODO: try using GetBinCenter() and looping over bins
                    // in combination with Fill method

                    double fval = 0.;

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
                    t_param_2 = h_mps->GetYaxis()->GetBinCenter(h_mps->GetNbinsY() - n_2);

                    //std::cout << "t_param_1=" << t_param_1 << " t_param_2=" << t_param_2 << std::endl;

                    params[param_1_ix] = t_param_1;
                    params[param_2_ix] = t_param_2;

                    TH1D *junk1, *junk2, *junk3, *junk4;
                    TString savename;
                    savename.Form("%s_%d_%d.png", h_mps_name.Data(), n_1, n_2);
                    //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", true);
                    //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", mode_fake_data);
                    //draw_channel(1, params, -1.0, "NOSAVE");

                    //std::cin.get();

                    //logLikelihood(n_params, nullptr, fval, params, 0);
                    fval = theFCN.operator()(params);
                    //std::cout << "fval(" << params[param_1_ix] << "," << params[param_2_ix] << ")=" << fval << std::endl;
                   //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", true);
                    //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", mode_fake_data);

                    if(fval < min)
                    {
                        min = fval;
                        min_x = t_param_1;
                        min_y = t_param_2;
                    }

                    if(fval < min_stripe)
                    {
                        min_stripe = fval;
                        min_stripe_y = t_param_2;
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
                    //std::cout << c_param << " / " << n_param_max << std::endl;
                }
                std::cout << c_param << " / " << n_param_max << std::endl;
                std::cout << "min_stripe=" << min_stripe << " min_stripe_x=" << t_param_1 << " min_stripe_y=" << min_stripe_y << std::endl;
            }
#endif

            } // if(0)

            TFile *f = new TFile("h_mps_1_0_2020-07-15_singleenergy.root", "recreate");
            h_mps->Write();
            f->Close();

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

            h_mps_contour->Draw("cont2same");
            h_mps = nullptr;
            c_mps->SaveAs("mps_2020-07-15.png");
            c_mps->SaveAs("mps_2020-07-15.pdf");


    
            ///////////////////////////////////////////////////////////////////////////
            // draw the minimum and draw the point (0,1)
            ///////////////////////////////////////////////////////////////////////////

            if(1)
            {
                params[param_1_ix] = 0.0;
                params[param_2_ix] = 1.0;

                double fval;
                //logLikelihood(n_params, nullptr, fval, params, 0);
                fval = theFCN.operator()(params);
                std::cout << "fval(" << params[param_1_ix] << "," << params[param_2_ix] << ")=" << fval << std::endl;

                TH1D *junk1, *junk2, *junk3, *junk4;
                TString savename;
                savename.Form("%s_%d_%d_HSD.png", h_mps_name.Data(), 1, 0);
                //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", true);
                draw(params, param_errs, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", mode_fake_data);
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

                TH1D *junk1, *junk2, *junk3, *junk4;
                TString savename;
                savename.Form("%s_%d_%d_minuit_1_minimum.png", h_mps_name.Data(), 1, 0);
                //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", true);
                draw(params, param_errs, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", mode_fake_data);
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

                TH1D *junk1, *junk2, *junk3, *junk4;
                TString savename;
                savename.Form("%s_%d_%d_mps_measured_minimum.png", h_mps_name.Data(), 1, 0);
                //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", true);
                draw(params, param_errs, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", mode_fake_data);
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
                draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", mode_fake_data);




                params[param_1_ix] = 0.296;
                params[param_2_ix] = 1.5;

                //logLikelihood(n_params, nullptr, fval, params, 0);
                fval = theFCN.operator()(params);
                std::cout << "fval(" << params[param_1_ix] << "," << params[param_2_ix] << ")=" << fval << std::endl;

                //TH1D *junk1, *junk2, *junk3, *junk4;
                //TString savename;
                savename.Form("%s_%d_%d_expected_scaled_SSD_minimum.png", h_mps_name.Data(), 1, 0);
                //draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", true);
                draw(params, nullptr, fval, junk1, junk2, junk3, junk4, std::string(savename), ".", mode_fake_data);
                #endif
            }
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

#endif // NEWLOGLIKFITTER_CHISQUARETEST_H
