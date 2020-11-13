#ifndef NEWLOGLIKFITTER_MPS_CALCULATE_H
#define NEWLOGLIKFITTER_MPS_CALCULATE_H



void newloglikfitter_mps_calculate
(
    const int number_job_id,
    const std::string &output_name,
    const int start_index,
    const int stop_index
)
{

    // TODO: NOTE: have to change value of "stop_index" as well
    // TODO: this should be part of the scriptX.txt file
    // the script files should contain all data about how the experiment
    // should be run
    const int n_param_xy = 301; //51;//301; // 1001
    int n_param_1 = n_param_xy;
    int n_param_2 = n_param_xy;
    //int n_param_max = n_param_1 * n_param_2;
    int n_param_max = 0;
    if(MODE_PARALLEL == 1)
    {
        n_param_max = (stop_index - start_index) * n_param_2;
    }
    else
    {
        n_param_max = n_param_1 * n_param_2;
    }
    int c_param = 0;


    ///////////////////////////////////////////////////////////////////////////
    // SET param_1 RANGE
    ///////////////////////////////////////////////////////////////////////////

    double param_1_min = 0.0;
    double param_1_max = 0.0;

    // param 1 is gA
    // custom range
//    param_1_min = 0.1; //-0.4; //-0.5; //1.0; //-0.5;
//    param_1_max = 1.7; //0.6; //1.6; //0.5; //2.5; //5.0; //2.5;
    // after changing the psiN0, psiN2 values...
//    param_1_min = -0.3;
//    param_1_max = 1.7; // adjusted for SYS1 energy offset was 1.5
    //param_1_min = -0.4;
    //param_1_max = 1.6; TODO
    // fake data values
//    if(g_mode_fake_data == true)
//    {
////        param_1_min = -0.4;
////        param_1_max = 0.6;
//        param_1_min = -0.7;
//        param_1_max = 0.7;
//    }
    
    // with systematics
//    if(V_ENABLE_SYSALL == true)
//    {
    // TODO: change depending on systematic
    // NOTE: no longer do this because both contours on single plot
    if(g_mode_fake_data == false)
    {
        param_1_min = -0.1;
        param_1_max = 1.9;
    }
    else if(g_mode_fake_data == true)
    {
        //param_1_min = -0.5;
        //param_1_max = 0.7;
        //param_1_min = -0.6;
        //param_1_max = 0.8;
        param_1_min = -0.5;
        param_1_max = 0.7;

        // fit about xi_31=0.7
        param_1_min = -0.1;
        param_1_max = 1.9;
    }
//    }

    // hack to get HSD
    //param_1_min = -0.1;
    //param_1_max = +0.1;

    
    ///////////////////////////////////////////////////////////////////////////
    // SET param_2 RANGE
    ///////////////////////////////////////////////////////////////////////////

    double param_2_min = 0.0;
    double param_2_max = 0.0;
    
    // param 2 is 150Nd amplitude
    // custom range
//    param_2_min = 0.8; //1.1; //0.0; //0.0;
//    param_2_max = 2.6; //2.6; //1.8; //2.0; //2.0; //4.0;
    // after changing the psiN0, psiN2 values...
//    param_2_min = 0.95;
//    param_2_max = 1.3;
    //param_2_min = 0.0;
    //param_2_max = 3.0;  //TODO
    // fake data values
//    if(g_mode_fake_data == true)
//    {
////        param_2_min = 0.2;
////        param_2_max = 1.8
//        param_2_max = 1.15;
//        param_2_min = 0.85;
//    }
    
    // with systematics
    if(g_mode_fake_data == false)
    {
        param_2_min = 0.95;
        param_2_max = 1.45;
    }
    else if(g_mode_fake_data == true)
    {
        //param_2_min = 0.75;
        //param_2_max = 1.25;
        //param_2_min = 0.65;
        //param_2_max = 1.40;
        param_2_min = 0.75;
        param_2_max = 1.25;

        // fit about xi31=0.7
        param_2_min = 0.9;
        param_2_max = 1.4;
    }

    // hack to get HSD
    //param_1_min = 0.9999;
    //param_1_max = 1.0001;


    ///////////////////////////////////////////////////////////////////////////
    // mps ITERATION
    ///////////////////////////////////////////////////////////////////////////

    // mps (after) name
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


    // mps_before name
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






    ///////////////////////////////////////////////////////////////////////////
    // INIT OUTPUT FILE
    ///////////////////////////////////////////////////////////////////////////

    std::string os_fname = std::string("mps_")
                         + std::string("JID") + std::to_string(number_job_id)
                         + std::string(".log");
    std::ofstream os(os_fname, std::ios::out | std::ios::app);

    // TODO: other systematics should be able to be saved in different
    // combinations without overwriting
    std::string output_name_append;
    if(g_mode_fake_data == false)
    {
        output_name_append += "_data";
    }
    else if(g_mode_fake_data == true)
    {
        output_name_append += "_fake";
    }
    if(V_ENABLE_SYSALL == false)
    {
        output_name_append += "_STAT";
    }
    else if(V_ENABLE_SYSALL == true)
    {
        output_name_append += "_STATSYS";
    }

    std::string ofs_resultsmatrix_before_fname =
        output_name + output_name_append + "_before" + "_"
        + "JID" + std::to_string(number_job_id)
        + ".txt";

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

    if(!ofs_resultsmatrix_before.is_open())
    {
        std::cout << "Error: could not open " << ofs_resultsmatrix_before_fname << std::endl;
        return;
    }

    if(!ofs_resultsmatrix_after.is_open())
    {
        std::cout << "Error: could not open " << ofs_resultsmatrix_after_fname << std::endl;
        return;
    }

    ///////////////////////////////////////////////////////////////////////////
    // start of header
    ///////////////////////////////////////////////////////////////////////////

    // before
    ofs_resultsmatrix_before << "[ g_mode_fake_data ] " << g_mode_fake_data << std::endl;
    // after
    ofs_resultsmatrix_after << "[ g_mode_fake_data ] " << g_mode_fake_data << std::endl;

    // before
    ofs_resultsmatrix_before << "[ V_ENABLE_SYSALL ] " << V_ENABLE_SYSALL << " ";
    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        ofs_resultsmatrix_before << V_ENABLE_SYSn[i] << " ";
    }
    ofs_resultsmatrix_before << std::endl;

    // after
    ofs_resultsmatrix_after << "[ V_ENABLE_SYSALL ] " << V_ENABLE_SYSALL << " ";
    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        ofs_resultsmatrix_after << V_ENABLE_SYSn[i] << " ";
    }
    ofs_resultsmatrix_after << std::endl;

    // before: save minimum systematic fit points
    ofs_resultsmatrix_before << "[ ENABLE_MIN_POINT_SYSn ] ";
    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        ofs_resultsmatrix_before << ENABLE_MIN_POINT_SYSn[i] << " ";
    }
    ofs_resultsmatrix_before << std::endl;


    // after: save minimum systematic fit points
    ofs_resultsmatrix_after << "[ ENABLE_MIN_POINT_SYSn ] ";
    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        ofs_resultsmatrix_after << ENABLE_MIN_POINT_SYSn[i] << " ";
    }
    ofs_resultsmatrix_after << std::endl;
    
#if 0
    // before 
    ofs_resultsmatrix_before << "[ min_point ] " << min_point[0] << " " << min_point[1] << std::endl;
    ofs_resultsmatrix_before << "[ min_point_fake_data ] " << min_point_fake_data[0] << " " << min_point_fake_data[1] << std::endl;
    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        ofs_resultsmatrix_before << "[ min_point_sys" << i << "_l ]"
                                 << " " << min_point_sysn_l[i][0]
                                 << " " << min_point_sysn_l[i][1] << std::endl;
        ofs_resultsmatrix_before << "[ min_point_sys" << i << "_h ]"
                                 << " " << min_point_sysn_h[i][0]
                                 << " " << min_point_sysn_h[i][1] << std::endl;
    }

    // after
    ofs_resultsmatrix_after << "[ min_point ] " << min_point[0] << " " << min_point[1] << std::endl;
    ofs_resultsmatrix_after << "[ min_point_fake_data ] " << min_point_fake_data[0] << " " << min_point_fake_data[1] << std::endl;
    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        ofs_resultsmatrix_after << "[ min_point_sys" << i << "_l ]"
                                << " " << min_point_sysn_l[i][0]
                                << " " << min_point_sysn_l[i][1] << std::endl;
        ofs_resultsmatrix_after << "[ min_point_sys" << i << "_h ]"
                                << " " << min_point_sysn_h[i][0]
                                << " " << min_point_sysn_h[i][1] << std::endl;
    }
#endif

    // before
    ofs_resultsmatrix_before << "[ n_param_1 ] " << n_param_1 << " " << param_1_min << " " << param_1_max << std::endl;
    ofs_resultsmatrix_before << "[ n_param_2 ] " << n_param_2 << " " << param_2_min << " " << param_2_max << std::endl;

    // after
    ofs_resultsmatrix_after << "[ n_param_1 ] " << n_param_1 << " " << param_1_min << " " << param_1_max << std::endl;
    ofs_resultsmatrix_after << "[ n_param_2 ] " << n_param_2 << " " << param_2_min << " " << param_2_max << std::endl;

    ///////////////////////////////////////////////////////////////////////////
    // end of header
    ///////////////////////////////////////////////////////////////////////////


    // minimum point found for entire mps iteration space
    // after minuit2 fit
    double min = std::numeric_limits<double>::infinity();
    double min_x = -1.0; //-0.085;
    double min_y = -1.0; //0.87;
    
    // minimum point found for entire mps iteration range
    // before minuit2 fit
    double min_before = std::numeric_limits<double>::infinity();
    double min_x_before = -1.0; //-0.085;
    double min_y_before = -1.0; //0.87;

    if(1)
    {
        std::cout << "START_INDEX=" << start_index << " STOP_INDEX=" << stop_index << std::endl;
        std::cout << "n_param_1=" << n_param_1 << " n_param_2=" << n_param_2 << std::endl;
        
        std::chrono::system_clock::time_point master_start_time = std::chrono::high_resolution_clock::now();

        // modify parameters
        //for(int n_1 = 0; n_1 <= n_param_1; ++ n_1)
        for(int n_1 = 0; n_1 < n_param_1; ++ n_1)
        {
            
            std::chrono::system_clock::time_point outloop_start_time = std::chrono::high_resolution_clock::now();

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

            double min_stripe = std::numeric_limits<double>::infinity();
            double min_stripe_y = 0.0;

            double t_param_1 = (((double)n_1 + 0.5) / (double)n_param_1) * (param_1_max - param_1_min) + param_1_min;

            //for(int n_2 = 0; n_2 <= n_param_2; ++ n_2)
            for(int n_2 = 0; n_2 < n_param_2; ++ n_2)
            {

                std::chrono::system_clock::time_point start_time = std::chrono::high_resolution_clock::now();

                double t_param_2 = (((double)n_2 + 0.5) / (double)n_param_2) * (param_2_max - param_2_min) + param_2_min;

                // assuming that xi31 and 150Nd amplitude are free
                // this may break if the parameter_names.lst file is changed
                gNumberFreeParams = 0;

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

                // these do not change between before/after
                int nch = theFCN.nch;
                //int nfp = g_pg.get_number_free_params();
                int nfp = gNumberFreeParams;
                int ndf = nch - nfp;
 
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

                // hack - trigger when min point is inside of current bin
                const double Nd150_A_width = (param_2_max - param_2_min) / (double)(n_param_2);
                const double xi_31_width = (param_1_max - param_1_min) / (double)(n_param_1);
                const double Nd150_A_bin_low_edge = Nd150_A_value - 0.5 * Nd150_A_width;
                const double Nd150_A_bin_high_edge = Nd150_A_value + 0.5 * Nd150_A_width;
                const double xi_31_bin_low_edge = xi_31_value - 0.5 * xi_31_width;
                const double xi_31_bin_high_edge = xi_31_value + 0.5 * xi_31_width;
                double Nd150_A_min_point = 1.3;
                double xi_31_min_point = 0.5;
                /*
                if(g_mode_fake_data == true)
                {
                    Nd150_A_min_point = min_point_fake_data[1];
                    xi_31_min_point = min_point_fake_data[0];
                }
                else
                {
                    Nd150_A_min_point = min_point[1];
                    xi_31_min_point = min_point[0];
                }
                */
                if((Nd150_A_bin_low_edge <= Nd150_A_min_point) &&
                   (Nd150_A_min_point < Nd150_A_bin_high_edge) &&
                   (xi_31_bin_low_edge <= xi_31_min_point) &&
                   (xi_31_min_point < xi_31_bin_high_edge))
                {
                    std::cout << "n_1=" << n_1 << " n_2=" << n_2 << std::endl;
                    std::cout << "Nd150_A_min_point=" << Nd150_A_min_point << std::endl;
                    std::cout << "xi_31_min_point=" << xi_31_min_point << std::endl;
                    std::cout << "Nd150_A_bin_low_edge=" << Nd150_A_bin_low_edge << std::endl;
                    std::cout << "Nd150_A_bin_high_edge=" << Nd150_A_bin_high_edge << std::endl;
                    std::cout << "xi_31_bin_low_edge=" << xi_31_bin_low_edge << std::endl;
                    std::cout << "xi_31_bin_high_edge=" << xi_31_bin_high_edge << std::endl;
                    DRAW_V_PHYS_ = true;
                    draw_V_PHYS_(".");
                    DRAW_V_PHYS_ = false;
                }
                double fval_after = theFCN.operator()(params_after);
                if((Nd150_A_bin_low_edge <= Nd150_A_min_point) &&
                   (Nd150_A_min_point < Nd150_A_bin_high_edge) &&
                   (xi_31_bin_low_edge <= xi_31_min_point) &&
                   (xi_31_min_point < xi_31_bin_high_edge))
                {
                    std::cout << "fval=" << fval_after << std::endl;
                }
                
                // these do not change between before/after
                nch = theFCN.nch;
                //nfp = g_pg.get_number_free_params();
                nfp = gNumberFreeParams;
                ndf = nch - nfp;


                os << "fval_before=" << fval_before << std::endl;
                //double fval_after = theFCN.operator()(params_after); // TODO: this produces a different result to above?
                os << "fval_after=" << fval_after << std::endl;


                /*
                if((n_2 == 0) || (n_2 == n_param_2 - 1))
                {
                    std::cout << "t_param_1=" << t_param_1
                              << " t_param_2=" << t_param_2
                              << " fval_before=" << fval_before
                              << " fval_after=" << fval_after
                              << std::endl;
                }
                */


                //////////////////////////
                // SET MINIMUM TRACKERS //
                //////////////////////////

                // min stripe after
                if(fval_after < min_stripe)
                {
                    min_stripe = fval_after;
                    min_stripe_y = t_param_2;
                }

                // min after
                if(fval_after < min)
                {
                    min = fval_after;
                    min_x = t_param_1;
                    min_y = t_param_2;
                }

                // min before
                if(fval_before < min_before)
                {
                    min_before = fval_before;
                    min_x_before = t_param_1;
                    min_y_before = t_param_2;
                }


                ///////////////////////////////////////////////////////////////
                // SAVE TO FILE
                ///////////////////////////////////////////////////////////////
                // TODO: save min stripe, min, min_before

                ofs_resultsmatrix_before << n_1 << " "
                                         << n_2 << " "
                                         << t_param_1 << " "
                                         << t_param_2 << " "
                                         << fval_before << " ";
                ofs_resultsmatrix_before << params_before.size() << " ";
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
                ofs_resultsmatrix_after << params_after.size() << " ";
                for(int pix = 0; pix < params_after.size(); ++ pix)
                {
                    ofs_resultsmatrix_after << params_after.at(pix) << " "
                                            << param_errs_after.at(pix) << " ";
                }
                ofs_resultsmatrix_after << std::endl;

                //std::cout << "n_1=" << n_1 << " n_2=" << n_2 << std::endl;
                //std::cout << "t_param_1=" << t_param_1 << " t_param_2=" << t_param_2 << std::endl;
                //std::cout << "fval=" << fval_after << std::endl;
                //std::cout << std::endl;

                /*
                std::string mps_output_name_before = "mps_output_singleenergy_before_"
                                            + std::to_string(bin_ix)
                                            + "_"
                                            + std::to_string(bin_iy)
                                            + ".png";
                draw(params_before, param_errs_before, fval_before,
                     number_job_id,
                     mps_output_name_before,
                     "mps_output_before",
                     false, 1);
                */

                /*
                std::string mps_output_name_after = "mps_output_singleenergy_after_"
                                            + std::to_string(bin_ix)
                                            + "_"
                                            + std::to_string(bin_iy)
                                            + ".png";
                draw(params_after, param_errs_after, fval_after,
                     number_job_id,
                     mps_output_name_after,
                     "mps_output_after",
                     false, 1);
                */

                ++ c_param;
                //std::cout << c_param << " / " << n_param_max << std::endl;

                std::chrono::system_clock::time_point end_time = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> runtime_sec = end_time - start_time;
                //if((n_2 == 0)) || (n_2 == n_param_2 - 1))
                if(n_2 == 0) // || (n_2 == n_param_2 - 1))
                {
                    std::cout << "Exec time: " << runtime_sec.count() << " s (n_2=" << n_2 << ")" << std::endl;
                }


            } // for n_2

            std::cout << c_param << " / " << n_param_max << std::endl;
            std::cout << "min_stripe=" << min_stripe << " min_stripe_x=" << t_param_1 << " min_stripe_y=" << min_stripe_y << std::endl;

            std::chrono::system_clock::time_point outloop_end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> outloop_runtime_sec = outloop_end_time - outloop_start_time;
            std::cout << "Outer loop exec time: " << outloop_runtime_sec.count() << " s" << std::endl;

            std::chrono::duration<double> total_runtime_sec = outloop_end_time - master_start_time;
            //double fraction_complete = (double)n_1 / (double)n_param_1;
            //double fraction_togo = 1.0 - fraction_complete;
            // remaining time = run time / number of runs completed * number of runs remaining
            int number_of_runs_completed = 0;
            int number_of_runs_remaining = 0;
            if(MODE_PARALLEL == 1)
            {
                number_of_runs_completed = (n_1 + 1) - start_index;
                number_of_runs_remaining = stop_index - (n_1 + 1);
            }
            else
            {
                number_of_runs_completed = (n_1 + 1);
                number_of_runs_remaining = n_param_1 - (n_1 + 1);
            }
            double remaining_runtime_sec = total_runtime_sec.count() / (double)number_of_runs_completed * (double)(number_of_runs_remaining);
            double remaining_runtime_hr = remaining_runtime_sec / 3600.0;
            std::cout << "ETA: " << remaining_runtime_hr << " h" << std::endl;


        } // for n_1

    }  // if(1)

    os.close();

    ofs_resultsmatrix_before.close();
    ofs_resultsmatrix_after.close();

}


#endif // NEWLOGLIKFITTER_MPS_CALCULATE_H
