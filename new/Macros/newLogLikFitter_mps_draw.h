#ifndef NEWLOGLIKFITTER_MPS_DRAW_H
#define NEWLOGLIKFITTER_MPS_DRAW_H



#include "newLogLikFitter_mps_draw_data.h"
#include "newLogLikFitter_mps_draw_helper.h"


// TODO: not all statemachine variables are included in the output file
// for example V_ENABLE_SYSx and V_ENABLE_SYSALL are missing
// these should be saved/loaded to/from file as well
// create a struct to manage the output file data
// and read-write functions
// this struct should be a copy of all current program data in memory rather
// than a reference

// NOTE: which output file opened is controlled by:
//
// - before_after_string: switches between "before" and "after"
// - mode_fake_data: switches between "data" and "fakedata"
// - V_ENABLE_SYSALL: switches between "STAT" and "STATSYS"
#if 0
void newloglikfitter_mps_draw_loader
(
    const int number_job_id,
    const std::string &output_name,
    std::string &before_after_string,
    const bool filename_enable_sysall,
    const bool filename_mode_fake_data,
    TH2D *&h_mps,
    int &n_param_1, double& param_1_min, double& param_1_max,
    int &n_param_2, double& param_2_min, double& param_2_max,
    bool *enable_min_point_sysn,
    bool &mode_fake_data,
    double *min_point,
    double *min_point_fake_data,
    double *min_point_sysn_l, double *min_point_sysn_h,
    double &min, double &min_x, double &min_y
)
{

    bool ok = false;
    if(before_after_string == "before")
    {
        ok = true;
    }
    else if(before_after_string == "after")
    {
        ok = true;
    }
    if(ok == false)
    {
        std::cout << "ERROR: " << __func__ << std::endl;
        return;
    }




    ///////////////////////////////////////////////////////////////////////////
    // OPEN FILE
    ///////////////////////////////////////////////////////////////////////////

    std::string output_name_append;


    // fakedata / data
    if(filename_mode_fake_data == false)
    {
        output_name_append += "_data";
    }
    else if(filename_mode_fake_data == true)
    {
        output_name_append += "_fake";
    }

    // SYS / STATSYS
    if(filename_enable_sysall == false)
    {
        output_name_append += "_STAT";
    }
    else if(filename_enable_sysall == true)
    {
        output_name_append += "_STATSYS";
    }



    std::string ofs_resultsmatrix_fname =
        output_name + output_name_append + "_" + before_after_string + "_"
        + "JID" + std::to_string(number_job_id)
        + ".txt";

    std::ifstream ofs_resultsmatrix(ofs_resultsmatrix_fname);

    std::cout << "*****************************************************" << std::endl;
    std::cout << "*****************************************************" << std::endl;
    std::cout << "loading data from " << ofs_resultsmatrix_fname << std::endl;
    std::cout << "*****************************************************" << std::endl;
    std::cout << "*****************************************************" << std::endl;

    if(!ofs_resultsmatrix.is_open())
    {
        std::cout << "ERROR: " << __func__ << " could not open file " << ofs_resultsmatrix_fname << std::endl;
    }





    ///////////////////////////////////////////////////////////////////////////
    // READ "BEFORE" / "AFTER"
    ///////////////////////////////////////////////////////////////////////////

    {
        //ofs_resultsmatrix >> n_param_1 >> param_1_min >> param_1_max;
        std::stringstream ss;
        std::string s;
        std::getline(ofs_resultsmatrix, s);
        //std::cout << "s=" << s << std::endl;
        ss << s;
        ss >> n_param_1 >> param_1_min >> param_1_max;
    }
    
    {
        //ofs_resultsmatrix >> n_param_2 >> param_2_min >> param_2_max;
        std::stringstream ss;
        std::string s;
        std::getline(ofs_resultsmatrix, s);
        //std::cout << "s=" << s << std::endl;
        ss << s;
        ss >> n_param_2 >> param_2_min >> param_2_max;
    }
    //int n_param_max = n_param_1 * n_param_2;

    // read file header data
    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        enable_min_point_sysn[i] = false;
    }

    {
        std::stringstream ss;
        std::string s;
        std::getline(ofs_resultsmatrix, s);
        ss << s;
        for(int i = 0; i < N_SYSTEMATICS; ++ i)
        {
            ss >> enable_min_point_sysn[i];
        }
    }

    {
        std::stringstream ss;
        std::string s;
        std::getline(ofs_resultsmatrix, s);
        ss << s;
        ss >> mode_fake_data;
    }
    
    {
        std::stringstream ss;
        std::string s;
        std::getline(ofs_resultsmatrix, s);
        ss << s;
        ss >> min_point[0] >> min_point[1];
    }

    {
        std::stringstream ss;
        std::string s;
        std::getline(ofs_resultsmatrix, s);
        ss << s;
        ss >> min_point_fake_data[0] >> min_point_fake_data[1];
    }

    {
        for(int i = 0; i < N_SYSTEMATICS; ++ i)
        {
            {
                std::stringstream ss;
                std::string s;
                std::getline(ofs_resultsmatrix, s);
                ss << s;
                ss >> min_point_sysn_l[i * 2 + 0] >> min_point_sysn_l[i * 2 + 1];
            }

            {
                std::stringstream ss;
                std::string s;
                std::getline(ofs_resultsmatrix, s);
                ss << s;
                ss >> min_point_sysn_h[i * 2 + 0] >> min_point_sysn_h[i * 2 + 1];
            }
        }
    }



/*
    std::cout << n_param_1 << " " << param_1_min << " " << param_1_max << std::endl;
    std::cout << n_param_2 << " " << param_2_min << " " << param_2_max << std::endl;
    std::cout << ENABLE_MIN_POINT_SYS1 << " "
              << ENABLE_MIN_POINT_SYS2 << " "
              << ENABLE_MIN_POINT_SYS3 << " "
              << ENABLE_MIN_POINT_SYS4 << " "
              << ENABLE_MIN_POINT_SYS5 << std::endl;
    std::cout << g_mode_fake_data << std::endl;
    std::cout << min_point[0] << " " << min_point[1] << std::endl;
    std::cout << min_point_fake_data[0] << " " << min_point_fake_data[1] << std::endl;
    std::cout << min_point_sys1_l[0] << " " << min_point_sys1_l[1] << std::endl;
    std::cout << min_point_sys1_h[0] << " " << min_point_sys1_h[1] << std::endl;
    std::cout << min_point_sys2_l[0] << " " << min_point_sys2_l[1] << std::endl;
    std::cout << min_point_sys2_h[0] << " " << min_point_sys2_h[1] << std::endl;
    std::cout << min_point_sys3_l[0] << " " << min_point_sys3_l[1] << std::endl;
    std::cout << min_point_sys3_h[0] << " " << min_point_sys3_h[1] << std::endl;
    std::cout << min_point_sys4_l[0] << " " << min_point_sys4_l[1] << std::endl;
    std::cout << min_point_sys4_h[0] << " " << min_point_sys4_h[1] << std::endl;
    std::cout << min_point_sys5_l[0] << " " << min_point_sys5_l[1] << std::endl;
    std::cout << min_point_sys5_h[0] << " " << min_point_sys5_h[1] << std::endl;
*/

    ///////////////////////////////////////////////////////////////////////////
    // INIT TH2 OBJECTS
    ///////////////////////////////////////////////////////////////////////////

    // MPS BEFORE/AFTER //

    TString h_mps_name_base;
    // SYS / STATSYS
    if(filename_enable_sysall == false)
    {
        h_mps_name_base += "_STAT";
    }
    else if(filename_enable_sysall == true)
    {
        h_mps_name_base += "_STATSYS";
    }

    // fakedata / data
    if(filename_mode_fake_data == false)
    {
        h_mps_name_base += "_data";
    }
    else if(filename_mode_fake_data == true)
    {
        h_mps_name_base += "_fake";
    }
    TString h_mps_name = h_mps_name_base + "_" + before_after_string;

    h_mps = new TH2D(h_mps_name, h_mps_name,
                           n_param_1, param_1_min, param_1_max,
                           n_param_2, param_2_min, param_2_max); 

    h_mps->SetContour(1000);
    h_mps->GetXaxis()->SetTitle("^{150}Nd Amplitude Scale Factor");
    h_mps->GetYaxis()->SetTitle("#xi_{31}^{2#nu#beta#beta}");


    ///////////////////////////////////////////////////////////////////////////
    // READ FILE BY LINE
    ///////////////////////////////////////////////////////////////////////////


    // minimum point found from scan
    min = std::numeric_limits<double>::infinity();
    min_x = -1.0;
    min_y = -1.0;
 
    // minimum stripe value for each stripe in scan
    double min_stripe = std::numeric_limits<double>::infinity();
    double min_stripe_y = 0.0;


    // reading
    std::size_t line_count = 1;
    int n_1_last = -1;
    int n_2_last = -1;
    double t_param_1, t_param_2;
    while(!ofs_resultsmatrix.eof())
    {
        //std::cin.get();
        std::stringstream ss;
        std::string s;
        std::getline(ofs_resultsmatrix, s);
        //std::cout << "s=" << s << std::endl;
        ++ line_count;

        ss << s;

        int n_1, n_2;
        double fval;
        std::vector<double> params;
        std::vector<double> param_errs;
        std::size_t params_size = 0;

        ss >> n_1 >> n_2;
        ss >> t_param_1 >> t_param_2;
        ss >> fval;
        ss >> params_size;
        //std::cout << "n_1=" << n_1 << " n_2=" << n_2 << std::endl;
        //std::cout << "t_param_1=" << t_param_1 << " t_param_2=" << t_param_2 << std::endl;
        //std::cout << "fval=" << fval << std::endl;
        //std::cout << "params_size=" << params_size << std::endl;
        //std::cout << std::endl;
        //for(;;)
        //{


            for(std::size_t i = 0; i < params_size; ++ i)
            {
                double tmp1, tmp2;
                ss >> tmp1 >> tmp2;
                params.push_back(tmp1);
                param_errs.push_back(tmp2);
            }
            /*
            try
            {
                if(ss.peek() == std::char_traits<wchar_t>::eof())
                {
                    break;
                }
                double tmp1, tmp2;
                ss >> tmp1 >> tmp2;
                params.push_back(tmp1);
                param_errs.push_back(tmp2);
            }
            catch(...)
            {
                break;
            }
            */
        //}
        //std::cout << "line: " << line_count << " params.size()=" << params.size() << " param_errs.size()=" << param_errs.size() << std::endl;

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

        int bin_ix = n_1 + 1;
        int bin_iy = n_2 + 1;

        if(fval < min_stripe)
        {
            min_stripe = fval;
            min_stripe_y = t_param_2;
        }

        if(fval < min)
        {
            min = fval;
            min_x = t_param_1;
            min_y = t_param_2;
        }

        //std::cout << "ix=" << bin_ix << " iy=" << bin_iy << " fval=" << fval << std::endl;
        h_mps->SetBinContent(bin_ix, bin_iy, fval);
        //std::cout << std::endl;

        ++ line_count;
        //std::cout << "line_count=" << line_count << std::endl;
    }
    std::cout << "min_stripe=" << min_stripe << " min_stripe_x=" << t_param_1 << " min_stripe_y=" << min_stripe_y << std::endl;
    std::cout << "read: " << ofs_resultsmatrix_fname << " -> done" << std::endl;

    // return results


    ofs_resultsmatrix.close();






}
#endif




///////////////////////////////////////////////////////////////////////////////
// new drawing function, without colz or z-axis
///////////////////////////////////////////////////////////////////////////////

void newloglikfitter_mps_draw_systematics
(
    const int number_job_id,
    const std::string &output_name,
    const int start_index,
    const int stop_index  // TODO: these may no longer make sense
)
{

    //const bool mode_fake_data_flag = false;
    const bool mode_fake_data_flag = g_mode_fake_data;

    TString c_mps_name_base = "c_mps_final_after";
    // fakedata / data
    if(mode_fake_data_flag == false)
    {
        c_mps_name_base += "_data";
    }
    else if(mode_fake_data_flag == true)
    {
        c_mps_name_base += "_fake";
    }
    TString c_mps_name = c_mps_name_base;
    
    std::cout << "rendering: " << c_mps_name << std::endl;

    // TODO: NOTE: have to change value of "stop_index" as well
    int c_param = 0;



    ///////////////////////////////////////////////////////////////////////////
    // OPEN FILE
    ///////////////////////////////////////////////////////////////////////////

    std::string after_string = "after";

    mpsdrawdata mps_draw_data_after_sysall;
    mps_draw_data_after_sysall.read(output_name, after_string, true, mode_fake_data_flag);

    mpsdrawdata mps_draw_data_after_sysnone;
    mps_draw_data_after_sysnone.read(output_name, after_string, false, mode_fake_data_flag);



    TCanvas *c_mps_after = nullptr;
    TString after_string_TString = TString(after_string);
    newloglikfitter_mps_draw_helper
    (
        c_mps_after,
        mps_draw_data_after_sysall,
        mps_draw_data_after_sysnone
    );


/*
    double min = min_after_sysnone;
    double min_x = min_x_after_sysnone;
    double min_y = min_y_after_sysnone;

    c_mps_after->cd();
    std::cout << "min=" << min << " min_x=" << min_x << " min_y=" << min_y << std::endl;
    //double clevels[3] = {min + 1.0, min + 2.0, min + 3.0};
    double clevels[3] = {min + 2.30, min + 4.61, min + 9.21};
    //double clevels[3] = {2.30, 4.61, 9.21}; // true minimum is 0.0 for HSD
    h_mps_after_sysnone->SetLineColor(kCyan); //kCyan/kGreen
//        h_mps_after_sysnone->SetLineStyle(0); // in drawmps.C but not used
    h_mps_after_sysnone->SetContour(3, clevels);
    h_mps_after_sysnone->Draw("cont3same");
*/


    //TString c_fname_png = c_mps_name + datetimestamp_TString + ".png";
    //TString c_fname_pdf = c_mps_name + datetimestamp_TString + ".pdf";
    TString c_fname = c_mps_name + "_"
                   + "JID" + std::to_string(number_job_id);// + "_"
                   //+ datetimestamp_TString;
    TString c_fname_png = c_fname + ".png";
    TString c_fname_pdf = c_fname + ".pdf";
    TString c_fname_eps = c_fname + ".eps";
    std::cout << "*****************************************************" << std::endl;
    std::cout << "c_fname=" << c_fname << std::endl;
    std::cout << "is the filename legal?" << std::endl;
    std::cout << "*****************************************************" << std::endl;
    c_mps_after->SaveAs(c_fname_png);
    c_mps_after->SaveAs(c_fname_pdf);
    c_mps_after->SaveAs(c_fname_eps);
    //h_mps = nullptr;
    
}




///////////////////////////////////////////////////////////////////////////////
// new drawing function, draws with colz corresponding to STATSYS data
///////////////////////////////////////////////////////////////////////////////

void newloglikfitter_mps_draw_systematics_colz
(
    const int number_job_id,
    const std::string &output_name,
    const int start_index,
    const int stop_index  // TODO: these may no longer make sense
)
{

    //const bool mode_fake_data_flag = false;
    const bool mode_fake_data_flag = g_mode_fake_data;

    TString c_mps_name_base = "c_mps_final_colz_after";
    // fakedata / data
    if(mode_fake_data_flag == false)
    {
        c_mps_name_base += "_data";
    }
    else if(mode_fake_data_flag == true)
    {
        c_mps_name_base += "_fake";
    }
    TString c_mps_name = c_mps_name_base;
    
    std::cout << "rendering: " << c_mps_name << std::endl;

    // TODO: NOTE: have to change value of "stop_index" as well
    int c_param = 0;



    ///////////////////////////////////////////////////////////////////////////
    // OPEN FILE
    ///////////////////////////////////////////////////////////////////////////

    std::string after_string = "after";

    mpsdrawdata mps_draw_data_after_sysall;
    mps_draw_data_after_sysall.read(output_name, after_string, true, mode_fake_data_flag);

    mpsdrawdata mps_draw_data_after_sysnone;
    mps_draw_data_after_sysnone.read(output_name, after_string, false, mode_fake_data_flag);



    TCanvas *c_mps_after = nullptr;
    TString after_string_TString = TString(after_string);
    newloglikfitter_mps_draw_helper_colz
    (
        c_mps_after,
        mps_draw_data_after_sysall,
        mps_draw_data_after_sysnone
    );


/*
    double min = min_after_sysnone;
    double min_x = min_x_after_sysnone;
    double min_y = min_y_after_sysnone;

    c_mps_after->cd();
    std::cout << "min=" << min << " min_x=" << min_x << " min_y=" << min_y << std::endl;
    //double clevels[3] = {min + 1.0, min + 2.0, min + 3.0};
    double clevels[3] = {min + 2.30, min + 4.61, min + 9.21};
    //double clevels[3] = {2.30, 4.61, 9.21}; // true minimum is 0.0 for HSD
    h_mps_after_sysnone->SetLineColor(kCyan); //kCyan/kGreen
//        h_mps_after_sysnone->SetLineStyle(0); // in drawmps.C but not used
    h_mps_after_sysnone->SetContour(3, clevels);
    h_mps_after_sysnone->Draw("cont3same");
*/


    //TString c_fname_png = c_mps_name + datetimestamp_TString + ".png";
    //TString c_fname_pdf = c_mps_name + datetimestamp_TString + ".pdf";
    TString c_fname = c_mps_name + "_"
                   + "JID" + std::to_string(number_job_id);// + "_"
                   //+ datetimestamp_TString;
    TString c_fname_png = c_fname + ".png";
    TString c_fname_pdf = c_fname + ".pdf";
    TString c_fname_eps = c_fname + ".eps";
    std::cout << "*****************************************************" << std::endl;
    std::cout << "c_fname=" << c_fname << std::endl;
    std::cout << "is the filename legal?" << std::endl;
    std::cout << "*****************************************************" << std::endl;
    c_mps_after->SaveAs(c_fname_png);
    c_mps_after->SaveAs(c_fname_pdf);
    c_mps_after->SaveAs(c_fname_eps);
    //h_mps = nullptr;
    
}





#if 0
void newloglikfitter_mps_draw(
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
    int c_param = 0;



    ///////////////////////////////////////////////////////////////////////////
    // OPEN FILE
    ///////////////////////////////////////////////////////////////////////////

    std::string before_string = "before";
    const bool filename_enable_sysall_after = V_ENABLE_SYSALL;
    const bool filename_mode_fake_data_after = g_mode_fake_data;
    TH2D *h_mps_before = nullptr;
    int n_param_1_before; double param_1_min_before; double param_1_max_before;
    int n_param_2_before; double param_2_min_before; double param_2_max_before;
    bool enable_min_point_sysn_before[N_SYSTEMATICS];
    bool mode_fake_data_before;
    double min_point_before[2];
    double min_point_fake_data_before[2];
    double min_point_sysn_l_before[N_SYSTEMATICS][2];
    double min_point_sysn_h_before[N_SYSTEMATICS][2];
    double min_before;
    double min_x_before;
    double min_y_before;
    newloglikfitter_mps_draw_loader(
        number_job_id,
        output_name,
        before_string,
        filename_enable_sysall_after,
        filename_mode_fake_data_after,
        h_mps_before,
        n_param_1_before, param_1_min_before, param_1_max_before,
        n_param_2_before, param_2_min_before, param_2_max_before,
        enable_min_point_sysn_before,
        mode_fake_data_before,
        min_point_before,
        min_point_fake_data_before,
        &min_point_sysn_l_before[0][0], &min_point_sysn_h_before[0][0],
        min_before, min_x_before, min_y_before
    );

    std::string after_string = "after";
    const bool filename_enable_sysall_before = V_ENABLE_SYSALL;
    const bool filename_mode_fake_data_before = g_mode_fake_data;
    TH2D *h_mps_after = nullptr;
    int n_param_1_after; double param_1_min_after; double param_1_max_after;
    int n_param_2_after; double param_2_min_after; double param_2_max_after;
    bool enable_min_point_sysn_after[N_SYSTEMATICS];
    bool mode_fake_data_after;
    double min_point_after[2];
    double min_point_fake_data_after[2];
    double min_point_sysn_l_after[N_SYSTEMATICS][2];
    double min_point_sysn_h_after[N_SYSTEMATICS][2];
    double min_after;
    double min_x_after;
    double min_y_after;
    newloglikfitter_mps_draw_loader(
        number_job_id,
        output_name,
        after_string,
        filename_enable_sysall_before,
        filename_mode_fake_data_before,
        h_mps_after,
        n_param_1_after, param_1_min_after, param_1_max_after,
        n_param_2_after, param_2_min_after, param_2_max_after,
        enable_min_point_sysn_after,
        mode_fake_data_after,
        min_point_after,
        min_point_fake_data_after,
        &min_point_sysn_l_after[0][0], &min_point_sysn_h_after[0][0],
        min_after, min_x_after, min_y_after
    );


    TCanvas *c_mps_before = nullptr;
    TString before_string_TString = TString(before_string);
    TMarker *mark_min_point_sysn_l_before[N_SYSTEMATICS] = 
    {
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr
    };
    TMarker *mark_min_point_sysn_h_before[N_SYSTEMATICS] = 
    {
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr
    };
    TLine *line_min_point_sysn_l_before[N_SYSTEMATICS] = 
    {
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr
    };
    TLine *line_min_point_sysn_h_before[N_SYSTEMATICS] = 
    {
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr
    };
    /*
    newloglikfitter_draw_mps_helper
    (
        number_job_id,
        before_string_TString,
        c_mps_before,
        h_mps_before,
        n_param_1_before, param_1_min_before, param_1_max_before,
        n_param_2_before, param_2_min_before, param_2_max_before,
        enable_min_point_sysn_before,
        mode_fake_data_before,
        min_point_before,
        min_point_fake_data_before,
        &min_point_sysn_l_before[0][0], &min_point_sysn_h_before[0][0],
        mark_min_point_sysn_l_before, mark_min_point_sysn_h_before,
        line_min_point_sysn_l_before, line_min_point_sysn_h_before,
        min_before, min_x_before, min_y_before
    );
    */
    TCanvas *c_mps_after = nullptr;
    TString after_string_TString = TString(after_string);
    TMarker *mark_min_point_sysn_l_after[N_SYSTEMATICS] = 
    {
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr
    };
    TMarker *mark_min_point_sysn_h_after[N_SYSTEMATICS] = 
    {
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr
    };
    TLine *line_min_point_sysn_l_after[N_SYSTEMATICS] = 
    {
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr
    };
    TLine *line_min_point_sysn_h_after[N_SYSTEMATICS] = 
    {
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr
    };
    /*
    newloglikfitter_draw_mps_helper
    (
        number_job_id,
        after_string_TString,
        c_mps_after,
        h_mps_after,
        n_param_1_after, param_1_min_after, param_1_max_after,
        n_param_2_after, param_2_min_after, param_2_max_after,
        enable_min_point_sysn_after,
        mode_fake_data_after,
        min_point_after,
        min_point_fake_data_after,
        &min_point_sysn_l_after[0][0], &min_point_sysn_h_after[0][0],
        mark_min_point_sysn_l_after, mark_min_point_sysn_h_after,
        line_min_point_sysn_l_after, line_min_point_sysn_h_after,
        min_after, min_x_after, min_y_after
    );
    */









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
#endif



#endif // NEWLOGLIKFITTER_MPS_DRAW_H
