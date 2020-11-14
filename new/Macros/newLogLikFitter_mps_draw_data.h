#ifndef NEWLOGLIKFITTER_MPS_DRAW_DATA_H
#define NEWLOGLIKFITTER_MPS_DRAW_DATA_H



// NOTE: this always reads from JID0 file
// JID0 now reserved for concatinated output
// script index starts from 1 when in parallel mode
class mpsdrawdata
{

    public:

    mpsdrawdata()
        : number_job_id{0}
        , output_name{""}
        , before_after_flag{"after"}
        , enable_sysall_flag{true}
        , mode_fake_data_flag{true}
        , h_mps{nullptr}
        , n_param_1{0}, param_1_min{0.0}, param_1_max{0.0}
        , n_param_2{0}, param_2_min{0.0}, param_2_max{0.0}
        //, min_point{0.0, 0.0}
        //, min_point_fake_data{0.0, 0.0}
        , enable_min_point_sysn(N_SYSTEMATICS, false)
        //, min_point_sysn_l(N_SYSTEMATICS, {0.0, 0.0})
        //, min_point_sysn_h(N_SYSTEMATICS, {0.0, 0.0})
        //, min_point_sysn_l(N_SYSTEMATICS)
        //, min_point_sysn_h(N_SYSTEMATICS)
        , file_read_mode_fake_data{true}
        , file_read_enable_sysall{true}
        , min{0.0}, min_x{0.0}, min_y{0.0}
        , mark_min_point_sysn_l(N_SYSTEMATICS, nullptr)
        , mark_min_point_sysn_h(N_SYSTEMATICS, nullptr)
        , line_min_point_sysn_l(N_SYSTEMATICS, nullptr)
        , line_min_point_sysn_h(N_SYSTEMATICS, nullptr) 
    {

    }

    virtual
    ~mpsdrawdata()
    {

    }


    int read(
        const std::string& output_name,
        const std::string& before_after_flag,
        const bool enable_sysall_flag,
        const bool mode_fake_data_flag
        )
    {
        this->output_name = output_name;
        this->before_after_flag = before_after_flag;
        this->enable_sysall_flag = enable_sysall_flag;
        this->mode_fake_data_flag = mode_fake_data_flag;


        std::string output_name_append;
        if(mode_fake_data_flag == false)
        {
            output_name_append += "_data";
        }
        else if(mode_fake_data_flag == true)
        {
            output_name_append += "_fake";
        }

        if(enable_sysall_flag == false)
        {
            output_name_append += "_STAT";
        }
        else if(enable_sysall_flag == true)
        {
            output_name_append += "_STATSYS";
        }

    
        std::string output_name_append_2;
        if(mode_fake_data_flag == false)
        {
            output_name_append_2 += "_data";
        }
        else if(mode_fake_data_flag == true)
        {
            output_name_append_2 += "_fake";
        }
    

        std::string ifs_min_point_fname =
            "min_point" + output_name_append_2 + "_";
        std::ifstream ifs_min_point(ifs_min_point_fname);
        double d0, d1;
        ifs_min_point >> d0 >> d1 >> min_fval;
        ifs_min_point.close();

        std::string ifs_min_point_SYSALL_fname =
            "min_point_SYSALL" + output_name_append_2 + "_";
        std::ifstream ifs_min_point_SYSALL(ifs_min_point_SYSALL_fname);
        /*double*/ d0, d1;
        ifs_min_point_SYSALL >> d0 >> d1 >> min_fval_SYSALL;
        ifs_min_point_SYSALL.close();

        std::string ifs_resultsmatrix_fname =
            output_name + output_name_append + "_" + before_after_flag + "_"
            + "JID" + std::to_string(number_job_id)
            + ".txt";

        std::ifstream ifs_resultsmatrix(ifs_resultsmatrix_fname);

        std::cout << "*****************************************************" << std::endl;
        std::cout << "*****************************************************" << std::endl;
        std::cout << "loading data from " << ifs_resultsmatrix_fname << std::endl;
        std::cout << "*****************************************************" << std::endl;
        std::cout << "*****************************************************" << std::endl;

        if(!ifs_resultsmatrix.is_open())
        {
            std::cout << "ERROR: " << __func__ << " could not open file " << ifs_resultsmatrix_fname << std::endl;
            return -1;
        }



        ///////////////////////////////////////////////////////////////////////////
        // READ "BEFORE" / "AFTER"
        ///////////////////////////////////////////////////////////////////////////


        // new header reading code
        // is inside the line-by-line loop 


        // set to "uninit" values
        got_file_read_mode_fake_data = false;
        got_V_ENABLE_SYSALL = false;
        got_V_ENABLE_SYSn = false;
        got_enable_min_point_sysn = false;
        got_min_point = false;
        got_min_point_fake_data = false;
        for(int i = 0; i < N_SYSTEMATICS; ++ i)
        {
            got_min_point_sysn_l[i] = false;
            got_min_point_sysn_h[i] = false;
        }
        got_n_param_1 = false;
        got_n_param_2 = false;




        ///////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////
        // read min points from file
        ///////////////////////////////////////////////////////////////////////
        
        min_point_load("min_point_data", min_point_data, min_point_data_err, min_point_data_fval);  
        min_point_load("min_point_fake", min_point_fake, min_point_fake_err, min_point_fake_fval);  
        min_point_load("min_point_data_SYSALL", min_point_data_SYSALL, min_point_data_SYSALL_err, min_point_data_SYSALL_fval);  
        min_point_load("min_point_fake_SYSALL", min_point_fake_SYSALL, min_point_fake_SYSALL_err, min_point_fake_SYSALL_fval);  
        min_point_load("min_point_data_HSD_CH0", min_point_data_HSD_CH0, min_point_data_HSD_CH0_err, min_point_data_HSD_CH0_fval);  
        min_point_load("min_point_data_HSD", min_point_data_HSD, min_point_data_HSD_err, min_point_data_HSD_fval);  
        min_point_load("min_point_data_SSD", min_point_data_SSD, min_point_data_SSD_err, min_point_data_SSD_fval);  
        for(int i = 0; i < N_SYSTEMATICS; ++ i)
        {
            TString fname_l;
            fname_l.Form("min_point_fake_SYS%dL", i);
            min_point_load(std::string(fname_l), min_point_fake_sysn_l[i], min_point_fake_sysn_l_err[i], min_point_fake_sysn_l_fval[i]);
            TString fname_h;
            fname_h.Form("min_point_fake_SYS%dH", i);
            min_point_load(std::string(fname_h), min_point_fake_sysn_h[i], min_point_fake_sysn_h_err[i], min_point_fake_sysn_h_fval[i]);
        }

        std::cout << "min_point_data=" << min_point_data[0] << " " << min_point_data[1] << std::endl;
        std::cout << "min_point_fake=" << min_point_fake[0] << " " << min_point_fake[1] << std::endl;
        std::cout << "min_point_SYSALL_data=" << min_point_data_SYSALL[0] << " " << min_point_data_SYSALL[1] << std::endl;
        std::cout << "min_point_SYSALL_fake=" << min_point_fake_SYSALL[0] << " " << min_point_fake_SYSALL[1] << std::endl;
        std::cout << "min_point_data_CH0=" << min_point_data_HSD_CH0[0] << " " << min_point_data_HSD_CH0[1] << std::endl;
        std::cout << "min_point_data_HSD=" << min_point_data_HSD[0] << " " << min_point_data_HSD[1] << std::endl;
        std::cout << "min_point_data_SSD=" << min_point_data_SSD[0] << " " << min_point_data_SSD[1] << std::endl;
        for(int i = 0; i < N_SYSTEMATICS; ++ i)
        {
            std::cout << "min_point_fake_SYS" << i << "L=" << min_point_fake_sysn_l[i][0] << " " << min_point_fake_sysn_l[i][1] << std::endl;
            std::cout << "min_point_fake_SYS" << i << "H=" << min_point_fake_sysn_h[i][0] << " " << min_point_fake_sysn_h[i][1] << std::endl;
        }

        ///////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////
        










    #if 0
        // read file header data

        // read statemachine mode (data/fakedata)
        ifs_resultsmatrix >> file_read_mode_fake_data;

        // read statemachine mode (systematics)
        ifs_resultsmatrix >> V_ENABLE_SYSALL;
        for(int i = 0; i < N_SYSTEMATICS; ++ i)
        {
            ifs_resultsmatrix >> V_ENABLE_SYSn[i];
        }

        // read min point enable flags
        for(int i = 0; i < N_SYSTEMATICS; ++ i)
        {
            bool tmp;
            ifs_resultsmatrix >> tmp;
            enable_min_point_sysn[i] = tmp;
        }
        
        // read min point, min point fakedata
        ifs_resultsmatrix >> min_point[0] >> min_point[1];
        ifs_resultsmatrix >> min_point_fake_data[0] >> min_point_fake_data[1];
        
        // read min point for each systematic
        for(int i = 0; i < N_SYSTEMATICS; ++ i)
        {
            ifs_resultsmatrix >> min_point_sysn_l[i][0] >> min_point_sysn_l[i][1];
            ifs_resultsmatrix >> min_point_sysn_h[i][0] >> min_point_sysn_h[i][1];
        }

        // read mps parameters
        ifs_resultsmatrix >> n_param_1 >> param_1_min >> param_1_max;
        ifs_resultsmatrix >> n_param_2 >> param_2_min >> param_2_max;
    #endif
        ///////////////////////////////////////////////////////////////////////
        // end of header
        ///////////////////////////////////////////////////////////////////////


        ///////////////////////////////////////////////////////////////////////
        // INIT TH2 OBJECTS
        ///////////////////////////////////////////////////////////////////////

        // MPS BEFORE/AFTER //

        TString h_mps_name_base;
        // SYS / STATSYS
        if(enable_sysall_flag == false)
        {
            h_mps_name_base += "_STAT";
        }
        else if(enable_sysall_flag == true)
        {
            h_mps_name_base += "_STATSYS";
        }

        // fakedata / data
        if(mode_fake_data_flag == false)
        {
            h_mps_name_base += "_data";
        }
        else if(mode_fake_data_flag == true)
        {
            h_mps_name_base += "_fake";
        }
        TString h_mps_name = h_mps_name_base + "_" + before_after_flag;

        //h_mps = new TH2D(h_mps_name, h_mps_name,
        //                       n_param_1, param_1_min, param_1_max,
        //                       n_param_2, param_2_min, param_2_max); 

        //h_mps->SetContour(1000);
        //h_mps->GetXaxis()->SetTitle("^{150}Nd Amplitude Scale Factor");
        //h_mps->GetYaxis()->SetTitle("#xi_{31}^{2#nu#beta#beta}");


        ///////////////////////////////////////////////////////////////////////
        // READ FILE BY LINE
        ///////////////////////////////////////////////////////////////////////


        // minimum point found from scan
        min = std::numeric_limits<double>::infinity();
        min_x = -1.0;
        min_y = -1.0;
     
        // minimum stripe value for each stripe in scan
        double min_stripe = std::numeric_limits<double>::infinity();
        double min_stripe_y = 0.0;


        // reading
        std::size_t line_count = 0;
        int n_1_last = -1;
        int n_2_last = -1;
        double t_param_1, t_param_2;
        while(!ifs_resultsmatrix.eof())
        {

            // lines start from 1
            ++ line_count;

            //++ line_count;


            ///////////////////////////////////////////////////////////////////
            // new code to read file header
            // the file header can appear anywhere in the file due to the
            // concatination process
            //
            ///////////////////////////////////////////////////////////////////

            std::string token;
            ifs_resultsmatrix >> token;
            if(token == "[")
            {
                const std::string TOKEN1 = "min_point_sys";
                ifs_resultsmatrix >> token;
                if(token == "g_mode_fake_data")
                {
                    ifs_resultsmatrix >> token;
                    if(token == "]")
                    {
                        // read statemachine mode (data/fakedata)
                        if(got_file_read_mode_fake_data == false)
                        {
                            ifs_resultsmatrix >> file_read_mode_fake_data;
                            got_file_read_mode_fake_data = true;
                        }
                        else
                        {
                            std::string dummy;
                            ifs_resultsmatrix >> dummy;
                        }
                    }
                    else
                    {
                        std::cout << "ERROR: EXPECTED \"]\"" << std::endl;
                    }
                }
                else if(token == "V_ENABLE_SYSALL")
                {
                    ifs_resultsmatrix >> token;
                    if(token == "]")
                    {
                        // read statemachine mode (systematics)
                        if(got_V_ENABLE_SYSALL == false)
                        {
                            ifs_resultsmatrix >> V_ENABLE_SYSALL;
                            for(int i = 0; i < N_SYSTEMATICS; ++ i)
                            {
                                ifs_resultsmatrix >> V_ENABLE_SYSn[i];
                            }
                            got_V_ENABLE_SYSALL = true;
                        }
                        else
                        {
                            std::string dummy;
                            ifs_resultsmatrix >> dummy;
                            for(int i = 0; i < N_SYSTEMATICS; ++ i)
                            {
                                ifs_resultsmatrix >> dummy;
                            }
                        }
                    }
                    else
                    {
                        std::cout << "ERROR: EXPECTED \"]\"" << std::endl;
                    }
                }
                else if(token == "ENABLE_MIN_POINT_SYSn")
                {
                    ifs_resultsmatrix >> token;
                    if(token == "]")
                    {
                        // read min point enable flags
                        if(got_enable_min_point_sysn == false)
                        {
                            for(int i = 0; i < N_SYSTEMATICS; ++ i)
                            {
                                bool tmp;
                                ifs_resultsmatrix >> tmp;
                                enable_min_point_sysn[i] = tmp;
                            }
                            got_enable_min_point_sysn = true;
                        }
                        else
                        {
                            std::string dummy;
                            for(int i = 0; i < N_SYSTEMATICS; ++ i)
                            {
                                ifs_resultsmatrix >> dummy;
                            }
                        }
                    }
                    else
                    {
                        std::cout << "ERROR: EXPECTED \"]\"" << std::endl;
                    }
                }
                else if(token == "min_point")
                {
                    // ignore this
                    ifs_resultsmatrix >> token;
                    if(token == "]")
                    {
                        // read min point
                        std::string dummy;
                        ifs_resultsmatrix >> dummy >> dummy;
                    }
                    else
                    {
                        std::cout << "ERROR: EXPECTED \"]\"" << std::endl;
                    }

                    #if 0
                    ifs_resultsmatrix >> token;
                    if(token == "]")
                    {
                        // read min point
                        if(got_min_point == false)
                        {
                            ifs_resultsmatrix >> min_point[0] >> min_point[1];
                            got_min_point = true;
                        }
                        else
                        {
                            std::string dummy;
                            ifs_resultsmatrix >> dummy >> dummy;
                        }
                    }
                    else
                    {
                        std::cout << "ERROR: EXPECTED \"]\"" << std::endl;
                    }
                    #endif
                }
                else if(token == "min_point_fake_data")
                {
                    ifs_resultsmatrix >> token;
                    if(token == "]")
                    {
                        // read min point fakedata
                        std::string dummy;
                        ifs_resultsmatrix >> dummy >> dummy;
                    }
                    else
                    {
                        std::cout << "ERROR: EXPECTED \"]\"" << std::endl;
                    }

                    #if 0
                    ifs_resultsmatrix >> token;
                    if(token == "]")
                    {
                        // read min point fakedata
                        if(got_min_point_fake_data == false)
                        {
                            ifs_resultsmatrix >> min_point_fake_data[0] >> min_point_fake_data[1];
                            got_min_point_fake_data = true;
                        }
                        else
                        {
                            std::string dummy;
                            ifs_resultsmatrix >> dummy >> dummy;
                        }
                    }
                    else
                    {
                        std::cout << "ERROR: EXPECTED \"]\"" << std::endl;
                    }
                    #endif
                }
                else if(token.substr(0, TOKEN1.size()) == TOKEN1)
                {
                    ifs_resultsmatrix >> token;
                    if(token == "]")
                    {
                        // read min point fakedata
                        std::string dummy;
                        ifs_resultsmatrix >> dummy >> dummy;
                    }
                    else
                    {
                        std::cout << "ERROR: EXPECTED \"]\"" << std::endl;
                    }

                    #if 0
                    std::size_t find_pos = token.rfind("_"); // find last underscore position
                    std::string number_str = token.substr(TOKEN1.size(), find_pos - TOKEN1.size());
                    std::string underscore_to_end = token.substr(find_pos);
                    int index = std::stoi(number_str);
                    ifs_resultsmatrix >> token;
                    if(token == "]")
                    {
                        if(underscore_to_end == "_l")
                        {
                            if(got_min_point_sysn_l[index] == false)
                            {
                                ifs_resultsmatrix >> min_point_sysn_l[index][0] >> min_point_sysn_l[index][1];
                                got_min_point_sysn_l[index] = true;
                            }
                            else
                            {
                                std::string dummy;
                                ifs_resultsmatrix >> dummy >> dummy;
                            }
                        }
                        else if(underscore_to_end == "_h")
                        {
                            if(got_min_point_sysn_h[index] == false)
                            {
                                ifs_resultsmatrix >> min_point_sysn_h[index][0] >> min_point_sysn_h[index][1];
                                got_min_point_sysn_h[index] = true;
                            }
                            else
                            {
                                std::string dummy;
                                ifs_resultsmatrix >> dummy >> dummy;
                            }
                        }
                        else
                        {
                            std::cout << "ERROR: underscore_to_end=" << underscore_to_end << " NOT RECOGNIZED" << std::endl;
                            std::cout << "waiting for user input" << std::endl;
                            std::cin.get();
                        }
                    }
                    else
                    {
                        std::cout << "ERROR: EXPECTED \"]\"" << std::endl;
                    }
                    #endif
                }
                else if(token == "n_param_1")
                {
                    ifs_resultsmatrix >> token;
                    if(token == "]")
                    {
                        // read min point fakedata
                        if(got_n_param_1 == false)
                        {
                            // read mps parameters
                            ifs_resultsmatrix >> n_param_1 >> param_1_min >> param_1_max;
                            got_n_param_1 = true;

                            if(got_n_param_1 && got_n_param_2 && (h_mps == nullptr))
                            {
                                h_mps = new TH2D(h_mps_name, h_mps_name,
                                                       n_param_1, param_1_min, param_1_max,
                                                       n_param_2, param_2_min, param_2_max); 

                                h_mps->SetContour(1000);
                                h_mps->GetXaxis()->SetTitle("^{150}Nd Amplitude Scale Factor");
                                h_mps->GetYaxis()->SetTitle("#xi_{31}^{2#nu#beta#beta}");
                            }
                        }
                        else
                        {
                            std::string dummy;
                            ifs_resultsmatrix >> dummy >> dummy >> dummy;
                        }
                    }
                    else
                    {
                        std::cout << "ERROR: EXPECTED \"]\"" << std::endl;
                    }

                }
                else if(token == "n_param_2")
                {
                    ifs_resultsmatrix >> token;
                    if(token == "]")
                    {
                        // read min point fakedata
                        if(got_n_param_2 == false)
                        {
                            // read mps parameters
                            ifs_resultsmatrix >> n_param_2 >> param_2_min >> param_2_max;
                            got_n_param_2 = true;

                            if(got_n_param_1 && got_n_param_2 && (h_mps == nullptr))
                            {
                                h_mps = new TH2D(h_mps_name, h_mps_name,
                                                       n_param_1, param_1_min, param_1_max,
                                                       n_param_2, param_2_min, param_2_max); 

                                h_mps->SetContour(1000);
                                h_mps->GetXaxis()->SetTitle("^{150}Nd Amplitude Scale Factor");
                                h_mps->GetYaxis()->SetTitle("#xi_{31}^{2#nu#beta#beta}");
                            }
                        }
                        else
                        {
                            std::string dummy;
                            ifs_resultsmatrix >> dummy >> dummy >> dummy;
                        }
                    }
                    else
                    {
                        std::cout << "ERROR: EXPECTED \"]\"" << std::endl;
                    }
                }
                else
                {
                    std::cout << "ERROR: token=" << token << " NOT RECOGNIZED" << std::endl;
                    std::cout << "waiting for user input" << std::endl;
                    std::cin.get();
                }

            }
            else
            {
                int n_1;
                try
                {
                    n_1 = std::stoi(token);
                }
                catch(...)
                {
                    std::cout << "stoi ERROR line: " << line_count << std::endl;
                    std::cout << "token=" << token << std::endl;
                    std::cout << "token.size()=" << token.size() << std::endl;
                    for(int i = 0; i < token.size(); ++ i)
                    {
                        std::cout << "(int)token[" << i << "]=" << (int)token[i] << std::endl;
                    }
                    std::cout << "cin.get()" << std::endl;
                    break;
                    std::cin.get();
                }
                //int n_1, n_2;
                int n_2;
                double fval;
                std::vector<double> params;
                std::vector<double> param_errs;
                std::size_t params_size = 0;

                // read head of line
                //ifs_resultsmatrix >> n_1 >> n_2;
                ifs_resultsmatrix >> n_2;
                ifs_resultsmatrix >> t_param_1 >> t_param_2;
                ifs_resultsmatrix >> fval;
                ifs_resultsmatrix >> params_size;
                params.reserve(params_size);
                param_errs.reserve(params_size);

                // read remainder of line
                for(std::size_t i = 0; i < params_size; ++ i)
                {
                    double tmp1, tmp2;
                    ifs_resultsmatrix >> tmp1 >> tmp2;
                    params.push_back(tmp1);
                    param_errs.push_back(tmp2);
                }
                
                //std::cout << "n_1=" << n_1 << " n_2=" << n_2 << std::endl;
                //std::cout << "t_param_1=" << t_param_1 << " t_param_2=" << t_param_2 << std::endl;
                //std::cout << "fval=" << fval << std::endl;
                //std::cout << "params_size=" << params_size << std::endl;
                //std::cout << std::endl;
                //for(;;)
                //{
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
                        std::cout << "n_1=" << n_1 << " min_stripe=" << min_stripe << " min_stripe_x=" << t_param_1 << " min_stripe_y=" << min_stripe_y << std::endl;

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

                Double_t content_check = h_mps->GetBinContent(bin_ix, bin_iy);
                if(content_check == 0.0)
                {
                    //std::cout << "ix=" << bin_ix << " iy=" << bin_iy << " fval=" << fval << std::endl;
                    h_mps->SetBinContent(bin_ix, bin_iy, fval);
                    //std::cout << std::endl;

                    // have to hope that 0.0 does not appear as a real value
                    // in the mps or the error message in below block will not
                    // be triggered
                }
                else
                {
                    std::cout << "line_count=" << line_count << " - ERROR: " << "mps already contains non-zero value at bin_ix=" << bin_ix << " bin_iy=" << bin_iy << " value=" << content_check << " new value=" << fval << std::endl;
                }

                // duplicated?
                //++ line_count;
                //std::cout << "line_count=" << line_count << std::endl;

            }



        }
        std::cout << "min_stripe=" << min_stripe << " min_stripe_x=" << t_param_1 << " min_stripe_y=" << min_stripe_y << std::endl;
        std::cout << "read: " << ifs_resultsmatrix_fname << " -> done" << std::endl;

        // return results


        ifs_resultsmatrix.close();


        ///////////////////////////////////////////////////////////////////////
        // init contour clone
        ///////////////////////////////////////////////////////////////////////


        ///////////////////////////////////////////////////////////////////////
        // init min point markers
        ///////////////////////////////////////////////////////////////////////

        // 0.1 MeV (ignore)
        // 1.2 %
        // 5.55 %
        // 0.5 % (ignore)
        // 3 keV 
        // foil thickness
        // energy loss
        // brem (ignore)
        
        Int_t markerstylen[N_SYSTEMATICS] = 
        {
            50, // "+"
            22, // triangle up
            43, // star
            52, // "x"
            20, // circle
            33, // diamond
            23, // triangle down
            21, // square
            33, // diamond
            23, // triangle down
            34  // plus
        };
        Double_t markersizen[N_SYSTEMATICS] = 
        {
            2.0, 1.5, 2.0, 2.0, 1.0, 2.0, 1.5, 1.0, 2.0, 1.5, 1.5
        };
        Color_t markercolorn[N_SYSTEMATICS] = 
        {
            //kRed, kOrange, kGreen, kBlue, kMagenta, kViolet, kViolet + 10
            kBlack,         // "+"
            kRed,           // triangle up
            kMagenta,       // star
            kBlack,         // "x"
            kBlue,          // circle
            kViolet + 1,    // diamond
            kGreen,         // triangle down
            kBlack,         // square
            kAzure + 1,     // diamond
            kViolet + 10,   // triangle down
            kAzure + 1      // plus
        };

        for(int i = 0; i < N_SYSTEMATICS; ++ i)
        {
            // SYS n
            min_point_marker_helper(
                mark_min_point_sysn_l[i],
                mark_min_point_sysn_h[i],
                line_min_point_sysn_l[i],
                line_min_point_sysn_h[i],
                min_point_data,
                min_point_fake,
                min_point_fake_sysn_l[i],
                min_point_fake_sysn_h[i],
                ENABLE_MIN_POINT_SYSn[i],
                markerstylen[i], markercolorn[i], markersizen[i], markercolorn[i]
                );
        }


        return 0;
    }


    int number_job_id;
    std::string output_name;
    std::string before_after_flag;
    bool enable_sysall_flag;
    bool mode_fake_data_flag;
    TH2D *h_mps;
    int n_param_1; double param_1_min; double param_1_max;
    int n_param_2; double param_2_min; double param_2_max;
    //double min_point[2];
    //double min_point_fake_data[2];
    std::vector<bool> enable_min_point_sysn;
    //std::vector<double[2]> min_point_sysn_l;
    //std::vector<double[2]> min_point_sysn_h;
    bool file_read_mode_fake_data;
    bool file_read_enable_sysall;
    double min; double min_x; double min_y;

    // track whether each line of "header" has been read
    // (header can appear multiple times)
    bool got_file_read_mode_fake_data;
    bool got_V_ENABLE_SYSALL;
    bool got_V_ENABLE_SYSn;
    bool got_enable_min_point_sysn;
    bool got_min_point;
    bool got_min_point_fake_data;
    bool got_min_point_sysn_l[N_SYSTEMATICS];
    bool got_min_point_sysn_h[N_SYSTEMATICS];
    bool got_n_param_1;
    bool got_n_param_2;

    std::vector<TMarker *> mark_min_point_sysn_l;
    std::vector<TMarker *> mark_min_point_sysn_h;
    std::vector<TLine *> line_min_point_sysn_l;
    std::vector<TLine *> line_min_point_sysn_h;

    double min_fval;
    double min_fval_SYSALL;

};


#endif // NEWLOGLIKFITTER_MPS_DRAW_DATA_H
