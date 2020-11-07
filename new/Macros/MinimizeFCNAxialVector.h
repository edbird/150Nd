#ifndef MINIMIZEFCNAXIALVECTOR_H
#define MINIMIZEFCNAXIALVECTOR_H



// Root headers, minuit
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"


#define VECTOR_RANGE_CHECK 1
#define MEASURE_FUNCTION_CALL_TIME 0



//Double_t getNumberMC1D(const Int_t channel, const Int_t bin_ix, const std::vector<double> &p);
//Double_t getNumberMC2D(const Int_t channel, const Int_t bin_ix, const Int_t bin_iy, const std::vector<double> &p);



// TODO: may need to change precision to float not double using
// MnMachinePrecision::setPrecision

class MinimizeFCNAxialVector : public ROOT::Minuit2::FCNBase
{


    public:


    MinimizeFCNAxialVector()
        : error_def{1.0} // LL // TODO: check as I am actually minimizing chisquare
    {
    }

    virtual
    ~MinimizeFCNAxialVector()
    {
    }

    virtual
    double Up() const
    {
        return error_def;
    }

    virtual
    double operator()(const std::vector<double> &param) const
    {


        #if MEASURE_FUNCTION_CALL_TIME
        std::chrono::system_clock::time_point start_time = std::chrono::high_resolution_clock::now();
        #endif


        debuglevel = 1;
        


        nch = 0;

        /*
        std::cout << "operator()-> param: ";
        for(int pix = 0; pix < param.size(); ++ pix)
        {
            std::cout << param[pix] << " ";
        }
        std::cout << std::endl;
        */
       
        //Int_t nPar = param.size();

        // TODO: need to re-enable to use MODE fake data
        /*
        if(g_mode_fake_data == true)
        {
            if(allFakeDataSamples1D == nullptr)
            {
                if(debuglevel > 0)
                {
                    std::cout << "calling build_fake_data()" << std::endl;
                }
                std::cin.get();
                build_fake_data();
            }
        }
        */
        // TODO: manually call this from main function (loadFiles)


        // save the walk
        ll_walk.push_back(std::make_pair(param.at(1), param.at(0)));


        // draw the output
        //TString fname;
        //fname.Form("lliter_%d", counter);
        //draw_channel(1, p, std::string(fname));





        const int xi_31_ext_param_number = g_pg.get_xi_31_ext_param_number(); 


        if(debuglevel >= 2)
        {
            std::cout << std::scientific;
            std::cout << "logLikelihood"
                      << " param[0]=" << param[0]
                      << " param[" << xi_31_ext_param_number << "]="
                      << param[xi_31_ext_param_number] << std::endl;
        }


        /*
        if(systematic_energy_offset != systematic_energy_offset_last)
        {
            systematic_energy_offset_last = systematic_energy_offset;
            
            if(debuglevel >= 2)
            {
                std::cout << __func__ << " rebuild_150Nd_data() -> systematic_energy_offset=" << systematic_energy_offset << std::endl;
            }

            rebuild_150Nd_data();
        }
        */


        // TODO: will not work if parameter number changes
        // 2020-06-17
        //if(p[1] != last_xi_31_parameter_value)
        //if(p[1] != paramLastValueMap[1])
        //if(param[axial_vector_parameter_0_param_number] != paramLastValueMap[axial_vector_parameter_0_param_number])
        if(param[xi_31_ext_param_number] != g_pg.file_params.at(xi_31_ext_param_number).paramLastValue)
        {

            // TODO: rebuild nd150 xi_31 paramter histogram here
            //std::cout << "rebuilding 150 Nd MC" << std::endl;

            const double xi_31{param[xi_31_ext_param_number]};
            if(debuglevel >= 3)
            {
                std::cout << "rebuild 150Nd MC" << std::endl;
                std::cout << "param[" << xi_31_ext_param_number << "]=" << param[xi_31_ext_param_number] << " != " << g_pg.file_params.at(xi_31_ext_param_number).paramLastValue << std::endl;
            }


            // 2020-10-03:
            // TODO: this function does not take into account any systematics
            // is this correct?
            // yes this is correct, because the systematics are applied ONLY
            // when constructing the FAKE DATA from all the MC input samples
            // (read from TTrees)
            #if MEASURE_FUNCTION_CALL_TIME 
            std::chrono::system_clock::time_point start_time_rebuild_150Nd_MC = std::chrono::high_resolution_clock::now();
            #endif 
            rebuild_150Nd_MC(xi_31, xi_31_baseline);
            #if MEASURE_FUNCTION_CALL_TIME 
            std::chrono::system_clock::time_point end_time_rebuild_150Nd_MC = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> runtime_microsec_rebuild_150Nd_MC = end_time_rebuild_150Nd_MC - start_time_rebuild_150Nd_MC;
            std::cout << "rebuild_150Nd_MC, time=" << 1.0e+06 * runtime_microsec_rebuild_150Nd_MC.count() << " microsecond" << std::endl;
            #endif 


    /*
            TCanvas *ctmp = new TCanvas("ctmp", "ctmp");
            hSingleEnergy->Draw();
            TString fname;
            fname.Form("ctmp_%d.png", counter);
            ctmp->SaveAs(fname);
            ++ counter;
    */


            ///////////////////////////////////////////////////////////////////////
            // set last value
            ///////////////////////////////////////////////////////////////////////
            

            //last_xi_31_parameter_value = xi_31;
            // TODO: check that I set .paramLastValue in g_pg.file_params later at end of this function
            // if so can remove this statement
            g_pg.file_params.at(xi_31_ext_param_number).paramLastValue = xi_31;
        }
        else
        {
            if(debuglevel >= 4)
            {
                std::cout << "dont rebuild 150Nd MC" << std::endl;
            }
        }

        // TODO: add check here to see if any disabled parameters are accessed


        bool new_method = true;
        //bool new_method = false;
        double chi2_total = 0.0;
        if(new_method)
        {

            ///////////////////////////////////////////////////////////////////////
            // loglikelihood/chisquare function
            // new method using V Matrix
            ///////////////////////////////////////////////////////////////////////


    
            ///////////////////////////////////////////////////////////////////
            // CERN ROOT MathMore Matrix Lib Objects
            ///////////////////////////////////////////////////////////////////

        
            #if MEASURE_FUNCTION_CALL_TIME 
            std::chrono::system_clock::time_point start_time_2 = std::chrono::high_resolution_clock::now();
            #endif

            //if(recalculate_V_PHYS_MATHMORE == true)
            //{
                if(V_ENABLE_BIN_1D_P1[0] == nullptr)
                {
                    for(int ch = 0; ch < number1DHists; ++ ch)
                    {
                        //delete V_ENABLE_BIN_1D_P1[ch];
                        //delete V_ENABLE_BIN_1D_P2[ch];
                        V_ENABLE_BIN_1D_P1[ch] = new std::vector<bool>;
                        V_ENABLE_BIN_1D_P2[ch] = new std::vector<bool>;
                        const Int_t NUM_BINS_XY = 50;
                        V_ENABLE_BIN_1D_P1[ch]->reserve(NUM_BINS_XY);
                        V_ENABLE_BIN_1D_P2[ch]->reserve(NUM_BINS_XY); 
                    }
                }
            //}


            if(recalculate_V_PHYS_MATHMORE == true)
            {
                // realloc these each time
                for(int ch = 0; ch < number1DHists; ++ ch)
                {
                    if(channel_enable_1D[ch] == 1)
                    {
                        V_ENABLE_BIN_1D_P1[ch]->clear();
                        V_ENABLE_BIN_1D_P2[ch]->clear();
                    }
                    // elements are set using push_back(), so call clear
                    // in initialization
                }
            }


            //check_alloc_V_PHYS();
            //check_alloc_V_PHYS_data(); // not used
            //if(recalculate_V_PHYS_xD_Px_MATHMORE == true)
            //{
            //    //zero_V_PHYS();
            //    zero_V_PHYS_data(); // not used
            //}

            check_alloc_V_PHYS_STAT_data();
            set_V_PHYS_STAT_data(); // does not do anything, done elsewhere
            //if(recalculate_V_PHYS_xD_Px_MATHMORE == true)
            //if(true) // must be done each time? since M will change
            // however check function: set_V_MATRIX
            // elements set there so only need this on first run?
            // this will do for now 
            //if(recalculate_V_PHYS_xD_Px_MATHMORE == true)
            if(false)
            {
                // Set to zero

                // NOTE: 2020-09-16 this is a 2D object not a 1D object in a 
                // 2D data structure, hence require double loop
                for(int ch = 0; ch < number1DHists; ++ ch)
                {
                    /*
                    for(Int_t bj = 1; bj <= V_PHYS_STAT_1D_P1[ch]->GetNbinsY(); ++ bj)
                    {
                        for(Int_t bi = 1; bi <= V_PHYS_STAT_1D_P1[ch]->GetNbinsX(); ++ bi)
                        {
                            const Double_t zero = 0.0;
                            V_PHYS_STAT_1D_P1[ch]->SetBinContent(bi, bj, zero);
                        }
                    }
                    for(Int_t bj = 1; bj <= V_PHYS_STAT_1D_P1[ch]->GetNbinsY(); ++ bj)
                    {
                        for(Int_t bi = 1; bi <= V_PHYS_STAT_1D_P2[ch]->GetNbinsX(); ++ bi)
                        {
                            const Double_t zero = 0.0;
                            V_PHYS_STAT_1D_P2[ch]->SetBinContent(bi, bj, zero);
                        }
                    }
                    */
                    const Int_t NUM_BINS_XY = 50;
                    for(int i = 0; i < NUM_BINS_XY * NUM_BINS_XY; ++ i)
                    {
                        const double zero = 0.0;
                        #if VECTOR_RANGE_CHECK
                        V_PHYS_STAT_1D_P1_data[ch]->at(i) = zero;
                        V_PHYS_STAT_1D_P2_data[ch]->at(i) = zero;
                        #else
                        V_PHYS_STAT_1D_P1_data[ch]->operator[](i) = zero;
                        V_PHYS_STAT_1D_P2_data[ch]->operator[](i) = zero;
                        #endif
                    }
                }

                //  Set to zero
                /*
                for(int ch = 0; ch < number2DHists; ++ ch)
                {
                    for(Int_t bj = 1; bj <= V_PHYS_STAT_2D_P1[ch]->GetNbinsY(); ++ bj)
                    {
                        for(Int_t bi = 1; bi <= V_PHYS_STAT_2D_P1[ch]->GetNbinsX(); ++ bi)
                        {
                            const Double_t zero = 0.0;
                            V_PHYS_STAT_2D_P1[ch]->SetBinContent(bi, bj, zero);
                        }
                    }
                    for(Int_t bj = 1; bj <= V_PHYS_STAT_2D_P2[ch]->GetNbinsY(); ++ bj)
                    {
                        for(Int_t bi = 1; bi <= V_PHYS_STAT_2D_P2[ch]->GetNbinsX(); ++ bi)
                        {
                            const Double_t zero = 0.0;
                            V_PHYS_STAT_2D_P2[ch]->SetBinContent(bi, bj, zero);
                        }
                    }
                }
                */
            }


            for(int i = 0; i < N_SYSTEMATICS; ++ i)
            {
                // SYSn
                check_alloc_V_PHYS_SYSx_data_helper(V_PHYS_SYSn_1D_P1_data[i], V_PHYS_SYSn_1D_P2_data[i]);

                // SYSn
                set_V_PHYS_SYSx_data_helper(V_PHYS_SYSn_1D_P1_data[i],
                                            V_PHYS_SYSn_1D_P2_data[i],
                                            systematic_n_V_MATRIX_coeff_1D_P1[i],
                                            systematic_n_V_MATRIX_coeff_1D_P2[i]);
            }
           

            check_alloc_V_PHYS_SYSALL_data();
            set_V_PHYS_SYSALL_data();
            // TODO: move if statements from set functions to here


            #if MEASURE_FUNCTION_CALL_TIME 
            std::chrono::system_clock::time_point start_time_E = std::chrono::high_resolution_clock::now();
            #endif 
            check_alloc_D();
            #if MEASURE_FUNCTION_CALL_TIME 
            std::chrono::system_clock::time_point end_time_E = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> runtime_microsec_E = end_time_E - start_time_E;
            std::cout << "check_alloc_D(), time=" << 1.0e+06 * runtime_microsec_E.count() << " microsecond" << std::endl;
            #endif 



            #if MEASURE_FUNCTION_CALL_TIME 
            std::chrono::system_clock::time_point start_time_F = std::chrono::high_resolution_clock::now();
            #endif
            check_alloc_M();
            #if MEASURE_FUNCTION_CALL_TIME 
            std::chrono::system_clock::time_point end_time_F = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> runtime_microsec_F = end_time_F - start_time_F;
            std::cout << "check_alloc_M(), time=" << 1.0e+06 * runtime_microsec_F.count() << " microsecond" << std::endl;
            #endif



            #if MEASURE_FUNCTION_CALL_TIME 
            std::chrono::system_clock::time_point start_time_G = std::chrono::high_resolution_clock::now();
            #endif
            check_alloc_D_minus_M();
            #if MEASURE_FUNCTION_CALL_TIME 
            std::chrono::system_clock::time_point end_time_G = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> runtime_microsec_G = end_time_G - start_time_G;
            std::cout << "check_alloc_D_minus_M(), time=" << 1.0e+06 * runtime_microsec_G.count() << " microsecond" << std::endl;
            #endif



            #if MEASURE_FUNCTION_CALL_TIME
            std::chrono::system_clock::time_point end_time_2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> runtime_microsec_2 = end_time_2 - start_time_2;
            std::cout << "Done init, time=" << 1.0e+06 * runtime_microsec_2.count() << " microsecond" << std::endl;
            #endif



            ///////////////////////////////////////////////////////////////////////
            // set values for M, D and D_minus_M

            // can D change?
            // yes - if systematic parameters change and we are in fake data mode
            // just assume it always changes for now
            
            // set the contents of D
            set_D();
            set_M(param);
            // draw
            /*
            TCanvas *c_D = new TCanvas("c_D", "c_D", 10000, 10000);
            D->Draw("colz");
            c_D->SaveAs("debug_c_D.png");
            TCanvas *c_M = new TCanvas("c_M", "c_M", 10000, 10000);
            M->Draw("colz");
            c_M->SaveAs("debug_c_M.png");
            */

            set_D_minus_M();
            // draw
            /*
            TCanvas *c_D_minus_M = new TCanvas("c_D_minus_M", "c_D_minus_M", 10000, 10000);
            D_minus_M->Draw("colz");
            c_D_minus_M->SaveAs("debug_c_D_minus_M.png");
            */

            
            #if MEASURE_FUNCTION_CALL_TIME 
            std::chrono::system_clock::time_point start_time_H = std::chrono::high_resolution_clock::now();
            #endif
            set_V_MATRIX();
            #if MEASURE_FUNCTION_CALL_TIME 
            std::chrono::system_clock::time_point end_time_H = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> runtime_microsec_H = end_time_H - start_time_H;
            std::cout << "set_V_MATRIX(), time=" << 1.0e+06 * runtime_microsec_H.count() << " microsecond" << std::endl;
            #endif

            double chi2_P1 = 0.0;
            double chi2_P2 = 0.0;
            int nch_P1 = 0;
            int nch_P2 = 0; // overshadows previous definitions

            #if MEASURE_FUNCTION_CALL_TIME 
            std::chrono::system_clock::time_point start_time_I = std::chrono::high_resolution_clock::now();
            #endif
            calculate_chi2_P1(chi2_P1, nch_P1);
            #if MEASURE_FUNCTION_CALL_TIME 
            std::chrono::system_clock::time_point end_time_I = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> runtime_microsec_I = end_time_I - start_time_I;
            std::cout << "calculate_chi2_P1(), time=" << 1.0e+06 * runtime_microsec_I.count() << " microsecond" << std::endl;
            #endif

            #if MEASURE_FUNCTION_CALL_TIME 
            std::chrono::system_clock::time_point start_time_J = std::chrono::high_resolution_clock::now();
            #endif
            calculate_chi2_P2(chi2_P2, nch_P2);
            #if MEASURE_FUNCTION_CALL_TIME 
            std::chrono::system_clock::time_point end_time_J = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> runtime_microsec_J = end_time_J - start_time_J;
            std::cout << "calculate_chi2_P2(), time=" << 1.0e+06 * runtime_microsec_J.count() << " microsecond" << std::endl;
            #endif

            chi2_total = chi2_P1 + chi2_P2;
            nch = nch_P1 + nch_P2;
            //std::cout << "channel=" << channel << " P1 chi2=" << chi2_1D_P1 << std::endl;
            //std::cout << "channel=" << channel << " P2 chi2=" << chi2_1D_P2 << std::endl;
            
            //std::cout << "finished matrix calculations" << std::endl;
            
            // penalty terms section
            double penalty_term = 0.0;

            #if MEASURE_FUNCTION_CALL_TIME 
            std::chrono::system_clock::time_point start_time_K = std::chrono::high_resolution_clock::now();
            #endif
            calculate_penalty_term(penalty_term, param);
            #if MEASURE_FUNCTION_CALL_TIME 
            std::chrono::system_clock::time_point end_time_K = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> runtime_microsec_K = end_time_K - start_time_K;
            std::cout << "calculate_penalty_term(), time=" << 1.0e+06 * runtime_microsec_K.count() << " microsecond" << std::endl;
            #endif
            //std::cout << "penalty_term=" << penalty_term << std::endl;

            chi2_total += penalty_term;

            //std::cout << "chi2_total=" << chi2_total << std::endl;
            //std::cin.get();
            #if MEASURE_FUNCTION_CALL_TIME
            std::chrono::system_clock::time_point end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> runtime_microsec = end_time - start_time;
            std::cout << "operator() call time: " << 1.0e+06 * runtime_microsec.count() << " microsecond" << std::endl;
            #endif


            //recalculate_V_PHYS_MATHMORE = true; //?
            recalculate_V_PHYS_SYS = false;

            return chi2_total;

        }
        else // if not new_method
        {

            double fval = oldmethod(param);

            return fval;

        } // old method


    }


    protected:

    
    //void check_alloc_V_PHYS() const;
    //void zero_V_PHYS() const;
    
    void check_alloc_V_PHYS_data() const;
    
    void check_alloc_V_PHYS_STAT_data() const;
    //void check_alloc_V_PHYS_SYS1_data() const;
    //void check_alloc_V_PHYS_SYS2_data() const;
    //void check_alloc_V_PHYS_SYS3_data() const;
    //void check_alloc_V_PHYS_SYS4_data() const;
    void check_alloc_V_PHYS_SYSALL_data() const;
    
    void set_V_PHYS_STAT_data() const;
    //void set_V_PHYS_SYS1_data() const;
    //void set_V_PHYS_SYS2_data() const;
    //void set_V_PHYS_SYS3_data() const;
    //void set_V_PHYS_SYS4_data() const;
    void set_V_PHYS_SYSALL_data() const;


    ///////////////////////////////////////////////////////////////////////////
    // helper functions
    //
    void check_alloc_V_PHYS_SYSx_data_helper(
        std::vector<double>**,
        std::vector<double>**
    ) const;
    
    void set_V_PHYS_SYSx_data_helper(
        std::vector<double>**,
        std::vector<double>**,
        std::vector<double>**,
        std::vector<double>**
    ) const;
    //
    ///////////////////////////////////////////////////////////////////////////


    void check_alloc_D() const;
    void check_alloc_M() const;
    void check_alloc_D_minus_M() const;

    void zero_V_PHYS_data() const;
    
    void set_D() const;
    void set_M(const std::vector<double> &param) const;
    void set_D_minus_M() const;
    //void set_D_minus_M(const std::vector<double> &param) const;
    
    void set_V_MATRIX() const;
    
    void calculate_chi2_P1(double &chi2_P1, int &nch_P1) const;
    void calculate_chi2_P2(double &chi2_P2, int &nch_P2) const;
    void calculate_penalty_term(double &penalty_term, const std::vector<double> &param) const;

    double oldmethod(const std::vector<double> &param) const;

    mutable int debuglevel;

    double error_def;

    public:

    mutable int nch;

};


#include "MinimizeFCNAxialVector_oldmethod.h"
#include "MinimizeFCNAxialVector_check_alloc.h"
//#include "MinimizeFCNAxialVector_V_PHYS_STAT.h"
#include "MinimizeFCNAxialVector_V_PHYS_SYS.h"
#include "MinimizeFCNAxialVector_set_D.h"
#include "MinimizeFCNAxialVector_set_M.h"
#include "MinimizeFCNAxialVector_set_D_minus_M.h"
#include "MinimizeFCNAxialVector_set_V_MATRIX.h"
#include "MinimizeFCNAxialVector_calculate_chi2.h"
#include "MinimizeFCNAxialVector_calculate_penalty_term.h"



///////////////////////////////////////////////////////////////////////////////
// get number MC 1D
///////////////////////////////////////////////////////////////////////////////
#if 0
Double_t getNumberMC1D(
    const Int_t channel,
    const Int_t bin_ix,
    const std::vector<double> &p
    //const Double_t *const p
    )
{


    //std::cout << "getNumberMC1D() called with channel=" << channel << " binx=" << binx << std::endl;
    //std::cout << "printing contents of parameter array" << std::endl;
    //for(int i = 0; i < numberParams; ++ i)
    //{
    //    std::cout << "i=" << i << " p[i]=" << p[i] << std::endl;
    //}


    // std::cout <<"getting number of1D MC... "  <<channel << std::endl;

    double nMC = 0.;

    // (1) grab a hist from the sample list of this channel
    // (2) figure out which parameter it corresponds to
    TH1D *tmpHist;
    //int which_param;

    //std::cout << "getting number of MC... "  << channel << std::endl;

    //std::cout << "channel=" << channel << std::endl;
    //std::cout << allMCSamples1D[channel]->GetEntries() << std::endl;


    // copied from above
    //for(int k = 0; k < allMCSamples1D[channel]->GetEntries(); k++)
    for(int j = 0; j < allMCSamples1D[channel]->GetEntries(); j++)
    {

        // disabled params do not appear in this array

        tmpHist = (TH1D*)allMCSamples1D[channel]->At(j);

        TString tmpName = tmpHist->GetName();

        int which_param = -1;
        bool found_param = false;

        //std::cout << "NEW CODE" << std::endl;
        //try
        //{
        // TODO: remove TString
        {
            std::string tmp_hist_name(tmpName);
            auto i_start = tmp_hist_name.find('_') + 1;
            auto i_end = tmp_hist_name.rfind('_');
            if(i_end - i_start > 0)
            {
                std::string tmp_sample_name = tmp_hist_name.substr(i_start, i_end - i_start);
                //std::cout << "tmp_sample_name=" << tmp_sample_name << std::endl;
                int param_number = -1;
                bool success = g_pg.convert_MC_name_to_param_number(tmp_sample_name, param_number);
                //if(MCNameToParamNumberMap.count(tmp_sample_name) > 0)
                if(success == true)
                {
                    //int paramNumber = MCNameToParamNumberMap.at(tmp_sample_name);
                    // TODO: removed std::string, change tmpName type to be std::string from TString
                
                    //std::cout << "paramNumber=" << paramNumber << " -> tmp_sample_name=" << tmp_sample_name << " ~> tmpName=" << tmpName << std::endl;                    
                    //which_param = paramNumber;
                    //which_param = paramNumberToMinuitParamNumberMap.at(paramNumber);
                    // TODO: might still need the conversion from param number to minuit param number
                    found_param = true;

                    //std::cout << tmp_sample_name << " -> " << paramNumber << std::endl;

                    //std::cout << "DEBUG: found parameter with minuit (internal) number: " << which_param << std::endl;
                    //std::cin.get();
                }
                else
                {
                   std::cout << "ERROR: could not find " << tmp_sample_name << " in MCNameToParamNumberMap" << std::endl;
                }
            }
        }
        /*
        }
        catch(std::exception &e)
        {
            std::cout << "e.what(): " << e.what() << std::endl;
            std::cout << "tmpName=" << tmpName << std::endl;
            std::cout << "contents of map" << std::endl;
            for(auto it = MCNameToParamNumberMap.cbegin(); it != MCNameToParamNumberMap.cend(); ++ it)
            {
                std::cout << it->first << " -> " << it->second << std::endl;
            }
        }
        std::cin.get();
        */

        //std::cout << "bin_ix=" << bin_ix << " tmpHist->GetName()=" << tmpHist->GetName() << " which_param=" << which_param << std::endl;

        /*
        if(foundParam)
        {
            //std::cout << "adding to nMC with index of which_param=" << which_param << std::endl;
            // TODO: think this is collecting the wrong parameter? or is it?
            nMC += p[which_param] * tmpHist->GetBinContent(binx);
        }
        else
        {
            std::cout << "error could not find histogram: " << tmpName << std::endl;
        }
        */

        if(found_param == true)
        {

            // check here to see if param is disabled
            // TODO
            //
            // since which_param must be an enabled parameter, no longer need this
            // leave for now as a check?
            //if(std::find(enabled_params.begin(), enabled_params.end(), which_param) == enabled_params.end())
            /*
            if(std::find(enabled_params.begin(), enabled_params.end(), minuitParamNumberToParamNumberMap.at(which_param)) == enabled_params.end())
            {
                std::cout << "ERROR: which_param=" << which_param << " - parameter is DISABLED" << std::endl;
                std::cin.get();
            }
            */
            // TODO: enabled checks

            if(which_param == 1)
            {
                std::cout << "Error: which_param == 1 !" << std::endl;
                std::cin.get();
            }

            nMC += p[which_param] * (double)tmpHist->GetBinContent(bin_ix);

        }
        else
        {
            std::cout << "error could not find histogram: tmpName=" << tmpName << std::endl;
        } 

    }

    //std::cin.get();

    return nMC;

}
#endif


#if 0
Double_t getNumberMC2D(
    const Int_t channel,
    const Int_t bin_ix,
    const Int_t bin_iy,
    //const Double_t *const p
    const std::vector<double> &p
    )
{

    double nMC = 0.;

    // (1) grab a hist from the sample list of this channel
    // (2) figure out which parameter it corresponds to
    TH1D *tmpHist;
    //int which_param;

    //std::cout << "getting number of MC... "  << channel << std::endl;
    //std::cout << allMCSamples[channel]->GetEntries() << std::endl;

    // copied from above
    //for(int k = 0; k < allMCSamples1D[channel]->GetEntries(); k++)
    for(int j = 0; j < allMCSamples2D[channel]->GetEntries(); j++)
    {

        // get histogram name
        tmpHist = (TH1D*)allMCSamples2D[channel]->At(j);
        TString tmpName = tmpHist->GetName();

        // set paramter to default value
        int which_param = -1;
        bool found_param = false;

        {
            std::string tmp_hist_name(tmpName);
            auto i_start = tmp_hist_name.find('_') + 1;
            auto i_end = tmp_hist_name.rfind('_');
            if(i_end - i_start > 0)
            {
                // strip off WHAT - TODO comment
                std::string tmp_sample_name = tmp_hist_name.substr(i_start, i_end - i_start);
                //std::cout << "tmp_sample_name=" << tmp_sample_name << std::endl;
                //std::cin.get(); // check still working

                if(MCNameToParamNumberMap.count(tmp_sample_name) > 0)
                {
                    int paramNumber = MCNameToParamNumberMap.at(tmp_sample_name);
                
                    // convert to minuit (internal) parameter number
                    which_param = paramNumberToMinuitParamNumberMap.at(paramNumber);
                    found_param = true;

                    //std::cout << "DEBUG: found parameter with minuit (internal) number: " << which_param << std::endl;
                    //std::cin.get();
                }
                else
                {
                   std::cout << "ERROR: could not find " << tmp_sample_name << " in MCNameToParamNumberMap" << std::endl;
                }
            }
        }

        if(found_param == true)
        {

            // check here to see if param is disabled
            // TODO
            //
            // since which_param must be an enabled parameter, no longer need this
            // leave for now as a check?
            //if(std::find(enabled_params.begin(), enabled_params.end(), which_param) == enabled_params.end())
            if(std::find(enabled_params.begin(), enabled_params.end(), minuitParamNumberToParamNumberMap.at(which_param)) == enabled_params.end())
            {
                //std::cout << "parameter number " << param_number << " is disabled" << std::endl;
                std::cout << "ERROR: which_param=" << which_param << " - parameter is DISABLED" << std::endl;
                std::cin.get();
            }

            if(which_param == 1)
            {
                std::cout << "Error: which_param == 1 !" << std::endl;
                std::cin.get();
            }

            // add to number of MC
            nMC += p[which_param] * (double)tmpHist->GetBinContent(bin_ix, bin_iy);
        }
        else
        {
            std::cout << "error could not find histogram: tmpName=" << tmpName << std::endl;
        }

    } // j

    return nMC;


}
#endif





#endif // MINIMIZEFCNAXIALVECTOR_H
