#ifndef MINIMIZEFCNAXIALVECTOR_V_PHYS_SYS_H
#define MINIMIZEFCNAXIALVECTOR_V_PHYS_SYS_H


///////////////////////////////////////////////////////////////////////////////
// V_PHYS_SYSn
///////////////////////////////////////////////////////////////////////////////

void
MinimizeFCNAxialVector::check_alloc_V_PHYS_SYSx_data_helper(
    std::vector<double> *V_PHYS_SYSx_1D_P1_data[],
    std::vector<double> *V_PHYS_SYSx_1D_P2_data[]
    ) const
{
    //std::cout << __func__ << std::endl;

    if(V_PHYS_SYSx_1D_P1_data[0] == nullptr)
    {
        //std::cout << "Alloc V_PHYS_SYSx (P1)" << std::endl;

        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;
            V_PHYS_SYSx_1D_P1_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
        }
    }

    if(V_PHYS_SYSx_1D_P2_data[0] == nullptr)
    {
        //std::cout << "Alloc V_PHYS_SYSx (P2)" << std::endl;

        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;
            V_PHYS_SYSx_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
        }
    }

}


void
MinimizeFCNAxialVector::set_V_PHYS_SYSx_data_helper(
    std::vector<double> *V_PHYS_SYSx_1D_P1_data[],
    std::vector<double> *V_PHYS_SYSx_1D_P2_data[],
    std::vector<double> *systematic_X_V_MATRIX_coeff_1D_P1[],
    std::vector<double> *systematic_X_V_MATRIX_coeff_1D_P2[]
    ) const
{

    if(recalculate_V_PHYS_SYS == true)
    {

        //std::cout << __func__ << std::endl;

        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;

            // initialize elements of V_PHYS_SYS1_*
            int channel = ch;

            // TODO: symmetry optimization
            //for(Int_t biny{0}; biny < M_1D_P1_data[channel]->size(); ++ biny)
            for(Int_t biny{0}; biny < NUM_BINS_XY; ++ biny)
            {
                //for(Int_t binx{0}; binx < M_1D_P1_data[channel]->size(); ++ binx)
                //for(Int_t binx{0}; binx < M_1D_P1_data[channel]->size(); ++ binx)
                for(Int_t binx{0}; binx < NUM_BINS_XY; ++ binx)
                {

                    // P1
                    {
                        #if VECTOR_RANGE_CHECK
                        double coeff_x = systematic_X_V_MATRIX_coeff_1D_P1[channel]->at(binx);
                        double coeff_y = systematic_X_V_MATRIX_coeff_1D_P1[channel]->at(biny);
                        V_PHYS_SYSx_1D_P1_data[channel]->at(biny * 50 + binx) = coeff_x * coeff_y;
                        #else
                        double coeff_x = systematic_X_V_MATRIX_coeff_1D_P1[channel]->operator[](binx);
                        double coeff_y = systematic_X_V_MATRIX_coeff_1D_P1[channel]->operator[](biny);
                        V_PHYS_SYSx_1D_P1_data[channel]->operator[](biny * 50 + binx) = coeff_x * coeff_y;
                        #endif
                    }

                    // P2
                    {
                        #if VECTOR_RANGE_CHECK
                        double coeff_x = systematic_X_V_MATRIX_coeff_1D_P2[channel]->at(binx);
                        double coeff_y = systematic_X_V_MATRIX_coeff_1D_P2[channel]->at(biny);
                        V_PHYS_SYSx_1D_P2_data[channel]->at(biny * 50 + binx) = coeff_x * coeff_y;
                        #else
                        double coeff_x = systematic_X_V_MATRIX_coeff_1D_P2[channel]->operator[](binx);
                        double coeff_y = systematic_X_V_MATRIX_coeff_1D_P2[channel]->operator[](biny);
                        V_PHYS_SYSx_1D_P2_data[channel]->operator[](biny * 50 + binx) = coeff_x * coeff_y;
                        #endif
                    }


                } // binx
            } // biny
        } // channel
    }

}



///////////////////////////////////////////////////////////////////////////////
// V_PHYS_SYSALL
///////////////////////////////////////////////////////////////////////////////

void
MinimizeFCNAxialVector::check_alloc_V_PHYS_SYSALL_data() const
{
    //std::cout << __func__ << std::endl;


    if(V_PHYS_SYSALL_1D_P1_data[0] == nullptr)
    {
        //std::cout << "Alloc V_PHYS_SYSALL" << std::endl;

        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;
            V_PHYS_SYSALL_1D_P1_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            V_PHYS_SYSALL_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            //std::cout << "ALLOC " << __func__ << std::endl;
        }
    }

}


void
MinimizeFCNAxialVector::set_V_PHYS_SYSALL_data() const
{
    
    if(recalculate_V_PHYS_SYS == true)
    {

        //std::cout << __func__ << std::endl;
        /*
        std::cout << "V_ENABLE_SYSALL=" << V_ENABLE_SYSALL << std::endl;
        std::cout << "V_ENABLE_SYS1=" << V_ENABLE_SYS1 << std::endl;
        std::cout << "V_ENABLE_SYS2=" << V_ENABLE_SYS2 << std::endl;
        std::cout << "V_ENABLE_SYS3=" << V_ENABLE_SYS3 << std::endl;
        std::cout << "V_ENABLE_SYS4=" << V_ENABLE_SYS4 << std::endl;
        */

        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;

            // initialize elements of V_PHYS_SYS4_*
            int channel = ch;

            // TODO: symmetry optimization
            for(Int_t biny{0}; biny < NUM_BINS_XY; ++ biny)
            {
                for(Int_t binx{0}; binx < NUM_BINS_XY; ++ binx)
                {

                    // P1
                    {
                        double sysn[N_SYSTEMATICS];
                        for(int i = 0; i < N_SYSTEMATICS; ++ i)
                        {
                            sysn[i] = 0.0;
                        }
                        double sysall = 0.0;

                        #if VECTOR_RANGE_CHECK
                            for(int i = 0; i < N_SYSTEMATICS; ++ i)
                            {
                                if(V_ENABLE_SYSn[i] == true)
                                {
                                    sysn[i] = V_PHYS_SYSn_1D_P1_data[i][channel]->at(biny * 50 + binx);
                                }
                            }
                        #else
                            for(int i = 0; i < N_SYSTEMATICS; ++ i)
                            {
                                if(V_ENABLE_SYSn[i] == true)
                                {
                                    sysn[i] = V_PHYS_SYSn_1D_P1_data[i][channel]->operator[](biny * 50 + binx);
                                }
                            }
                        #endif

                        for(int i = 0; i < N_SYSTEMATICS; ++ i)
                        {
                            sysall += sysn[i];
                        }

                        #if VECTOR_RANGE_CHECK
                            V_PHYS_SYSALL_1D_P1_data[channel]->at(biny * 50 + binx) = sysall;
                        #else
                            V_PHYS_SYSALL_1D_P1_data[channel]->operator[](biny * 50 + binx) = sysall;
                        #endif
                    }

                    // P2
                    {
                        double sysn[N_SYSTEMATICS];
                        for(int i = 0; i < N_SYSTEMATICS; ++ i)
                        {
                            sysn[i] = 0.0;
                        }
                        double sysall = 0.0;

                        #if VECTOR_RANGE_CHECK
                            for(int i = 0; i < N_SYSTEMATICS; ++ i)
                            {
                                if(V_ENABLE_SYSn[i] == true)
                                {
                                    sysn[i] = V_PHYS_SYSn_1D_P2_data[i][channel]->at(biny * 50 + binx);
                                }
                            }
                        #else
                            for(int i = 0; i < N_SYSTEMATICS; ++ i)
                            {
                                if(V_ENABLE_SYSn[i] == true)
                                {
                                    sysn[i] = V_PHYS_SYSn_1D_P2_data[i][channel]->operator[](biny * 50 + binx);
                                }
                            }
                        #endif
                            
                        for(int i = 0; i < N_SYSTEMATICS; ++ i)
                        {
                            sysall += sysn[i];
                        }

                        #if VECTOR_RANGE_CHECK
                            V_PHYS_SYSALL_1D_P2_data[channel]->at(biny * 50 + binx) = sysall;
                        #else
                            V_PHYS_SYSALL_1D_P2_data[channel]->operator[](biny * 50 + binx) = sysall;
                        #endif
                    }


                }
            }
        }

       
    }

}


#endif // MINIMIZEFCNAXIALVECTOR_V_PHYS_SYS_H
