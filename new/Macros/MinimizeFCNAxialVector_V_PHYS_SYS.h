#ifndef MINIMIZEFCNAXIALVECTOR_V_PHYS_SYS_H
#define MINIMIZEFCNAXIALVECTOR_V_PHYS_SYS_H



///////////////////////////////////////////////////////////////////////////////
// V_PHYS_SYS1
///////////////////////////////////////////////////////////////////////////////

void
MinimizeFCNAxialVector::check_alloc_V_PHYS_SYS1_data() const
{
    //std::cout << __func__ << std::endl;

    if(V_PHYS_SYS1_1D_P1_data[0] == nullptr)
    {
        std::cout << "Alloc V_PHYS_SYS1" << std::endl;

        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;
            V_PHYS_SYS1_1D_P1_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            V_PHYS_SYS1_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            std::cout << "ALLOC " << __func__ << std::endl;
        }
    }

}


void
MinimizeFCNAxialVector::set_V_PHYS_SYS1_data() const
{

    if(recalculate_V_PHYS_SYS == true)
    {

        std::cout << __func__ << std::endl;

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
                        double coeff_x = systematic_offset_V_MATRIX_coeff_1D_P1[channel]->at(binx);
                        double coeff_y = systematic_offset_V_MATRIX_coeff_1D_P1[channel]->at(biny);
                        V_PHYS_SYS1_1D_P1_data[channel]->at(biny * 50 + binx) = coeff_x * coeff_y;
                        #else
                        double coeff_x = systematic_offset_V_MATRIX_coeff_1D_P1[channel]->operator[](binx);
                        double coeff_y = systematic_offset_V_MATRIX_coeff_1D_P1[channel]->operator[](biny);
                        V_PHYS_SYS1_1D_P1_data[channel]->operator[](biny * 50 + binx) = coeff_x * coeff_y; // TODO: abs value here?
                        //V_PHYS_SYS1_1D_P1_data[channel]->operator[](biny * 50 + binx) = std::abs(coeff1) * std::abs(coeff2); // TODO: abs value here?
                        // NOTE: does not seem to work: no difference
                        #endif
                    }

                    // P2
                    {
                        #if VECTOR_RANGE_CHECK
                        double coeff_x = systematic_offset_V_MATRIX_coeff_1D_P2[channel]->at(binx);
                        double coeff_y = systematic_offset_V_MATRIX_coeff_1D_P2[channel]->at(biny);
                        V_PHYS_SYS1_1D_P2_data[channel]->at(biny * 50 + binx) = coeff_x * coeff_y;
                        #else
                        double coeff_x = systematic_offset_V_MATRIX_coeff_1D_P2[channel]->operator[](binx);
                        double coeff_y = systematic_offset_V_MATRIX_coeff_1D_P2[channel]->operator[](biny);
                        V_PHYS_SYS1_1D_P2_data[channel]->operator[](biny * 50 + binx) = coeff_x * coeff_y; // TODO: abs value here?
                        //V_PHYS_SYS1_1D_P2_data[channel]->operator[](biny * 50 + binx) = std::abs(coeff1) * std::abs(coeff2); // TODO: abs value here?
                        #endif
                    }


                } // binx
            } // biny
        } // channel
    }

}



///////////////////////////////////////////////////////////////////////////////
// V_PHYS_SYS2
///////////////////////////////////////////////////////////////////////////////

void
MinimizeFCNAxialVector::check_alloc_V_PHYS_SYS2_data() const
{
    //std::cout << __func__ << std::endl;


    if(V_PHYS_SYS2_1D_P1_data[0] == nullptr)
    {
        std::cout << "Alloc V_PHYS_SYS2" << std::endl;

        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;
            V_PHYS_SYS2_1D_P1_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            V_PHYS_SYS2_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            std::cout << "ALLOC " << __func__ << std::endl;
        }
    }

}


void
MinimizeFCNAxialVector::set_V_PHYS_SYS2_data() const
{
    
    if(recalculate_V_PHYS_SYS == true)
    {

        std::cout << __func__ << std::endl;

        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;

            // initialize elements of V_PHYS_SYS2_*
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
                        double coeff1 = systematic_scale_V_MATRIX_coeff_1D_P1[channel]->at(binx);
                        double coeff2 = systematic_scale_V_MATRIX_coeff_1D_P1[channel]->at(biny);
                        V_PHYS_SYS2_1D_P1_data[channel]->at(biny * 50 + binx) = coeff1 * coeff2;
                        #else
                        double coeff1 = systematic_scale_V_MATRIX_coeff_1D_P1[channel]->operator[](binx);
                        double coeff2 = systematic_scale_V_MATRIX_coeff_1D_P1[channel]->operator[](biny);
                        V_PHYS_SYS2_1D_P1_data[channel]->operator[](biny * 50 + binx) = coeff1 * coeff2;
                        #endif

                        //std::cout << "biny=" << biny << " binx=" << binx << " " << coeff1 << " * " << coeff2 << std::endl;
                    }

                    // P2
                    {
                        #if VECTOR_RANGE_CHECK
                        double coeff1 = systematic_scale_V_MATRIX_coeff_1D_P2[channel]->at(binx);
                        double coeff2 = systematic_scale_V_MATRIX_coeff_1D_P2[channel]->at(biny);
                        V_PHYS_SYS2_1D_P2_data[channel]->at(biny * 50 + binx) = coeff1 * coeff2;
                        #else
                        double coeff1 = systematic_scale_V_MATRIX_coeff_1D_P2[channel]->operator[](binx);
                        double coeff2 = systematic_scale_V_MATRIX_coeff_1D_P2[channel]->operator[](biny);
                        V_PHYS_SYS2_1D_P2_data[channel]->operator[](biny * 50 + binx) = coeff1 * coeff2;
                        #endif

                        //std::cout << "biny=" << biny << " binx=" << binx << " " << coeff1 << " * " << coeff2 << std::endl;
                        //std::cin.get();
                    }


                }
            }
        }

       
    }

}



///////////////////////////////////////////////////////////////////////////////
// V_PHYS_SYS3
///////////////////////////////////////////////////////////////////////////////

void
MinimizeFCNAxialVector::check_alloc_V_PHYS_SYS3_data() const
{
    //std::cout << __func__ << std::endl;


    if(V_PHYS_SYS3_1D_P1_data[0] == nullptr)
    {
        std::cout << "Alloc V_PHYS_SYS3" << std::endl;

        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;
            V_PHYS_SYS3_1D_P1_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            V_PHYS_SYS3_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            std::cout << "ALLOC " << __func__ << std::endl;
        }
    }

}


void
MinimizeFCNAxialVector::set_V_PHYS_SYS3_data() const
{
    
    if(recalculate_V_PHYS_SYS == true)
    {

        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;

            // initialize elements of V_PHYS_SYS3_*
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
                        double coeff1 = systematic_efficiency_V_MATRIX_coeff_1D_P1[channel]->at(binx);
                        double coeff2 = systematic_efficiency_V_MATRIX_coeff_1D_P1[channel]->at(biny);
                        V_PHYS_SYS3_1D_P1_data[channel]->at(biny * 50 + binx) = coeff1 * coeff2;
                        #else
                        double coeff1 = systematic_efficiency_V_MATRIX_coeff_1D_P1[channel]->operator[](binx);
                        double coeff2 = systematic_efficiency_V_MATRIX_coeff_1D_P1[channel]->operator[](biny);
                        V_PHYS_SYS3_1D_P1_data[channel]->operator[](biny * 50 + binx) = coeff1 * coeff2;
                        #endif
                    }

                    // P2
                    {
                        #if VECTOR_RANGE_CHECK
                        double coeff1 = systematic_efficiency_V_MATRIX_coeff_1D_P2[channel]->at(binx);
                        double coeff2 = systematic_efficiency_V_MATRIX_coeff_1D_P2[channel]->at(biny);
                        V_PHYS_SYS3_1D_P2_data[channel]->at(biny * 50 + binx) = coeff1 * coeff2;
                        #else
                        double coeff1 = systematic_efficiency_V_MATRIX_coeff_1D_P2[channel]->operator[](binx);
                        double coeff2 = systematic_efficiency_V_MATRIX_coeff_1D_P2[channel]->operator[](biny);
                        V_PHYS_SYS3_1D_P2_data[channel]->operator[](biny * 50 + binx) = coeff1 * coeff2;
                        #endif
                    }


                }
            }
        }

       
    }

}



///////////////////////////////////////////////////////////////////////////////
// V_PHYS_SYS4
///////////////////////////////////////////////////////////////////////////////

void
MinimizeFCNAxialVector::check_alloc_V_PHYS_SYS4_data() const
{
    //std::cout << __func__ << std::endl;


    if(V_PHYS_SYS4_1D_P1_data[0] == nullptr)
    {
        std::cout << "Alloc V_PHYS_SYS4" << std::endl;

        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;
            V_PHYS_SYS4_1D_P1_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            V_PHYS_SYS4_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            std::cout << "ALLOC " << __func__ << std::endl;
        }
    }

}


void
MinimizeFCNAxialVector::set_V_PHYS_SYS4_data() const
{

    if(recalculate_V_PHYS_SYS == true)
    {

        std::cout << __func__ << std::endl;

        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;

            // initialize elements of V_PHYS_SYS4_*
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
                        double coeff1 = systematic_enrichment_V_MATRIX_coeff_1D_P1[channel]->at(binx);
                        double coeff2 = systematic_enrichment_V_MATRIX_coeff_1D_P1[channel]->at(biny);
                        V_PHYS_SYS4_1D_P1_data[channel]->at(biny * 50 + binx) = coeff1 * coeff2;
                        #else
                        double coeff1 = systematic_enrichment_V_MATRIX_coeff_1D_P1[channel]->operator[](binx);
                        double coeff2 = systematic_enrichment_V_MATRIX_coeff_1D_P1[channel]->operator[](biny);
                        V_PHYS_SYS4_1D_P1_data[channel]->operator[](biny * 50 + binx) = coeff1 * coeff2;
                        #endif

                        // these bins are zero - the 150Nd does not extend to
                        // bin 34, only the BKGs do
                        /*
                        if((biny == 34) || (binx == 34))
                        {
                            std::cout << coeff1 << " " << coeff2 << std::endl;
                            std::cin.get();
                        }
                        */

                    }

                    // P2
                    {
                        #if VECTOR_RANGE_CHECK
                        double coeff1 = systematic_enrichment_V_MATRIX_coeff_1D_P2[channel]->at(binx);
                        double coeff2 = systematic_enrichment_V_MATRIX_coeff_1D_P2[channel]->at(biny);
                        V_PHYS_SYS4_1D_P2_data[channel]->at(biny * 50 + binx) = coeff1 * coeff2;
                        #else
                        double coeff1 = systematic_enrichment_V_MATRIX_coeff_1D_P2[channel]->operator[](binx);
                        double coeff2 = systematic_enrichment_V_MATRIX_coeff_1D_P2[channel]->operator[](biny);
                        V_PHYS_SYS4_1D_P2_data[channel]->operator[](biny * 50 + binx) = coeff1 * coeff2;
                        #endif

                        // these bins are zero - the 150Nd does not extend to
                        // bin 34, only the BKGs do
                        /*
                        if((biny == 34) || (binx == 34))
                        {
                            std::cout << "biny=" << biny << " binx=" << binx << " " << coeff1 << " " << coeff2 << std::endl;
                            std::cin.get();
                        }
                        */

                    }


                }
            }
        }

       
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
        std::cout << "Alloc V_PHYS_SYSALL" << std::endl;

        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;
            V_PHYS_SYSALL_1D_P1_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            V_PHYS_SYSALL_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            std::cout << "ALLOC " << __func__ << std::endl;
        }
    }

}


void
MinimizeFCNAxialVector::set_V_PHYS_SYSALL_data() const
{
    
    if(recalculate_V_PHYS_SYS == true)
    {

        std::cout << __func__ << std::endl;
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
                        double sys1 = 0.0;
                        double sys2 = 0.0;
                        double sys3 = 0.0;
                        double sys4 = 0.0;
                        double sysall = 0.0;
                        #if VECTOR_RANGE_CHECK
                            if(V_ENABLE_SYS1 == true)
                            {
                                sys1 = V_PHYS_SYS1_1D_P1_data[channel]->at(biny * 50 + binx);
                            }
                            if(V_ENABLE_SYS2 == true)
                            {
                                sys2 = V_PHYS_SYS2_1D_P1_data[channel]->at(biny * 50 + binx);
                            }
                            if(V_ENABLE_SYS3 == true)
                            {
                                sys3 = V_PHYS_SYS3_1D_P1_data[channel]->at(biny * 50 + binx);
                            }
                            if(V_ENABLE_SYS4 == true)
                            {
                                sys4 = V_PHYS_SYS4_1D_P1_data[channel]->at(biny * 50 + binx);
                            }

                            sysall = sys1 + sys2 + sys3 + sys4;
                            V_PHYS_SYSALL_1D_P1_data[channel]->at(biny * 50 + binx) = sysall;
                        #else
                            if(V_ENABLE_SYS1 == true)
                            {
                                sys1 = V_PHYS_SYS1_1D_P1_data[channel]->operator[](biny * 50 + binx);
                            }
                            if(V_ENABLE_SYS2 == true)
                            {
                                sys2 = V_PHYS_SYS2_1D_P1_data[channel]->operator[](biny * 50 + binx);
                            }
                            if(V_ENABLE_SYS3 == true)
                            {
                                sys3 = V_PHYS_SYS3_1D_P1_data[channel]->operator[](biny * 50 + binx);
                            }
                            if(V_ENABLE_SYS4 == true)
                            {
                                sys4 = V_PHYS_SYS4_1D_P1_data[channel]->operator[](biny * 50 + binx);
                            }

                            sysall = sys1 + sys2 + sys3 + sys4;
                            V_PHYS_SYSALL_1D_P1_data[channel]->operator[](biny * 50 + binx) = sysall;
                        #endif

                    }

                    // P2
                    {
                        double sys1 = 0.0;
                        double sys2 = 0.0;
                        double sys3 = 0.0;
                        double sys4 = 0.0;
                        double sysall = 0.0;
                        #if VECTOR_RANGE_CHECK
                            if(V_ENABLE_SYS1 == true)
                            {
                                sys1 = V_PHYS_SYS1_1D_P2_data[channel]->at(biny * 50 + binx);
                            }
                            if(V_ENABLE_SYS2 == true)
                            {
                                sys2 = V_PHYS_SYS2_1D_P2_data[channel]->at(biny * 50 + binx);
                            }
                            if(V_ENABLE_SYS3 == true)
                            {
                                sys3 = V_PHYS_SYS3_1D_P2_data[channel]->at(biny * 50 + binx);
                            }
                            if(V_ENABLE_SYS4 == true)
                            {
                                sys4 = V_PHYS_SYS4_1D_P2_data[channel]->at(biny * 50 + binx);
                            }

                            sysall = sys1 + sys2 + sys3 + sys4;
                            V_PHYS_SYSALL_1D_P2_data[channel]->at(biny * 50 + binx) = sysall;
                        #else
                            if(V_ENABLE_SYS1 == true)
                            {
                                sys1 = V_PHYS_SYS1_1D_P2_data[channel]->operator[](biny * 50 + binx);
                            }
                            if(V_ENABLE_SYS2 == true)
                            {
                                sys2 = V_PHYS_SYS2_1D_P2_data[channel]->operator[](biny * 50 + binx);
                            }
                            if(V_ENABLE_SYS3 == true)
                            {
                                sys3 = V_PHYS_SYS3_1D_P2_data[channel]->operator[](biny * 50 + binx);
                            }
                            if(V_ENABLE_SYS4 == true)
                            {
                                sys4 = V_PHYS_SYS4_1D_P2_data[channel]->operator[](biny * 50 + binx);
                            }

                            sysall = sys1 + sys2 + sys3 + sys4; 
                            V_PHYS_SYSALL_1D_P2_data[channel]->operator[](biny * 50 + binx) = sysall;
                        #endif
                    }


                }
            }
        }

       
    }

}


#endif // MINIMIZEFCNAXIALVECTOR_V_PHYS_SYS_H
