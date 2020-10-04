#ifndef MINIMIZEFCNAXIALVECTOR_CHECK_ALLOC_H
#define MINIMIZEFCNAXIALVECTOR_CHECK_ALLOC_H

void
MinimizeFCNAxialVector::check_alloc_V_PHYS_data() const
{

    ///////////////////////////////////////////////////////////////////
    // V_PHYS
    ///////////////////////////////////////////////////////////////////

    if(V_PHYS_1D_P1_data[0] == nullptr)
    {
        std::cout << "Alloc V_PHYS" << std::endl;

        // 1D
        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;
            V_PHYS_1D_P1_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            V_PHYS_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            std::cout << "ALLOC " << __func__ << std::endl;
        }

        // 2D
        /*
        for(int ch = 0; ch < number2DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50 * 50;
            TString hname;

            hname.Form("V_PHYS_2D_P1_CH%d", ch);
            V_PHYS_2D_P1[ch] = new TH2D(hname, hname,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY);

            hname.Form("V_PHYS_2D_P2_CH%d", ch);
            V_PHYS_2D_P2[ch] = new TH2D(hname, hname,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY);
        }
        */
    }

}


void
MinimizeFCNAxialVector::check_alloc_V_PHYS_STAT_data() const
{

    ///////////////////////////////////////////////////////////////////
    // V_PHYS_STAT
    ///////////////////////////////////////////////////////////////////
    //if(V_PHYS_STAT_1D_P1[0] == nullptr)
    
    if(V_PHYS_STAT_1D_P1_data[0] == nullptr)
    {
        //std::cout << "Alloc V_PHYS_STAT" << std::endl;
        //V_PHYS_STAT = new TH2D("V_PHYS_STAT", "V_PHYS_STAT",
        //                  NUM_BINS_XY, 0.0, NUM_BINS_XY,
        //                  NUM_BINS_XY, 0.0, NUM_BINS_XY);
        std::cout << "Alloc V_PHYS_STAT" << std::endl;
        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;
            V_PHYS_STAT_1D_P1_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            V_PHYS_STAT_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            std::cout << "ALLOC " << __func__ << std::endl;
        }

        /*
        for(int ch = 0; ch < number2DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50 * 50;
            TString hname;

            hname.Form("V_PHYS_STAT_2D_P1_CH%d", ch);
            //V_PHYS_STAT_2D_P1[ch] = new TH2D(hname, hname,
            //                  NUM_BINS_XY, 0.0, NUM_BINS_XY,
            //                  NUM_BINS_XY, 0.0, NUM_BINS_XY);

            hname.Form("V_PHYS_STAT_2D_P2_CH%d", ch);
            //V_PHYS_STAT_2D_P2[ch] = new TH2D(hname, hname,
            //                  NUM_BINS_XY, 0.0, NUM_BINS_XY,
            //                  NUM_BINS_XY, 0.0, NUM_BINS_XY);

            V_PHYS_STAT_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
        }
        */
    }

}


void
MinimizeFCNAxialVector::check_alloc_V_PHYS_SYS1_data() const
{

    ///////////////////////////////////////////////////////////////////
    // V_PHYS_SYS1
    ///////////////////////////////////////////////////////////////////

    //if(V_PHYS_SYS1_1D_P1[0] == nullptr)
    if(V_PHYS_SYS1_1D_P1_data[0] == nullptr)
    {
        //std::cout << "Alloc V_PHYS_SYS" << std::endl;
        //V_PHYS_SYS = new TH2D("V_PHYS_SYS", "V_PHYS_SYS",
        //                  NUM_BINS_XY, 0.0, NUM_BINS_XY,
        //                  NUM_BINS_XY, 0.0, NUM_BINS_XY);

        // clear the contents of V_PHYS_STAT
        //for(Int_t bin_i{1}; bin_i <= V_PHYS_STAT->GetNbinsX(); ++ bin_i)
        //{
        //    for(Int_t bin_j{1}; bin_j <= V_PHYS_STAT->GetNbinsY(); ++ bin_j)
        //    {
        //        const Double_t zero = 0.0;
        //        V_PHYS_STAT->SetBinContent(bin_i, bin_j, zero);
        //    }
        //}
        std::cout << "Alloc V_PHYS_SYS1" << std::endl;
        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;
            V_PHYS_SYS1_1D_P1_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            V_PHYS_SYS1_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            std::cout << "ALLOC " << __func__ << std::endl;
        

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
                        double coeff1 = systematic_offset_V_MATRIX_coeff_1D_P1[channel]->at(binx);
                        double coeff2 = systematic_offset_V_MATRIX_coeff_1D_P1[channel]->at(biny);
                        V_PHYS_SYS1_1D_P1_data[channel]->at(biny * 50 + binx) = coeff1 * coeff2;
                        #else
                        double coeff1 = systematic_offset_V_MATRIX_coeff_1D_P1[channel]->operator[](binx);
                        double coeff2 = systematic_offset_V_MATRIX_coeff_1D_P1[channel]->operator[](biny);
                        V_PHYS_SYS1_1D_P1_data[channel]->operator[](biny * 50 + binx) = coeff1 * coeff2; // TODO: abs value here?
                        //V_PHYS_SYS1_1D_P1_data[channel]->operator[](biny * 50 + binx) = std::abs(coeff1) * std::abs(coeff2); // TODO: abs value here?
                        // NOTE: does not seem to work: no difference
                        #endif
                    }

                    // P2
                    {
                        #if VECTOR_RANGE_CHECK
                        double coeff1 = systematic_offset_V_MATRIX_coeff_1D_P2[channel]->at(binx);
                        double coeff2 = systematic_offset_V_MATRIX_coeff_1D_P2[channel]->at(biny);
                        V_PHYS_SYS1_1D_P2_data[channel]->at(biny * 50 + binx) = coeff1 * coeff2;
                        #else
                        double coeff1 = systematic_offset_V_MATRIX_coeff_1D_P2[channel]->operator[](binx);
                        double coeff2 = systematic_offset_V_MATRIX_coeff_1D_P2[channel]->operator[](biny);
                        V_PHYS_SYS1_1D_P2_data[channel]->operator[](biny * 50 + binx) = coeff1 * coeff2; // TODO: abs value here?
                        //V_PHYS_SYS1_1D_P2_data[channel]->operator[](biny * 50 + binx) = std::abs(coeff1) * std::abs(coeff2); // TODO: abs value here?
                        #endif
                    }


                }
            }
        }



        /*
        for(int ch = 0; ch < number2DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50 * 50;
            TString hname;

            hname.Form("V_PHYS_SYS1_2D_P1_CH%d", ch);
            V_PHYS_SYS1_2D_P1[ch] = new TH2D(hname, hname,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY);

            hname.Form("V_PHYS_SYS1_2D_P2_CH%d", ch);
            V_PHYS_SYS1_2D_P2[ch] = new TH2D(hname, hname,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY);
        }
        */
    }

}


void
MinimizeFCNAxialVector::check_alloc_V_PHYS_SYS2_data() const
{

    ///////////////////////////////////////////////////////////////////
    // V_PHYS_SYS2
    ///////////////////////////////////////////////////////////////////

    if(V_PHYS_SYS2_1D_P1_data[0] == nullptr)
    {
        //std::cout << "Alloc V_PHYS_SYS" << std::endl;
        //V_PHYS_SYS = new TH2D("V_PHYS_SYS", "V_PHYS_SYS",
        //                  NUM_BINS_XY, 0.0, NUM_BINS_XY,
        //                  NUM_BINS_XY, 0.0, NUM_BINS_XY);

        // clear the contents of V_PHYS_STAT
        //for(Int_t bin_i{1}; bin_i <= V_PHYS_STAT->GetNbinsX(); ++ bin_i)
        //{
        //    for(Int_t bin_j{1}; bin_j <= V_PHYS_STAT->GetNbinsY(); ++ bin_j)
        //    {
        //        const Double_t zero = 0.0;
        //        V_PHYS_STAT->SetBinContent(bin_i, bin_j, zero);
        //    }
        //}
        std::cout << "Alloc V_PHYS_SYS2" << std::endl;
        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;
            V_PHYS_SYS2_1D_P1_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            V_PHYS_SYS2_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            std::cout << "ALLOC " << __func__ << std::endl;
        

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
                    }


                }
            }
        }

       
    }

}


void
MinimizeFCNAxialVector::check_alloc_V_PHYS_SYS3_data() const
{

    ///////////////////////////////////////////////////////////////////
    // V_PHYS_SYS3
    ///////////////////////////////////////////////////////////////////

    if(V_PHYS_SYS3_1D_P1_data[0] == nullptr)
    {
        std::cout << "Alloc V_PHYS_SYS3" << std::endl;
        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;
            V_PHYS_SYS3_1D_P1_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            V_PHYS_SYS3_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            std::cout << "ALLOC " << __func__ << std::endl;

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


void
MinimizeFCNAxialVector::check_alloc_V_PHYS_SYS4_data() const
{

    ///////////////////////////////////////////////////////////////////
    // V_PHYS_SYS4
    ///////////////////////////////////////////////////////////////////

    if(V_PHYS_SYS4_1D_P1_data[0] == nullptr)
    {
        std::cout << "Alloc V_PHYS_SYS4" << std::endl;
        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;
            V_PHYS_SYS4_1D_P1_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            V_PHYS_SYS4_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            std::cout << "ALLOC " << __func__ << std::endl;

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
                    }


                }
            }
        }

       
    }

}


void MinimizeFCNAxialVector::check_alloc_D() const
{

    ///////////////////////////////////////////////////////////////////
    // V_MATRIX D (data)
    ///////////////////////////////////////////////////////////////////

    //if(D_1D_P1[0] == nullptr)
    if(D_1D_P1_data[0] == nullptr)
    {
        //std::cout << "Alloc D" << std::endl;
        //D = new TH2D("D", "D",
        //                  NUM_BINS_XY, 0.0, NUM_BINS_XY,
        //                  1, 0.0, 1.0);

        // clear the contents of D
        //for(Int_t bin_i{1}; bin_i <= D->GetNbinsX(); ++ bin_i)
        //{
        //    const Double_t zero = 0.0;
        //    D->SetBinContent(bin_i, 1, zero);
        //}
        std::cout << "Alloc V_D" << std::endl;
        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;
            D_1D_P1_data[ch] = new std::vector<double>(NUM_BINS_XY, 0.0);
            D_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY, 0.0);
            std::cout << "ALLOC " << __func__ << std::endl;
        }

        /*
        for(int ch = 0; ch < number2DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50 * 50;
            TString hname;

            hname.Form("D_2D_P1_CH%d", ch);
            D_2D_P1[ch] = new TH2D(hname, hname,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY);

            hname.Form("D_2D_P2_CH%d", ch);
            D_2D_P2[ch] = new TH2D(hname, hname,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY);
        }
        */
    }

}



void MinimizeFCNAxialVector::check_alloc_M() const
{

    ///////////////////////////////////////////////////////////////////
    // V_MATRIX M (MC)
    ///////////////////////////////////////////////////////////////////

    //if(M_1D_P1[0] == nullptr)
    if(M_1D_P1_data[0] == nullptr)
    {
        //std::cout << "Alloc M" << std::endl;
        //M = new TH2D("M", "M",
        //                  NUM_BINS_XY, 0.0, NUM_BINS_XY,
        //                  1, 0.0, 1.0);

        // clear the contents of M
        //for(Int_t bin_i{1}; bin_i <= M->GetNbinsX(); ++ bin_i)
        //{
        //    const Double_t zero = 0.0;
        //    M->SetBinContent(bin_i, 1, zero);
        //}

        // 1D
        std::cout << "Alloc V_M" << std::endl;
        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;
            M_1D_P1_data[ch] = new std::vector<double>(NUM_BINS_XY, 0.0);
            M_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY, 0.0);
            std::cout << "ALLOC " << __func__ << std::endl;
        }

        // 2D
        /*
        for(int ch = 0; ch < number2DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50 * 50;
            TString hname;

            hname.Form("M_2D_P1_CH%d", ch);
            M_2D_P1[ch] = new TH2D(hname, hname,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY);

            hname.Form("M_2D_P2_CH%d", ch);
            M_2D_P2[ch] = new TH2D(hname, hname,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY);
        }
        */
    }

}



void MinimizeFCNAxialVector::check_alloc_D_minus_M() const
{

    ///////////////////////////////////////////////////////////////////
    // V_MATRIX D_minus_M (Data - MC)
    ///////////////////////////////////////////////////////////////////

    //if(D_minus_M_1D_P1[0] == nullptr)
    if(D_minus_M_1D_P1_data[0] == nullptr)
    {
        //std::cout << "Alloc D_minus_M" << std::endl;
        //D_minus_M = new TH2D("D_minus_M", "D_minus_M",
        //                  NUM_BINS_XY, 0.0, NUM_BINS_XY,
        //                  1, 0.0, 1.0);

        // clear the contents of M
        //for(Int_t bin_i{1}; bin_i <= D_minus_M->GetNbinsX(); ++ bin_i)
        //{
        //    const Double_t zero = 0.0;
        //    D_minus_M->SetBinContent(bin_i, 1, zero);
        //}

        // 1D
        std::cout << "Alloc D_minus_M" << std::endl;
        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;
            D_minus_M_1D_P1_data[ch] = new std::vector<double>(NUM_BINS_XY, 0.0);
            D_minus_M_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY, 0.0);
            std::cout << "ALLOC " << __func__ << std::endl;
        }


        // 2D
        /*
        for(int ch = 0; ch < number2DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50 * 50;
            TString hname;

            hname.Form("D_minus_M_2D_P1_CH%d", ch);
            D_minus_M_2D_P1[ch] = new TH2D(hname, hname,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY);

            hname.Form("D_minus_M_2D_P2_CH%d", ch);
            D_minus_M_2D_P2[ch] = new TH2D(hname, hname,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY);
        }
        */

    }

}



void
MinimizeFCNAxialVector::zero_V_PHYS_data() const
{
    // Set to zero
    for(int ch = 0; ch < number1DHists; ++ ch)
    {
        if(channel_enable_1D[ch] == 1)
        {
            /*
            for(Int_t bi = 1; bi <= V_PHYS_1D_P1[ch]->GetNbinsX(); ++ bi)
            {
                const Double_t zero = 0.0;
                V_PHYS_1D_P1[ch]->SetBinContent(bi, 1, zero);
            }

            for(Int_t bi = 1; bi <= V_PHYS_1D_P2[ch]->GetNbinsX(); ++ bi)
            {
                const Double_t zero = 0.0;
                V_PHYS_1D_P2[ch]->SetBinContent(bi, 1, zero);
            }
            */

            const int NUM_BINS_XY = 50;
            for(int i = 0; i < NUM_BINS_XY * NUM_BINS_XY; ++ i)
            {
                const double zero = 0.0;
                #if VECTOR_RANGE_CHECK
                V_PHYS_1D_P1_data[ch]->at(i) = zero;
                V_PHYS_1D_P2_data[ch]->at(i) = zero;
                #else
                V_PHYS_1D_P1_data[ch]->operator[](i) = zero;
                V_PHYS_1D_P2_data[ch]->operator[](i) = zero;
                #endif
            }
        }
    }

    /*
    // Set to zero
    for(int ch = 0; ch < number2DHists; ++ ch)
    {
        if(channel_enable_2D[ch] == 1)
        {
            for(Int_t bj = 1; bj <= V_PHYS_2D_P1[ch]->GetNbinsY(); ++ bj)
            {
                for(Int_t bi = 1; bi <= V_PHYS_2D_P1[ch]->GetNbinsX(); ++ bi)
                {
                    const Double_t zero = 0.0;
                    V_PHYS_2D_P1[ch]->SetBinContent(bi, bj, zero);
                }
            }
            for(Int_t bj = 1; bj <= V_PHYS_2D_P2[ch]->GetNbinsY(); ++ bj)
            {
                for(Int_t bi = 1; bi <= V_PHYS_2D_P2[ch]->GetNbinsX(); ++ bi)
                {
                    const Double_t zero = 0.0;
                    V_PHYS_2D_P2[ch]->SetBinContent(bi, bj, zero);
                }
            }
        }
    }
    */
}



void
MinimizeFCNAxialVector::set_D() const
{

    // loop over all channels
    for(int channel = 0; channel < number1DHists; ++ channel)
    {
        if(debuglevel >= 4)
        {
            std::cout << "channel=" << channel << std::endl;
        }

        // check channel enabled
        if(channel_enable_1D[channel] == 0)
        {
            if(debuglevel >= 5)
            {
                std::cout << "1D: channel " << channel << " disabled, skip" << std::endl;
            }
            continue;
        }

        // TODO: this no longer works because there are multiple phases
        // in the tmpData1D arrays
        // copy code from below to find correct object
        // TODO: phase 1 and phase 2 data objects
        //Double_t *tmpData1D_P1 = nullptr; //(TH1D*)allDataSamples1D->At(channel);
        //Double_t *tmpData1D_P2 = nullptr;
        //Double_t *tmpFakeData1D_P1 = nullptr; //(TH1D*)allFakeDataSamples1D->At(channel);
        //Double_t *tmpFakeData1D_P2 = nullptr;

        std::string histname = std::string(channel_histname_1D[channel]);
        std::string search_object_P1;
        std::string search_object_P2;
        if(g_mode_fake_data == false)
        {
            search_object_P1 = histname + std::string(DataFile) + "_P1";
            search_object_P2 = histname + std::string(DataFile) + "_P2";
        }
        else if(g_mode_fake_data == true)
        {
            search_object_P1 = histname + std::string("fakedata") + "_P1";
            search_object_P2 = histname + std::string("fakedata") + "_P2";
        }
        TH1D *tmpDataHist1D_P1 = nullptr;
        TH1D *tmpDataHist1D_P2 = nullptr;
        
        if(debuglevel >= 6)
        {
            std::cout << "search_object_P1=" << search_object_P1
                      << " search_object_P2=" << search_object_P2 << std::endl;
        }

        if(g_mode_fake_data == false)
        {
            tmpDataHist1D_P1 = (TH1D*)allDataSamples1D->FindObject(search_object_P1.c_str());
            tmpDataHist1D_P2 = (TH1D*)allDataSamples1D->FindObject(search_object_P2.c_str());
        }
        else if(g_mode_fake_data == true)
        {
            tmpDataHist1D_P1 = (TH1D*)allFakeDataSamples1D->FindObject(search_object_P1.c_str());
            tmpDataHist1D_P2 = (TH1D*)allFakeDataSamples1D->FindObject(search_object_P2.c_str());
        }

        if(tmpDataHist1D_P1 == nullptr)
        {
            std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
            throw "problem";
        }
        if(tmpDataHist1D_P2 == nullptr)
        {
            std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
            throw "problem";
        }

        //tmpData1D_P1 = new Double_t[tmpDataHist1D_P1->GetNbinsX()]; //tmpHistP1->Clone("tmpData1D_P1");
        //tmpData1D_P2 = new Double_t[tmpDataHist1D_P2->GetNbinsX()]; //tmpHistP2->Clone("tmpData1D_P2");
        for(Int_t bin_x{0}; bin_x < tmpDataHist1D_P1->GetNbinsX(); ++ bin_x)
        {
            //tmpData1D_P1[bin_x] = tmpDataHist1D_P1->GetBinContent(bin_x);
            //Int_t super_index = channel * 2 * 50 + bin_x;
            //Int_t super_index = bin_x;
            Double_t content = tmpDataHist1D_P1->GetBinContent(bin_x + 1);
            //D_1D_P1[channel]->SetBinContent(super_index + 1, 1, content_output);
            #if VECTOR_RANGE_CHECK
            D_1D_P1_data[channel]->at(bin_x) = content;
            #else
            D_1D_P1_data[channel]->operator[](bin_x) =  content;
            #endif
        }
        
        //std::cout << "LIST OF P2 DATA channel=" << channel << std::endl;
        //std::cin.get();
        for(Int_t bin_x{0}; bin_x < tmpDataHist1D_P2->GetNbinsX(); ++ bin_x)
        {
            //tmpData1D_P2[bin_x] = tmpDataHist1D_P2->GetBinContent(bin_x);
            //Int_t super_index = channel * 2 * 50 + 50 + bin_x;
            //Int_t super_index = bin_x;
            Double_t content = tmpDataHist1D_P2->GetBinContent(bin_x + 1);
            //D_1D_P2[channel]->SetBinContent(super_index + 1, 1, content_output);
            #if VECTOR_RANGE_CHECK
            D_1D_P2_data[channel]->at(bin_x) = content;
            #else
            D_1D_P2_data[channel]->operator[](bin_x) = content;
            #endif
            //std::cout << "bin_x=" << bin_x + 1 << " content=" << content_output << std::endl;
        }
    }
}


void 
MinimizeFCNAxialVector::set_M(const std::vector<double> &param) const
{

    // loop over all channels
    for(int channel = 0; channel < number1DHists; ++ channel)
    {
        if(debuglevel >= 4)
        {
            std::cout << "channel=" << channel << std::endl;
        }

        // check channel enabled
        if(channel_enable_1D[channel] == 0)
        {
            if(debuglevel >= 5)
            {
                std::cout << "1D: channel " << channel << " disabled, skip" << std::endl;
            }
            continue;
        }

        // reset M
        //for(Int_t bin_x{1}; bin_x <= M_1D_P1[channel]->GetNbinsX(); ++ bin_x)
        for(Int_t bin_x{0}; bin_x < M_1D_P1_data[channel]->size(); ++ bin_x)
        {
            //M_1D_P1[channel]->SetBinContent(bin_x, 1, 0.0);
            #if VECTOR_RANGE_CHECK
            M_1D_P1_data[channel]->at(bin_x) = 0.0;
            #else
            M_1D_P1_data[channel]->operator[](bin_x) = 0.0;
            #endif
        }
        
        // reset M
        //for(Int_t bin_x{1}; bin_x <= M_1D_P2[channel]->GetNbinsX(); ++ bin_x)
        for(Int_t bin_x{0}; bin_x < M_1D_P2_data[channel]->size(); ++ bin_x)
        {
            //M_1D_P2[channel]->SetBinContent(bin_x, 1, 0.0);
            #if VECTOR_RANGE_CHECK
            M_1D_P2_data[channel]->at(bin_x) = 0.0;
            #else
            M_1D_P2_data[channel]->operator[](bin_x) = 0.0;
            #endif
        }

    // can M change?
    // yes - if xi_31 parameter changes, which is detectable using g_pg
    // as shown above
    // can it change another way?
    // it depends on all "amplitude" parameters, so need to rebuild if any
    // of these parameters changes
    // just assume it always changes for now


        /*
        Double_t *tmpTotalMC1D_P1 = new Double_t[tmpDataHist1D_P1->GetNbinsX()];
        Double_t *tmpTotalMC1D_P2 = new Double_t[tmpDataHist1D_P2->GetNbinsX()];
        for(Int_t bin_x{1}; bin_x <= tmpDataHist1D_P1->GetNbinsX(); ++ bin_x)
        {
            tmpTotalMC1D_P1[bin_x] = 0.0;
        }
        for(Int_t bin_x{1}; bin_x <= tmpDataHist1D_P2->GetNbinsX(); ++ bin_x)
        {
            tmpTotalMC1D_P2[bin_x] = 0.0;
        }
        */



        // loop over all the parameters
        std::map<int, file_parameter>::iterator it{g_pg.file_params.begin()};
        for(; it != g_pg.file_params.end(); ++ it)
        {
            int paramNumberInt = -1;

            int paramNumber = it->second.paramNumber;
            //std::cout << "paramNumber=" << paramNumber << std::endl;
            bool paramEnabled = it->second.paramEnabled;
            bool paramEnabledP1 = it->second.paramEnabledP1;
            bool paramEnabledP2 = it->second.paramEnabledP2;
            double paramInitValue = it->second.paramInitValue;
            double paramInitError = it->second.paramInitError;
            int paramConstraintMode = it->second.paramConstraintMode;

            if(debuglevel >= 5)
            {
                std::cout << "paramNumber=" << paramNumber << std::endl;
            }

            // NOTE: notes on the V MATRIX method and the enable/disable
            // of parameters
            //
            // 1: how can parameters be enabled/disabled
            //
            // a:
            // There are two sources, the first is the
            // parameter_names.lst file, which has options to
            // enable/disable the parameter totally, or enable/disable
            // it for P1 and P2
            // 
            // b:
            // The other source is the channel enable/disable switch
            // in the header file
            // This affects the V_CHEN matrix because it switches
            // on/off groups of BINS in the chisquare calculation
            //
            // By contrast, (a) switches on and off parameters
            // which does not affect the number of bins only the bin
            // content

            if(paramEnabled == false)
            {
                if(debuglevel > 0)
                {
                    //std::cout << __func__ << " param " << paramNumber <<  " is disabled (overall) skip" << std::endl;
                    //std::cin.get();
                }
                continue;
            }
            if((paramEnabledP1 || paramEnabledP2) == false)
            {
                // if both are false
                if(debuglevel > 0)
                {
                    std::cout << __func__ << " param " << paramNumber << " is disabled for P1 and P2, skip" << std::endl;
                    std::cin.get();
                }
                continue;
            }
            // TODO: this correctly ignores any parameter which is disabled such that
            // paramEnabled == false
            // however, it also ignores parameters which are set as disabled for P1
            // and P2, and when these are irrelevent due to the value of gEnablePhaseX
            // so... the internal and external index will not match
            // need to add some code to fix this when the parameters are read from
            // file, (probably)
            // unless I just ignore that here... perhaps paramEnabled dictates
            // whether parameter is drawn and the phase1/phase2 enable flag
            // is to decide whether minuit does the fit or not (in which case
            // the param may still contribute to chisquare but may not be minimized
            // by minuit)


            // looping over parameters...
            // get a paramNumber
            // corresponds to a paramNumberInternal (minuit param number)
            // which is the index of the vector param where the current
            // activity can be found
            // also get a paramName, which is something like ac228_bi212_tl208_int
            // which is not very useful
            // but get a list of mc names, which for this example will be
            // ac228_int_rot, bi212_int_rot, tl208_int_rot
            // the corresponding MC samples names will be
            // theHistogram + BkgFiles[i] + "_Px_fit"
            // aka:
            // hTotalE_ + bi212_int_rot + _P1_fit
            // hTotalE_ + bi212_int_rot + _P2_fit
            // now, the paramNumber can be used to look up the parameter in
            // g_pg, from which we can obtain paramEnabled, paramEnabledP1
            // and paramEnabledP2, as well as the constraint mode, constraint
            // value and error
            // so, have to strip off the "_fit" ending, and strip off the
            // histogram (channel) name from the front, aka remove "hTotalE_"
            // then left with something like
            // bi212_int_rot_P1 and a channel number which comes from the loop
            // over allDataSamples1D->GetEntries()
            // EXCEPT this will now be wrong, because we have data for P1 and P2
            // so what we want to loop over is an index channel which runs from
            // 0 to number1DHists
            // so we have a channel index which is used to select the channel
            // which selects the index of allMCSamples1D[channel]
            // we then have to somehow "get" from this TObject
            // the histogram with the name bi212_int_rot_P1
            // and bi212_int_rot_P2
            // once we have done that we switch on paramEnabled, and paramEnabledP1
            // and paramEnabledP2 to decide whether or not we skip (in the case
            // where paramEnabled == false), apply a scaling with param[index]
            // for P1 and P2 if they are enabled, or simply add a value
            // where the scaling factor is given by paramInitValue for P1 and P2
            // if they are disabled
            // TODO: does this last part make sense? no it doesn't. skip if
            // P1/P2 not enabled

            std::vector<std::string>::iterator mc_name_it{it->second.MCNameList.begin()};
            for(; mc_name_it != it->second.MCNameList.end(); ++ mc_name_it)
            {
                std::string mc_name = *mc_name_it;
                std::string histname = std::string(channel_histname_1D[channel]);
                std::string search_object_P1 = histname + mc_name + "_P1_fit";
                std::string search_object_P2 = histname + mc_name + "_P2_fit";
                TH1D *tmpHist1D_P1 = nullptr;
                TH1D *tmpHist1D_P2 = nullptr;

                if(debuglevel >= 6)
                {
                    std::cout << "search_object_P1=" << search_object_P1
                              << " search_object_P2=" << search_object_P2 << std::endl;
                }

                paramNumberInt = g_pg.ExtToIntParamNumberMap.at(paramNumber);
                if(debuglevel >= 5)
                {
                    std::cout << "paramNumber=" << paramNumber << " -> " << paramNumberInt << std::endl;
                }

                // phase 1
                if(paramEnabledP1 == true)
                {
                    tmpHist1D_P1 = (TH1D*)allMCSamples1D[channel]->FindObject(search_object_P1.c_str());

                    if(tmpHist1D_P1 == nullptr)
                    {
                        std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
                        throw "problem";
                    }

                    Double_t scale_factor_P1 = param.at(paramNumberInt);
                    if(debuglevel >= 6)
                    {
                        std::cout << "enabled P1: scale factor P1: " << scale_factor_P1 << std::endl;
                    }

                    for(Int_t bin_x{0}; bin_x < tmpHist1D_P1->GetNbinsX(); ++ bin_x)
                    {
                        //Double_t content_input = M_1D_P1[channel]->GetBinContent(bin_x + 1, 1);
                        #if VECTOR_RANGE_CHECK
                        Double_t content_input = M_1D_P1_data[channel]->at(bin_x);
                        #else
                        Double_t content_input = M_1D_P1_data[channel]->operator[](bin_x);
                        #endif
                        Double_t content_add = scale_factor_P1 * tmpHist1D_P1->GetBinContent(bin_x + 1);
                        Double_t content_output = content_input + content_add;
                        //M_1D_P1[channel]->SetBinContent(bin_x + 1, 1, content_output);
                        #if VECTOR_RANGE_CHECK
                        M_1D_P1_data[channel]->at(bin_x) = content_output;
                        #else
                        M_1D_P1_data[channel]->operator[](bin_x) = content_output;
                        #endif
                        //std::cout << "content_input=" << content_input << " content_output=" << content_output << " content_add=" << content_add << std::endl;
                    }
                }
                else
                {
                    if(debuglevel >= 4)
                    {
                        std::cout << "disabled P1" << std::endl;
                    }
                }

                // phase 2
                if(paramEnabledP2 == true)
                {
                    tmpHist1D_P2 = (TH1D*)allMCSamples1D[channel]->FindObject(search_object_P2.c_str());
                
                    if(tmpHist1D_P2 == nullptr)
                    {
                        std::cout << "ERROR: Could not find object " << search_object_P2 << std::endl;
                        throw "problem";
                    }

                    Double_t scale_factor_P2 = param.at(paramNumberInt);
                    if(debuglevel >= 4)
                    {
                        std::cout << "enabled P2: scale factor P2: " << scale_factor_P2 << std::endl;
                    }

                    //std::cout << "super_index (start) : " << channel * 2 * 50 + 50 + 0 << std::endl;
                    for(Int_t bin_x{0}; bin_x < tmpHist1D_P2->GetNbinsX(); ++ bin_x)
                    {
                        //Int_t super_index = channel * 2 * 50 + 50 + bin_x;
                        //Double_t content_input = M_1D_P2[channel]->GetBinContent(bin_x + 1, 1);
                        #if VECTOR_RANGE_CHECK
                        Double_t content_input = M_1D_P2_data[channel]->at(bin_x);
                        #else
                        Double_t content_input = M_1D_P2_data[channel]->operator[](bin_x);
                        #endif
                        Double_t content_add = scale_factor_P2 * tmpHist1D_P2->GetBinContent(bin_x + 1);
                        Double_t content_output = content_input + content_add;
                        //M_1D_P2[channel]->SetBinContent(bin_x + 1, 1, content_output);
                        #if VECTOR_RANGE_CHECK
                        M_1D_P2_data[channel]->at(bin_x) = content_output;
                        #else
                        M_1D_P2_data[channel]->operator[](bin_x) = content_output;
                        #endif
                        //std::cout << "debug: " << "super_index=" << super_index << " content_input=" << content_input << " content_add=" << content_add << " content_output=" << content_output << " M:" << M->GetBinContent(super_index + 1, 1) << std::endl;

                    }
                }
                else
                {
                    if(debuglevel >= 4)
                    {
                        std::cout << "disabled P2" << std::endl;
                    }
                }




                // TODO: what is parameter is NOT enabled

            } // mc sample name iterator

        } // file_param iterator
    }
}


#if 1
void
MinimizeFCNAxialVector::set_D_minus_M() const
{
    // M is set, D is set

    
    for(int channel = 0; channel < number1DHists; ++ channel)
    {
        if(channel_enable_1D[channel] == 1)
        {
            // set D_minus_M
            //for(Int_t binx{1}; binx <= D_1D_P1[channel]->GetNbinsX(); ++ binx)
            for(Int_t binx{0}; binx < D_1D_P1_data[channel]->size(); ++ binx)
            {
                // P1
                {
                    //Double_t content_D = D_1D_P1[channel]->GetBinContent(binx, 1);
                    //Double_t content_M = M_1D_P1[channel]->GetBinContent(binx, 1);
                    #if VECTOR_RANGE_CHECK
                    Double_t content_D = D_1D_P1_data[channel]->at(binx);
                    Double_t content_M = M_1D_P1_data[channel]->at(binx);
                    #else
                    Double_t content_D = D_1D_P1_data[channel]->operator[](binx);
                    Double_t content_M = M_1D_P1_data[channel]->operator[](binx);
                    #endif
                    /*
                    if(content_D != 0.0)
                    {
                        std::cout << "binx=" << binx << " ~> " << content_D << " " << content_M << std::endl;
                    }
                    */
                    Double_t content_D_minus_M = content_D - content_M;
                    //D_minus_M_1D_P1[channel]->SetBinContent(binx, 1, content_D_minus_M);
                    #if VECTOR_RANGE_CHECK
                    D_minus_M_1D_P1_data[channel]->at(binx) = content_D_minus_M;
                    #else
                    D_minus_M_1D_P1_data[channel]->operator[](binx) = content_D_minus_M;
                    #endif
                }

                // P2
                {
                    //Double_t content_D = D_1D_P2[channel]->GetBinContent(binx, 1);
                    //Double_t content_M = M_1D_P2[channel]->GetBinContent(binx, 1);
                    #if VECTOR_RANGE_CHECK
                    Double_t content_D = D_1D_P2_data[channel]->at(binx);
                    Double_t content_M = M_1D_P2_data[channel]->at(binx);
                    #else
                    Double_t content_D = D_1D_P2_data[channel]->operator[](binx);
                    Double_t content_M = M_1D_P2_data[channel]->operator[](binx);
                    #endif
                    /*
                    if(content_D != 0.0)
                    {
                        std::cout << "binx=" << binx << " ~> " << content_D << " " << content_M << std::endl;
                    }
                    */
                    Double_t content_D_minus_M = content_D - content_M;
                    //D_minus_M_1D_P2[channel]->SetBinContent(binx, 1, content_D_minus_M);
                    #if VECTOR_RANGE_CHECK
                    D_minus_M_1D_P2_data[channel]->at(binx) = content_D_minus_M;
                    #else
                    D_minus_M_1D_P2_data[channel]->operator[](binx) = content_D_minus_M;
                    #endif
                }
            }
        }
    }
}
#endif


void
MinimizeFCNAxialVector::set_V_MATRIX() const
{

    // TODO: speed up this code by only building one matrix object
    // V_PHYS, and building it element by element, using symmetry
    // optimizations if possible. this avoids copying several V_MATRIX
    // objects into one final object
    // 
    // then further optimize by storing only D - M and computing the
    // matrix multiplication in place
    // this might not be possible because I need to invert the matrix
    // unless the inversion can be done in place? (might be slower
    // even if I can invent such an algorithm)


    // A funny thing COULD happen here
    // (but it is not very likely)
    // If the MATHMORE matrix becomes singular due to adding the components
    // of all the V_PHYS matricies
    // But the STAT matrix is non-singular
    // then the algorithm will blow up
    //
    // To fix this: Remove all V_PHYS matricies, and enter values
    // straight into a single V_PHYS matrix, or MATHMORE matrix
    // (need to do former to count number of required bins?)

    for(int channel = 0; channel < number1DHists; ++ channel)
    {

        if(channel_enable_1D[channel] == 1)
        {

            if(recalculate_V_PHYS_xD_Px_MATHMORE == true)
            {


                // set V_PHYS_SYS1
                // NOTE: moved to init, because values do not change

                #if 0
                // TODO: symmetry optimization
                for(Int_t biny{0}; biny < M_1D_P1_data[channel]->size(); ++ biny)
                {
                    //for(Int_t binx{0}; binx < M_1D_P1_data[channel]->size(); ++ binx)
                    for(Int_t binx{0}; binx < M_1D_P1_data[channel]->size(); ++ binx)
                    {
 
                        // P1
                        {
                            #if VECTOR_RANGE_CHECK
                            double coeff1 = systematic_offset_V_MATRIX_coeff_1D_P1[channel]->at(binx);
                            double coeff2 = systematic_offset_V_MATRIX_coeff_1D_P1[channel]->at(biny);
                            V_PHYS_SYS1_1D_P1_data[channel]->at(biny * 50 + binx) = coeff1 * coeff2;
                            #else
                            double coeff1 = systematic_offset_V_MATRIX_coeff_1D_P1[channel]->operator[](binx);
                            double coeff2 = systematic_offset_V_MATRIX_coeff_1D_P1[channel]->operator[](biny);
                            V_PHYS_SYS1_1D_P1_data[channel]->operator[](biny * 50 + binx) = coeff1 * coeff2;
                            #endif
                        }

                        // P2
                        {
                            #if VECTOR_RANGE_CHECK
                            double coeff1 = systematic_offset_V_MATRIX_coeff_1D_P2[channel]->at(binx);
                            double coeff2 = systematic_offset_V_MATRIX_coeff_1D_P2[channel]->at(biny);
                            V_PHYS_SYS1_1D_P2_data[channel]->at(biny * 50 + binx) = coeff1 * coeff2;
                            #else
                            double coeff1 = systematic_offset_V_MATRIX_coeff_1D_P2[channel]->operator[](binx);
                            double coeff2 = systematic_offset_V_MATRIX_coeff_1D_P2[channel]->operator[](biny);
                            V_PHYS_SYS1_1D_P2_data[channel]->operator[](biny * 50 + binx) = coeff1 * coeff2;
                            #endif
                        }


                    }
                }
                #endif


                // set V_PHYS_STAT
                int counter_P1 = 0;
                int counter_P2 = 0;
                //for(Int_t binx{1}; binx <= M_1D_P1[channel]->GetNbinsX(); ++ binx)
                for(Int_t binx{0}; binx < M_1D_P1_data[channel]->size(); ++ binx)
                {

                    // P1
                    {
                        #if VECTOR_RANGE_CHECK
                        Double_t content_M = M_1D_P1_data[channel]->at(binx);
                        Double_t content_D = D_1D_P1_data[channel]->at(binx);
                        #else
                        Double_t content_M = M_1D_P1_data[channel]->operator[](binx);
                        Double_t content_D = D_1D_P1_data[channel]->operator[](binx);
                        #endif
                        //Double_t content_M = M_1D_P1[channel]->GetBinContent(binx, 1);
                        //Double_t content_D = D_1D_P1[channel]->GetBinContent(binx, 1);
                        //if(content == 0.0) continue;
                        //if(content < 0.0) continue;
                        //if(content <= 0.0)
                        //{
                            //V_PHYS_STAT_1D_P1->SetBinContent(binx, binx, 1.0);
                        //}
                        //else
                        //{

                        Double_t stat = 1.0; // default value
                        if(content_M <= 0.0)
                        {
                            V_ENABLE_BIN_1D_P1[channel]->push_back(false);
                        }
                        else
                        {
                            V_ENABLE_BIN_1D_P1[channel]->push_back(true);

                            Double_t sigma_M = std::sqrt(content_M);
                            Double_t sigma_D = std::sqrt(content_D);
                            stat = sigma_M * sigma_M;

                            ++ counter_P1;

                        }
                            #if VECTOR_RANGE_CHECK
                            V_PHYS_STAT_1D_P1_data[channel]->at(binx * 50 + binx) = stat;
                            #else
                            V_PHYS_STAT_1D_P1_data[channel]->operator[](binx * 50 + binx) = stat;
                            #endif

                        //////Double_t sigma_M = std::sqrt(content_M);
                        //////Double_t sigma_D = std::sqrt(content_D);
                        //Double stat = 1.0 / (sigma * sigma);
                        //Double_t stat = 1.0 / std::abs(content);
                        // TODO; re-enable
                        ////Double_t stat = sigma_M * sigma_M; // + sigma_D * sigma_D;
                        //V_PHYS_STAT_1D_P1->SetBinContent(binx, binx, stat);
                        ////if(stat != 0.0) stat = 1.0 / stat;
                        ////else stat = 1.0;
                        //////Double_t stat = sigma_M * sigma_M;

                        /*
                        if(std::isnan(stat))
                        {
                            std::cout << "binx=" << binx << " sigma_M=" << sigma_M << " content_M=" << content_M << std::endl;
                        }
                        */

                        //if(stat == 0.0)
                        //{
                        //    stat = 1.0;
                        //    V_ENABLE_BIN_1D_P1[channel]->push_back(false);
                        //}
                        //else
                        //{
                        //    V_ENABLE_BIN_1D_P1[channel]->push_back(true);
                        //    ++ counter_P1;
                        //}
                        //V_PHYS_STAT_1D_P1[channel]->SetBinContent(binx, binx, stat);
                        //#if VECTOR_RANGE_CHECK
                        //V_PHYS_STAT_1D_P1_data[channel]->at(binx * 50 + binx) = stat;
                        //#else
                        //V_PHYS_STAT_1D_P1_data[channel]->operator[](binx * 50 + binx) = stat;
                        //#endif
                        //}
                    }

                    // P2
                    {
                        //Double_t content_M = M_1D_P2[channel]->GetBinContent(binx, 1);
                        //Double_t content_D = D_1D_P2[channel]->GetBinContent(binx, 1);
                        #if VECTOR_RANGE_CHECK
                        Double_t content_M = M_1D_P2_data[channel]->at(binx);
                        Double_t content_D = D_1D_P2_data[channel]->at(binx);
                        #else
                        Double_t content_M = M_1D_P2_data[channel]->operator[](binx);
                        Double_t content_D = D_1D_P2_data[channel]->operator[](binx);
                        #endif

                        #if 0
                        //if(content == 0.0) continue;
                        //if(content < 0.0) continue;
                        //if(content <= 0.0)
                        //{
                            //V_PHYS_STAT_1D_P2->SetBinContent(binx, binx, 1.0);
                        //}
                        //else
                        //{
                        Double_t sigma_M = std::sqrt(content_M);
                        Double_t sigma_D = std::sqrt(content_D);
                        //Double stat = 1.0 / (sigma * sigma);
                        //Double_t stat = 1.0 / std::abs(content);
                        // TODO: re-enable
                        ////Double_t stat = sigma_M * sigma_M; // + sigma_D * sigma_D;
                        //V_PHYS_STAT_1D_P2->SetBinContent(binx, binx, stat);
                        ////if(stat != 0.0) stat = 1.0 / stat;
                        ////else stat = 1.0;
                        Double_t stat = sigma_M * sigma_M;
                        // TODO: this might not be the best way to detect signularity
                        // may be more sensible to use if D==0 or M==0, and then or sigmaM==0
                        if(stat == 0.0)
                        {
                            stat = 1.0;
                            V_ENABLE_BIN_1D_P2[channel]->push_back(false);
                        }
                        else
                        {
                            V_ENABLE_BIN_1D_P2[channel]->push_back(true);
                            ++ counter_P2;
                        }
                        //V_PHYS_STAT_1D_P2[channel]->SetBinContent(binx, binx, stat);
                        #if VECTOR_RANGE_CHECK
                        V_PHYS_STAT_1D_P2_data[channel]->at(binx * 50 + binx) =  stat;
                        #else
                        V_PHYS_STAT_1D_P2_data[channel]->operator[](binx * 50 + binx) =  stat;
                        #endif
                        //}
                        #endif

                        Double_t stat = 1.0;
                        // TODO: this might not be the best way to detect signularity
                        // may be more sensible to use if D==0 or M==0, and then or sigmaM==0
                        if(content_M <= 0.0)
                        {
                            V_ENABLE_BIN_1D_P2[channel]->push_back(false);
                        }
                        else
                        {
                            V_ENABLE_BIN_1D_P2[channel]->push_back(true);
                        
                            Double_t sigma_M = std::sqrt(content_M);
                            Double_t sigma_D = std::sqrt(content_D);
                            stat = sigma_M * sigma_M;

                            ++ counter_P2;

                        }
                            #if VECTOR_RANGE_CHECK
                            V_PHYS_STAT_1D_P2_data[channel]->at(binx * 50 + binx) =  stat;
                            #else
                            V_PHYS_STAT_1D_P2_data[channel]->operator[](binx * 50 + binx) =  stat;
                            #endif
                    }
                }
                //std::cout << "P1: number of enabled bins = " << counter_P1 << std::endl;
                //std::cout << "P2: number of enabled bins = " << counter_P2 << std::endl;
                ///////////////////////////////////////////////////////////////////
                // CERN ROOT MathMore Matrix Lib Objects
                ///////////////////////////////////////////////////////////////////


                //if(V_PHYS_1D_P1_MATHMORE[0] == nullptr)
                //if(recalculate_V_PHYS_xD_Px_MATHMORE == true)
                //{
                //std::cout << "Alloc V_PHYS_MATHMORE channel=" << channel << std::endl;

                //const Int_t NUM_BINS_XY = 50;

                //for(int ch = 0; ch < number1DHists; ++ ch)
                //{


                #if MEASURE_FUNCTION_CALL_TIME
                std::chrono::system_clock::time_point start_time = std::chrono::high_resolution_clock::now();
                #endif
                {
                    if(V_PHYS_1D_P1_MATHMORE[channel] != nullptr)
                    {
                        if(counter_P1 != V_PHYS_1D_P1_MATHMORE[channel]->GetNrows())
                        {
                            // nrows always == ncols
                            delete V_PHYS_1D_P1_MATHMORE[channel];
                            V_PHYS_1D_P1_MATHMORE[channel] = new TMatrixD(counter_P1, counter_P1); // NUM_BINS_XY, NUM_BINS_XY);
                        }
                    }
                    else
                    {
                    std::cout << "MATHMORE: P1: " << V_PHYS_1D_P1_MATHMORE[0] << std::endl;
                        V_PHYS_1D_P1_MATHMORE[channel] = new TMatrixD(counter_P1, counter_P1); // NUM_BINS_XY, NUM_BINS_XY);
                        std::cout << "ALLOC " << __func__ << std::endl;
                    }

                    if(V_PHYS_1D_P2_MATHMORE[channel] != nullptr)
                    {
                        if(counter_P2 != V_PHYS_1D_P2_MATHMORE[channel]->GetNrows())
                        {
                            delete V_PHYS_1D_P2_MATHMORE[channel];
                            V_PHYS_1D_P2_MATHMORE[channel] = new TMatrixD(counter_P2, counter_P2); //NUM_BINS_XY, NUM_BINS_XY);
                        }
                    }
                    else
                    {
                    std::cout << "MATHMORE: P2: " << V_PHYS_1D_P2_MATHMORE[0] << std::endl;
                        V_PHYS_1D_P2_MATHMORE[channel] = new TMatrixD(counter_P2, counter_P2); //NUM_BINS_XY, NUM_BINS_XY);
                        std::cout << "ALLOC " << __func__ << std::endl;
                    }
                    // TODO: must be able to optimize this
                }
                #if MEASURE_FUNCTION_CALL_TIME
                std::chrono::system_clock::time_point end_time = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> runtime_sec = end_time - start_time;
                std::cout << "V_PHYS_1D_Px_MATHMORE allocation time: " << 1.0e+06 * runtime_sec.count() << " microsecond" << std::endl;
                #endif

                //std::cout << "size: " << V_PHYS_1D_P1_MATHMORE[channel]->GetNrows() << " " << V_PHYS_1D_P1_MATHMORE[channel]->GetNcols() << std::endl;
                //std::cout << "size: " << V_PHYS_1D_P2_MATHMORE[channel]->GetNrows() << " " << V_PHYS_1D_P2_MATHMORE[channel]->GetNcols() << std::endl;
                //std::cin.get();
                //}

                /*
                for(int ch = 0; ch < number1DHists; ++ ch)
                {
                    for(Int_t iy = 0; iy < V_PHYS_1D_P1_MATHMORE[ch]->GetNrows(); ++ iy)
                    {
                        for(Int_t ix = 0; ix < V_PHYS_1D_P1_MATHMORE[ch]->GetNcols(); ++ ix)
                        {
                            V_PHYS_1D_P1_MATHMORE[ch]->operator[](iy).operator[](ix) = 0.0;
                            V_PHYS_1D_P2_MATHMORE[ch]->operator[](iy).operator[](ix) = 0.0;
                        }
                    }
                }
                */
                //}


                // TODO: fix problems here
                // if channel is disabled, ignore it
                // should be done in a big output loop?
                // DONE
                // there is a reference to channel here
                // and yet there is a loop over ch above
                // how can this make sense? is it a bug?
                // FIXED
                // TODO: ALL CHANNELS ARE CURRENTLY ENABLED?
   
                #define DRAWVPHYSMATRIX 1

                #if DRAWVPHYSMATRIX
                TH2D *draw_V_PHYS_STAT_P1 = nullptr;
                TH2D *draw_V_PHYS_SYS1_P1 = nullptr;
                TH2D *draw_V_PHYS_SYS2_P1 = nullptr;
                TH2D *draw_V_PHYS_SYS3_P1 = nullptr;
                TH2D *draw_V_PHYS_SYS4_P1 = nullptr;
                TH2D *draw_V_PHYS_P1 = nullptr;
                TH2D *draw_V_PHYSINV_P1 = nullptr;
                TH2D *draw_V_PHYS_STAT_P2 = nullptr;
                TH2D *draw_V_PHYS_SYS1_P2 = nullptr;
                TH2D *draw_V_PHYS_SYS2_P2 = nullptr;
                TH2D *draw_V_PHYS_SYS3_P2 = nullptr;
                TH2D *draw_V_PHYS_SYS4_P2 = nullptr;
                TH2D *draw_V_PHYS_P2 = nullptr;
                TH2D *draw_V_PHYSINV_P2 = nullptr;
                if(DRAW_V_PHYS_MATRIX == true)
                {
                    if(channel == 1)
                    {
                        draw_V_PHYS_STAT_P1 = new TH2D("draw_V_PHYS_STAT_P1", "draw_V_PHYS_STAT_P1",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_SYS1_P1 = new TH2D("draw_V_PHYS_SYS1_P1", "draw_V_PHYS_SYS1_P1",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_SYS2_P1 = new TH2D("draw_V_PHYS_SYS2_P1", "draw_V_PHYS_SYS2_P1",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_SYS3_P1 = new TH2D("draw_V_PHYS_SYS3_P1", "draw_V_PHYS_SYS3_P1",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_SYS4_P1 = new TH2D("draw_V_PHYS_SYS4_P1", "draw_V_PHYS_SYS4_P1",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_P1 =      new TH2D("draw_V_PHYS_P1", "draw_V_PHYS_P1",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYSINV_P1 =   new TH2D("draw_V_PHYSINV_P1", "draw_V_PHYSINV_P1",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_STAT_P2 = new TH2D("draw_V_PHYS_STAT_P2", "draw_V_PHYS_STAT_P2",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_SYS1_P2 = new TH2D("draw_V_PHYS_SYS1_P2", "draw_V_PHYS_SYS1_P2",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_SYS2_P2 = new TH2D("draw_V_PHYS_SYS2_P2", "draw_V_PHYS_SYS2_P2",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_SYS3_P2 = new TH2D("draw_V_PHYS_SYS3_P2", "draw_V_PHYS_SYS3_P2",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_SYS4_P2 = new TH2D("draw_V_PHYS_SYS4_P2", "draw_V_PHYS_SYS4_P2",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_P2 =      new TH2D("draw_V_PHYS_P2", "draw_V_PHYS_P2",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYSINV_P2 =   new TH2D("draw_V_PHYSINV_P2", "draw_V_PHYSINV_P2",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        std::cout << "ALLOC " << __func__ << std::endl;
                    }
                }
                #endif


                ///////////////////////////////////////////////////////////
                // Phase 1

                int j_counter = 0;
                //for(Int_t j = 0; j < V_PHYS_STAT_1D_P1[channel]->GetNbinsY(); ++ j)
                for(Int_t j = 0; j < 50; ++ j)
                {
                    if(V_ENABLE_BIN_1D_P1[channel]->at(j) == true)
                    {
                        // do nothing
                    }
                    else
                    {
                        continue;
                    }

                    int i_counter = 0;
                    for(Int_t i = 0; i < 50; ++ i)
                    {
                        if(V_ENABLE_BIN_1D_P1[channel]->at(i) == true)
                        {
                            // do nothing
                        }
                        else
                        {
                            continue;
                        }

                        //Double_t content = V_PHYS_STAT_1D_P1[channel]->GetBinContent(i + 1, j + 1);
                        #if VECTOR_RANGE_CHECK
                        //Double_t content = V_PHYS_STAT_1D_P1_data[channel]->at(i + j * 50);
                        double c1 = 0.0;
                        double c2 = 0.0;
                        double c3 = 0.0;
                        double c4 = 0.0;
                        double c5 = 0.0;
                        if(V_ENABLE_STAT == true)
                        {
                            c1 = V_PHYS_STAT_1D_P1_data[channel]->at(i + j * 50);
                        }
                        if(V_ENABLE_SYSALL == true)
                        {
                            if(V_ENABLE_SYS1 == true)
                            {
                                c2 = V_PHYS_SYS1_1D_P1_data[channel]->at(i + j * 50);
                            }
                            if(V_ENABLE_SYS2 == true)
                            {
                                c3 = V_PHYS_SYS2_1D_P1_data[channel]->at(i + j * 50);
                            }
                            if(V_ENABLE_SYS3 == true)
                            {
                                c4 = V_PHYS_SYS3_1D_P1_data[channel]->at(i + j * 50);
                            }
                            if(V_ENABLE_SYS4 == true)
                            {
                                c5 = V_PHYS_SYS4_1D_P1_data[channel]->at(i + j * 50);
                            }
                        }
                        Double_t content = c1 + c2 + c3 + c4 + c5;
                        #else
                        //Double_t content = V_PHYS_STAT_1D_P1_data[channel]->operator[](i + j * 50);
                        double c1 = 0.0;
                        double c2 = 0.0;
                        double c3 = 0.0;
                        double c4 = 0.0;
                        double c5 = 0.0;
                        if(V_ENABLE_STAT == true)
                        {
                            c1 = V_PHYS_STAT_1D_P1_data[channel]->operator[](i + j * 50);
                        }
                        if(V_ENABLE_SYSALL == true)
                        {
                            if(V_ENABLE_SYS1 == true)
                            {
                                c2 = V_PHYS_SYS1_1D_P1_data[channel]->operator[](i + j * 50);
                            }
                            if(V_ENABLE_SYS2 == true)
                            {
                                c3 = V_PHYS_SYS2_1D_P1_data[channel]->operator[](i + j * 50);
                            }
                            if(V_ENABLE_SYS3 == true)
                            {
                                c4 = V_PHYS_SYS3_1D_P1_data[channel]->operator[](i + j * 50);
                            }
                            if(V_ENABLE_SYS4 == true)
                            {
                                c5 = V_PHYS_SYS4_1D_P1_data[channel]->operator[](i + j * 50);
                            }
                        }
                        Double_t content = c1 + c2 + c3 + c4 + c5;
                        #endif
                        //std::cout << "i=" << i << " j=" << j << " i_counter=" << i_counter << " j_counter=" << j_counter << " content=" << content << std::endl;
                        V_PHYS_1D_P1_MATHMORE[channel]->operator[](j_counter).operator[](i_counter) = content;
                        //std::cout << "j=" << j << " i=" << i << " " << content << std::endl;

                        std::cout << "P1" << std::endl;
                        std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c1=" << c1 << std::endl;
                        std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c2=" << c2 << std::endl;
                        std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c3=" << c3 << std::endl;
                        std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c4=" << c4 << std::endl;
                        std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c5=" << c5 << std::endl;
                        std::cin.get();

                        #if DRAWVPHYSMATRIX
                        if(DRAW_V_PHYS_MATRIX == true)
                        {
                            if(channel == 1)
                            {
                                draw_V_PHYS_STAT_P1->SetBinContent(i + 1, j + 1, c1);
                                draw_V_PHYS_SYS1_P1->SetBinContent(i + 1, j + 1, c2);
                                draw_V_PHYS_SYS2_P1->SetBinContent(i + 1, j + 1, c3);
                                draw_V_PHYS_SYS3_P1->SetBinContent(i + 1, j + 1, c4);
                                draw_V_PHYS_SYS4_P1->SetBinContent(i + 1, j + 1, c5);
                                draw_V_PHYS_P1->SetBinContent(i + 1, j + 1, content);

                            /*
                                std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c1=" << c1 << std::endl;
                                std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c2=" << c2 << std::endl;
                                std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c3=" << c3 << std::endl;
                                std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c5=" << c5 << std::endl;
                                //std::cout << "channel=" << channel << " i=" << i << " j=" << j << " content=" << content << std::endl;

                                std::cout << "V_ENABLE_SYSALL=" << V_ENABLE_SYSALL << std::endl;
                                std::cout << "V_ENABLE_SYS1=" << V_ENABLE_SYS1 << std::endl;
                                std::cout << "V_ENABLE_SYS2=" << V_ENABLE_SYS2 << std::endl;
                                std::cout << "V_ENABLE_SYS3=" << V_ENABLE_SYS3 << std::endl;
                                std::cout << "V_ENABLE_SYS4=" << V_ENABLE_SYS4 << std::endl;
                                //std::cin.get();
                            */
                            }
                        }
                        #endif


                        ++ i_counter;
                    }

                    ++ j_counter;
                }

                ///////////////////////////////////////////////////////////
                // Phase 2

                j_counter = 0;
                //for(Int_t j = 0; j < V_PHYS_STAT_1D_P2[channel]->GetNbinsY(); ++ j)
                for(Int_t j = 0; j < 50; ++ j)
                {
                    if(V_ENABLE_BIN_1D_P2[channel]->at(j) == true)
                    {
                        // do nothing
                    }
                    else
                    {
                        continue;
                    }

                    int i_counter = 0;
                    for(Int_t i = 0; i < 50; ++ i)
                    {
                        if(V_ENABLE_BIN_1D_P2[channel]->at(i) == true)
                        {
                            // do nothing
                        }
                        else
                        {
                            continue;
                        }

                        //Double_t content = V_PHYS_STAT_1D_P2[channel]->GetBinContent(i + 1, j + 1);
                        #if VECTOR_RANGE_CHECK
                        //Double_t content = V_PHYS_STAT_1D_P2_data[channel]->at(i + j * 50);
                        double c1 = 0.0;
                        double c2 = 0.0;
                        double c3 = 0.0;
                        double c4 = 0.0;
                        double c5 = 0.0;
                        if(V_ENABLE_STAT == true)
                        {
                            c1 = V_PHYS_STAT_1D_P2_data[channel]->at(i + j * 50);
                        }
                        if(V_ENABLE_SYSALL == true)
                        {
                            if(V_ENABLE_SYS1 == true)
                            {
                                c2 = V_PHYS_SYS1_1D_P2_data[channel]->at(i + j * 50);
                            }
                            if(V_ENABLE_SYS2 == true)
                            {
                                c3 = V_PHYS_SYS2_1D_P2_data[channel]->at(i + j * 50);
                            }
                            if(V_ENABLE_SYS3 == true)
                            {
                                c4 = V_PHYS_SYS3_1D_P2_data[channel]->at(i + j * 50);
                            }
                            if(V_ENABLE_SYS4 == true)
                            {
                                c5 = V_PHYS_SYS4_1D_P2_data[channel]->at(i + j * 50);
                            }
                        }
                        double content = c1 + c2 + c3 + c4 + c5;
                        #else
                        //Double_t content = V_PHYS_STAT_1D_P2_data[channel]->operator[](i + j * 50);
                        double c1 = 0.0;
                        double c2 = 0.0;
                        double c3 = 0.0;
                        double c4 = 0.0;
                        double c5 = 0.0;
                        if(V_ENABLE_STAT == true)
                        {
                            c1 = V_PHYS_STAT_1D_P2_data[channel]->operator[](i + j * 50);
                        }
                        if(V_ENABLE_SYSALL == true)
                        {
                            if(V_ENABLE_SYS1 == true)
                            {
                                c2 = V_PHYS_SYS1_1D_P2_data[channel]->operator[](i + j * 50);
                            }
                            if(V_ENABLE_SYS2 == true)
                            {
                                c3 = V_PHYS_SYS2_1D_P2_data[channel]->operator[](i + j * 50);
                            }
                            if(V_ENABLE_SYS3 == true)
                            {
                                c4 = V_PHYS_SYS3_1D_P2_data[channel]->operator[](i + j * 50);
                            }
                            if(V_ENABLE_SYS4 == true)
                            {
                                c5 = V_PHYS_SYS4_1D_P2_data[channel]->operator[](i + j * 50);
                            }
                        }
                        double content = c1 + c2 + c3 + c4 + c5;
                        #endif
                        //std::cout << "i=" << i << " j=" << j << " i_counter=" << i_counter << " j_counter=" << j_counter << " content=" << content << std::endl;
                        V_PHYS_1D_P2_MATHMORE[channel]->operator[](j_counter).operator[](i_counter) = content;

                        std::cout << "P2" << std::endl;
                        std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c1=" << c1 << std::endl;
                        std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c2=" << c2 << std::endl;
                        std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c3=" << c3 << std::endl;
                        std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c4=" << c4 << std::endl;
                        std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c5=" << c5 << std::endl;
                        std::cin.get();

                        #if DRAWVPHYSMATRIX
                        if(DRAW_V_PHYS_MATRIX == true)
                        {
                            if(channel == 1)
                            {
                                draw_V_PHYS_STAT_P2->SetBinContent(i + 1, j + 1, c1);
                                draw_V_PHYS_SYS1_P2->SetBinContent(i + 1, j + 1, c2);
                                draw_V_PHYS_SYS2_P2->SetBinContent(i + 1, j + 1, c3);
                                draw_V_PHYS_SYS3_P2->SetBinContent(i + 1, j + 1, c4);
                                draw_V_PHYS_SYS4_P2->SetBinContent(i + 1, j + 1, c5);
                                draw_V_PHYS_P2->SetBinContent(i + 1, j + 1, content);

                            /*
                                std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c1=" << c1 << std::endl;
                                std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c2=" << c2 << std::endl;
                                std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c3=" << c3 << std::endl;
                                std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c4=" << c4 << std::endl;
                                std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c5=" << c5 << std::endl;
                                //std::cout << "channel=" << channel << " i=" << i << " j=" << j << " content=" << content << std::endl;
                            */
                            }
                        }
                        #endif


                        ++ i_counter;
                    }

                    ++ j_counter;
                }

                #if DRAWVPHYSMATRIX
                if(DRAW_V_PHYS_MATRIX == true)
                {
                    if(channel == 1)
                    {
                        TCanvas *canvas_V_PHYS_STAT_P1 = new TCanvas("canvas_V_PHYS_STAT_P1", "canvas_V_PHYS_STAT_P1");
                        draw_V_PHYS_STAT_P1->Draw("colz");
                        canvas_V_PHYS_STAT_P1->SaveAs("draw_V_PHYS_STAT_P1.png");

                        TCanvas *canvas_V_PHYS_SYS1_P1 = new TCanvas("canvas_V_PHYS_SYS1_P1", "canvas_V_PHYS_SYS1_P1");
                        draw_V_PHYS_SYS1_P1->Draw("colz");
                        canvas_V_PHYS_SYS1_P1->SaveAs("draw_V_PHYS_SYS1_P1.png");

                        TCanvas *canvas_V_PHYS_SYS2_P1 = new TCanvas("canvas_V_PHYS_SYS2_P1", "canvas_V_PHYS_SYS2_P1");
                        draw_V_PHYS_SYS2_P1->Draw("colz");
                        canvas_V_PHYS_SYS2_P1->SaveAs("draw_V_PHYS_SYS2_P1.png");

                        TCanvas *canvas_V_PHYS_SYS3_P1 = new TCanvas("canvas_V_PHYS_SYS3_P1", "canvas_V_PHYS_SYS3_P1");
                        draw_V_PHYS_SYS3_P1->Draw("colz");
                        canvas_V_PHYS_SYS3_P1->SaveAs("draw_V_PHYS_SYS3_P1.png");

                        TCanvas *canvas_V_PHYS_SYS4_P1 = new TCanvas("canvas_V_PHYS_SYS4_P1", "canvas_V_PHYS_SYS4_P1");
                        draw_V_PHYS_SYS4_P1->Draw("colz");
                        canvas_V_PHYS_SYS4_P1->SaveAs("draw_V_PHYS_SYS4_P1.png");

                        TCanvas *canvas_V_PHYS_P1 = new TCanvas("canvas_V_PHYS_P1", "canvas_V_PHYS_P1");
                        draw_V_PHYS_P1->Draw("colz");
                        canvas_V_PHYS_P1->SaveAs("draw_V_PHYS_P1.png");

                        TCanvas *canvas_V_PHYS_STAT_P2 = new TCanvas("canvas_V_PHYS_STAT_P2", "canvas_V_PHYS_STAT_P2");
                        draw_V_PHYS_STAT_P2->Draw("colz");
                        canvas_V_PHYS_STAT_P2->SaveAs("draw_V_PHYS_STAT_P2.png");

                        TCanvas *canvas_V_PHYS_SYS1_P2 = new TCanvas("canvas_V_PHYS_SYS1_P2", "canvas_V_PHYS_SYS1_P2");
                        draw_V_PHYS_SYS1_P2->Draw("colz");
                        canvas_V_PHYS_SYS1_P2->SaveAs("draw_V_PHYS_SYS1_P2.png");

                        TCanvas *canvas_V_PHYS_SYS2_P2 = new TCanvas("canvas_V_PHYS_SYS2_P2", "canvas_V_PHYS_SYS2_P2");
                        draw_V_PHYS_SYS2_P2->Draw("colz");
                        canvas_V_PHYS_SYS2_P2->SaveAs("draw_V_PHYS_SYS2_P2.png");

                        TCanvas *canvas_V_PHYS_SYS3_P2 = new TCanvas("canvas_V_PHYS_SYS3_P2", "canvas_V_PHYS_SYS3_P2");
                        draw_V_PHYS_SYS3_P2->Draw("colz");
                        canvas_V_PHYS_SYS3_P2->SaveAs("draw_V_PHYS_SYS3_P2.png");

                        TCanvas *canvas_V_PHYS_SYS4_P2 = new TCanvas("canvas_V_PHYS_SYS4_P2", "canvas_V_PHYS_SYS4_P2");
                        draw_V_PHYS_SYS4_P2->Draw("colz");
                        canvas_V_PHYS_SYS4_P2->SaveAs("draw_V_PHYS_SYS4_P2.png");

                        TCanvas *canvas_V_PHYS_P2 = new TCanvas("canvas_V_PHYS_P2", "canvas_V_PHYS_P2");
                        draw_V_PHYS_P2->Draw("colz");
                        canvas_V_PHYS_P2->SaveAs("draw_V_PHYS_P2.png");
                        std::cout << "ALLOC " << __func__ << std::endl;

                        //std::cout << "printed all canvas" << std::endl;
                        //std::cin.get();
                    }
                }
                #endif

                /*
                for(Int_t iy = 0; iy < V_PHYS_1D_P1_MATHMORE[channel]->GetNrows(); ++ iy)
                {
                    for(Int_t ix = 0; ix < V_PHYS_1D_P1_MATHMORE[channel]->GetNcols(); ++ ix)
                    {
                        Double_t vP1 = V_PHYS_1D_P1_MATHMORE[channel]->operator[](iy).operator[](ix);
                        if(vP1 != 0.0)
                        {
                            std::cout << "P1 ix=" << ix << " iy=" << iy << " " << vP1 << std::endl;
                        }
                    }
                }
                for(Int_t iy = 0; iy < V_PHYS_1D_P2_MATHMORE[channel]->GetNrows(); ++ iy)
                {
                    for(Int_t ix = 0; ix < V_PHYS_1D_P2_MATHMORE[channel]->GetNcols(); ++ ix)
                    {
                        Double_t vP2 = V_PHYS_1D_P2_MATHMORE[channel]->operator[](iy).operator[](ix);
                        if(vP2 != 0.0)
                        {
                            std::cout << "P2 ix=" << ix << " iy=" << iy << " " << vP2 << std::endl;
                        }
                    }
                }
                */
                //std::cout << "Matrix size: " << V_PHYS_1D_P1_MATHMORE[channel]->GetNrows() << " x " << V_PHYS_1D_P1_MATHMORE[channel]->GetNcols() << std::endl;
                //std::cout << "Matrix size: " << V_PHYS_1D_P2_MATHMORE[channel]->GetNrows() << " x " << V_PHYS_1D_P2_MATHMORE[channel]->GetNcols() << std::endl;
                #define MEASURE_FUNCTION_CALL_TIME_MATRIX_INVERT 0
                #if MEASURE_FUNCTION_CALL_TIME_MATRIX_INVERT
                std::cout << "Start Invert" << std::endl;
                std::chrono::system_clock::time_point start_time = std::chrono::high_resolution_clock::now();
                #endif
                //std::cout << "P1 V_MATRIX size: " << V_PHYS_1D_P1_MATHMORE[channel]->GetNrows() << " " << V_PHYS_1D_P1_MATHMORE[channel]->GetNcols() << std::endl;
                if(V_PHYS_1D_P1_MATHMORE[channel]->Determinant() == 0.0)
                {
                    std::cout << "P1 V_MATRIX size: " << V_PHYS_1D_P1_MATHMORE[channel]->GetNrows() << " " << V_PHYS_1D_P1_MATHMORE[channel]->GetNcols() << std::endl;
                    for(int j = 0; j < V_PHYS_1D_P1_MATHMORE[channel]->GetNcols(); ++ j)
                    {
                        for(int i = 0; i < V_PHYS_1D_P1_MATHMORE[channel]->GetNrows(); ++ i)
                        {
                            std::cout << "i=" << i << " j=" << j << " " << V_PHYS_1D_P1_MATHMORE[channel]->operator[](j).operator[](i) << std::endl;
                        }
                    }
                }

                //TMatrixD copy(*V_PHYS_1D_P1_MATHMORE[channel]);
                V_PHYS_1D_P1_MATHMORE[channel]->Invert();
                //TMatrixD multiplymat = copy * *V_PHYS_1D_P1_MATHMORE[channel];
                //std::cout << "Elements after mult" << std::endl;
                //    for(int j = 0; j < V_PHYS_1D_P1_MATHMORE[channel]->GetNcols(); ++ j)
                //    {
                //        for(int i = 0; i < V_PHYS_1D_P1_MATHMORE[channel]->GetNrows(); ++ i)
                //        {
                //            std::cout << "i=" << i << " j=" << j << " " << multiplymat.operator[](j).operator[](i) << std::endl;
                //        }
                //    }
                //    std::cin.get();
                //std::cout << "Next Invert" << std::endl;
                //std::cout << "P2 V_MATRIX size: " << V_PHYS_1D_P2_MATHMORE[channel]->GetNrows() << " " << V_PHYS_1D_P2_MATHMORE[channel]->GetNcols() << std::endl;
                V_PHYS_1D_P2_MATHMORE[channel]->Invert();
                #if MEASURE_FUNCTION_CALL_TIME_MATRIX_INVERT
                std::chrono::system_clock::time_point end_time = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> runtime_sec = end_time - start_time;
                std::cout << "Done Invert, time=" << 1.0e+06 * runtime_sec.count() << " microsecond" << std::endl;
                #endif




                j_counter = 0;
                for(Int_t j = 0; j < 50; ++ j)
                {
                    if(V_ENABLE_BIN_1D_P1[channel]->at(j) == true)
                    {
                        // do nothing
                    }
                    else
                    {
                        continue;
                    }

                    int i_counter = 0;
                    for(Int_t i = 0; i < 50; ++ i)
                    {
                        if(V_ENABLE_BIN_1D_P1[channel]->at(i) == true)
                        {
                            // do nothing
                        }
                        else
                        {
                            continue;
                        }

                        double content = V_PHYS_1D_P1_MATHMORE[channel]->operator[](j_counter).operator[](i_counter);

                        #if DRAWVPHYSMATRIX
                        if(DRAW_V_PHYS_MATRIX == true)
                        {
                            if(channel == 1)
                            {
                                draw_V_PHYSINV_P1->SetBinContent(i + 1, j + 1, content);
                                //std::cout << content << "\t";
                            }
                        }
                        #endif
                        ++ i_counter;
                    }
                    #if DRAWVPHYSMATRIX
                    if(DRAW_V_PHYS_MATRIX == true)
                    {
                        //std::cout << std::endl;
                    }
                    #endif

                    ++ j_counter;
                }


                j_counter = 0;
                for(Int_t j = 0; j < 50; ++ j)
                {
                    if(V_ENABLE_BIN_1D_P2[channel]->at(j) == true)
                    {
                        // do nothing
                    }
                    else
                    {
                        continue;
                    }

                    int i_counter = 0;
                    for(Int_t i = 0; i < 50; ++ i)
                    {
                        if(V_ENABLE_BIN_1D_P2[channel]->at(i) == true)
                        {
                            // do nothing
                        }
                        else
                        {
                            continue;
                        }

                        double content = V_PHYS_1D_P2_MATHMORE[channel]->operator[](j_counter).operator[](i_counter);

                        #if DRAWVPHYSMATRIX
                        if(DRAW_V_PHYS_MATRIX == true)
                        {
                            if(channel == 1)
                            {
                                draw_V_PHYSINV_P2->SetBinContent(i + 1, j + 1, content);
                                //std::cout << content << "\t";
                            }
                        }
                        #endif
                        ++ i_counter;
                    }
                        #if DRAWVPHYSMATRIX
                        if(DRAW_V_PHYS_MATRIX == true)
                        {
                            //std::cout << std::endl;
                        }
                        #endif

                    ++ j_counter;
                }


                #if DRAWVPHYSMATRIX
                if(DRAW_V_PHYS_MATRIX == true)
                {
                    if(channel == 1)
                    {
                        std::cout << "ALLOC " << __func__ << std::endl;
                        TCanvas *canvas_V_PHYSINV_P1 = new TCanvas("canvas_V_PHYSINV_P1", "canvas_V_PHYSINV_P1");
                        draw_V_PHYSINV_P1->Draw("colz");
                        canvas_V_PHYSINV_P1->SaveAs("draw_V_PHYSINV_P1.png");

                        TCanvas *canvas_V_PHYSINV_P2 = new TCanvas("canvas_V_PHYSINV_P2", "canvas_V_PHYSINV_P2");
                        draw_V_PHYSINV_P2->Draw("colz");
                        canvas_V_PHYSINV_P2->SaveAs("draw_V_PHYSINV_P2.png");

                        std::cout << "printed all canvas" << std::endl;
                        //std::cin.get();
                        
                        DRAW_V_PHYS_MATRIX = false;
                    }
                }
                #endif
                
                // disable this to recalculate the V MATRIX each loop
                // if the number of singular elements changes it may crash
                // not sure at the moment
                //recalculate_V_PHYS_xD_Px_MATHMORE = false;
            }
        }
    }
}



void
MinimizeFCNAxialVector::calculate_chi2_P1(double &chi2_P1, int &nch_P1) const
{
    chi2_P1 = 0.0;
    nch_P1 = 0;

    for(int channel = 0; channel < number1DHists; ++ channel)
    {
        if(channel_enable_1D[channel] == true)
        {
            //std::cout << "CALCULATE CHI2 FOR PHASE 1" << std::endl;

            // calculate (M - D) x V_PHYS x Transpose(M - D)
            //std::cout << "starting matrix calculations" << std::endl;
            double chi2_1D_P1 = 0.0;

            int l_counter = 0;
            //for(Int_t l{1}; l <= V_PHYS_STAT_1D_P1[channel]->GetNbinsY(); ++ l)
            for(Int_t l{0}; l < M_1D_P1_data[channel]->size(); ++ l)
            {
                if(V_ENABLE_BIN_1D_P1[channel]->at(l) == true)
                {
                    // do nothing
                }
                else
                {
                    continue;
                }

                int m_counter = 0;
                //for(Int_t m{1}; m <= V_PHYS_STAT_1D_P1[channel]->GetNbinsX(); ++ m)
                for(Int_t m{0}; m < M_1D_P1_data[channel]->size(); ++ m)
                // transpose, should be applied to (M-D)->GetNbinsY() but this object does not exist
                {
                    if(V_ENABLE_BIN_1D_P1[channel]->at(m) == true)
                    {
                        // do nothing
                    }
                    else
                    {
                        continue;
                    }

                    //double D_content_1 = D_1D_P1[channel]->GetBinContent(l, 1);
                    //double M_content_1 = M_1D_P1[channel]->GetBinContent(l, 1);
                    #if VECTOR_RANGE_CHECK
                    //double D_content_1 = D_1D_P1_data[channel]->at(l);
                    //double M_content_1 = M_1D_P1_data[channel]->at(l);
                    double D_minus_M_content_1 = D_minus_M_1D_P1_data[channel]->at(l);
                    #else
                    //double D_content_1 = D_1D_P1_data[channel]->operator[](l);
                    //double M_content_1 = M_1D_P1_data[channel]->operator[](l);
                    double D_minus_M_content_1 = D_minus_M_1D_P1_data[channel]->operator[](l);
                    #endif
                    //double delta_1 = D_content_1 - M_content_1;
                    double delta_1 = D_minus_M_content_1;
                    //double V_CHEN_content = V_CHEN->GetBinContent(k, l);
                    //double V_PHYS_STAT_content = V_PHYS_STAT_1D_P1[channel]->GetBinContent(l, m);

                    double V_PHYS_content = 0.0;
                    //std::cout << "l=" << l << " m=" << m << " l_counter=" << l_counter << " m_counter=" << m_counter << std::endl;
                    V_PHYS_content = V_PHYS_1D_P1_MATHMORE[channel]->operator[](m_counter).operator[](l_counter);
                    
                    //double D_content_2 = D_1D_P1[channel]->GetBinContent(m, 1);
                    //double M_content_2 = M_1D_P1[channel]->GetBinContent(m, 1); 
                    #if VECTOR_RANGE_CHECK
                    //double D_content_2 = D_1D_P1_data[channel]->at(m);
                    //double M_content_2 = M_1D_P1_data[channel]->at(m);
                    double D_minus_M_content_2 = D_minus_M_1D_P1_data[channel]->at(m);
                    #else
                    //double D_content_2 = D_1D_P1_data[channel]->operator[](m);
                    //double M_content_2 = M_1D_P1_data[channel]->operator[](m);
                    double D_minus_M_content_2 = D_minus_M_1D_P1_data[channel]->operator[](m);
                    #endif
                    //double delta_2 = D_content_2 - M_content_2;
                    double delta_2 = D_minus_M_content_2;
                    //double next = delta_1 * V_CHEN_content * V_PHYS_STAT_content * delta_2;
                    double next = delta_1 * V_PHYS_content * delta_2;

                    //if(next != 0.0)
                    //{
                    //    std::cout << "m=" << m << " l=" << l << " content=" << V_PHYS_content << " " << D_minus_M_content_1 << " " << D_minus_M_content_2 << " -> " << next << std::endl;
                    //}

                    if(std::isnan(next))
                    {
                        //std::cout << "NAN: next=" << next << " k=" << k << " l=" << l << " m=" << m << std::endl;
                        std::cout << "NAN: next=" << next << " l=" << l << " m=" << m << std::endl;
                        //std::cout << V_CHEN_content << std::endl;
                        std::cout << V_PHYS_content << std::endl;
                        std::cin.get();
                    }
                    else if(std::isinf(next))
                    {
                        //std::cout << "INF: next=" << next << " k=" << k << " l=" << l << " m=" << m << std::endl;
                        std::cout << "INF: next=" << next << " l=" << l << " m=" << m << std::endl;
                        //std::cout << V_CHEN_content << std::endl;
                        std::cout << V_PHYS_content << std::endl;
                        std::cin.get();
                    }
                    chi2_1D_P1 += next;
                    
                    ++ m_counter;
                }
                //std::cout << "after adding strip: " << chi2_1D_P1 << std::endl;

                ++ l_counter;
                ++ nch_P1;
            }

            chi2_P1 += chi2_1D_P1;
            //std::cout << "chi2_1D_P1=" << chi2_1D_P1 << std::endl;
        }

    }

}



void
MinimizeFCNAxialVector::calculate_chi2_P2(double &chi2_P2, int &nch_P2) const
{
    chi2_P2 = 0.0;
    nch_P2 = 0;

    for(int channel = 0; channel < number1DHists; ++ channel)
    {
        if(channel_enable_1D[channel] == true)
        {
            //std::cout << "CALCULATE CHI2 FOR PHASE 2" << std::endl;

            double chi2_1D_P2 = 0.0;
            
            int l_counter = 0;
            //for(Int_t l{1}; l <= V_PHYS_STAT_1D_P2[channel]->GetNbinsY(); ++ l)
            for(Int_t l{0}; l < 50; ++ l)
            {
                if(V_ENABLE_BIN_1D_P2[channel]->at(l) == true)
                {
                    // do nothing
                }
                else
                {
                    continue;
                }

                int m_counter = 0;
                //for(Int_t m{1}; m <= V_PHYS_STAT_1D_P2[channel]->GetNbinsX(); ++ m)
                for(Int_t m{0}; m < 50; ++ m)
                // transpose, should be applied to (M-D)->GetNbinsY() but this object does not exist
                {
                    if(V_ENABLE_BIN_1D_P2[channel]->at(m) == true)
                    {
                        // do nothing
                    }
                    else
                    {
                        continue;
                    }

                    //double D_content_1 = D_1D_P2[channel]->GetBinContent(l, 1);
                    //double M_content_1 = M_1D_P2[channel]->GetBinContent(l, 1);
                    #if VECTOR_RANGE_CHECK
                    //double D_content_1 = D_1D_P2_data[channel]->at(l);
                    //double M_content_1 = M_1D_P2_data[channel]->at(l);
                    double D_minus_M_content_1 = D_minus_M_1D_P2_data[channel]->at(l);
                    #else
                    //double D_content_1 = D_1D_P2_data[channel]->operator[](l);
                    //double M_content_1 = M_1D_P2_data[channel]->operator[](l);
                    double D_minus_M_content_1 = D_minus_M_1D_P2_data[channel]->operator[](l);
                    #endif
                    //double delta_1 = D_content_1 - M_content_1;
                    double delta_1 = D_minus_M_content_1;
                    //double V_CHEN_content = V_CHEN->GetBinContent(k, l);
                    //double V_PHYS_STAT_content = V_PHYS_STAT_1D_P2[channel]->GetBinContent(l, m);
                    double V_PHYS_content = 0.0;

                    //std::cout << "l=" << l << " m=" << m << " l_counter=" << l_counter << " m_counter=" << m_counter << std::endl;
                    V_PHYS_content = V_PHYS_1D_P2_MATHMORE[channel]->operator[](m_counter).operator[](l_counter);


                    //double D_content_2 = D_1D_P2[channel]->GetBinContent(m, 1);
                    //double M_content_2 = M_1D_P2[channel]->GetBinContent(m, 1);
                    #if VECTOR_RANGE_CHECK
                    //double D_content_2 = D_1D_P2_data[channel]->at(m);
                    //double M_content_2 = M_1D_P2_data[channel]->at(m);
                    double D_minus_M_content_2 = D_minus_M_1D_P2_data[channel]->at(m);
                    #else
                    //double D_content_2 = D_1D_P2_data[channel]->operator[](m);
                    //double M_content_2 = M_1D_P2_data[channel]->operator[](m);
                    double D_minus_M_content_2 = D_minus_M_1D_P2_data[channel]->operator[](m);
                    #endif
                    //double delta_2 = D_content_2 - M_content_2;
                    double delta_2 = D_minus_M_content_2;
                    //double next = delta_1 * V_CHEN_content * V_PHYS_STAT_content * delta_2;
                    double next = delta_1 * V_PHYS_content * delta_2;

                    //if(next != 0.0)
                    //{
                    //    std::cout << "m=" << m << " l=" << l << " content=" << V_PHYS_content << " " << D_minus_M_content_1 << " " << D_minus_M_content_2 << " -> " << next << std::endl;
                    //}

                    if(std::isnan(next))
                    {
                        //std::cout << "NAN: next=" << next << " k=" << k << " l=" << l << " m=" << m << std::endl;
                        std::cout << "NAN: next=" << next << " l=" << l << " m=" << m << std::endl;
                        //std::cout << V_CHEN_content << std::endl;
                        std::cout << V_PHYS_content << std::endl;
                        std::cin.get();
                    }
                    else if(std::isinf(next))
                    {
                        //std::cout << "INF: next=" << next << " k=" << k << " l=" << l << " m=" << m << std::endl;
                        std::cout << "INF: next=" << next << " l=" << l << " m=" << m << std::endl;
                        //std::cout << V_CHEN_content << std::endl;
                        std::cout << V_PHYS_content << std::endl;
                        std::cin.get();
                    }
                    chi2_1D_P2 += next;
                    
                    ++ m_counter;
                }
                //std::cout << "after adding strip: " << chi2_1D_P2 << std::endl;

                ++ l_counter;
                ++ nch_P2;
            }

            chi2_P2 += chi2_1D_P2;
            //std::cout << "chi2_1D_P2=" << chi2_1D_P2 << std::endl;
        }
    }

}


void
MinimizeFCNAxialVector::calculate_penalty_term(double& penalty_term_ret, const std::vector<double> &param) const
{
    // penalty terms section
    double penalty_term = 0.0;

    // loop over all the parameters
    std::map<int, file_parameter>::iterator it{g_pg.file_params.begin()};
    for(; it != g_pg.file_params.end(); ++ it)
    {
        int paramNumberInt = -1;

        int paramNumber = it->second.paramNumber;
        bool paramEnabled = it->second.paramEnabled;
        bool paramEnabledP1 = it->second.paramEnabledP1;
        bool paramEnabledP2 = it->second.paramEnabledP2;
        double paramInitValue = it->second.paramInitValue;
        double paramInitError = it->second.paramInitError;
        double paramConstraintValue = it->second.paramConstraintValue;
        double paramConstraintError = it->second.paramConstraintError;
        int paramConstraintMode = it->second.paramConstraintMode;

        // stop it crashing
        // (below)
        if(paramEnabled == false)
        {
            continue;
        }

        // ignore parameter if it is of background type and FORCE_BKG_HARD
        // is set
        if(FORCE_BKG_HARD == true)
        {
            if(paramNumber > 1)
            {
                continue;
            }
        }

        paramNumberInt = g_pg.ExtToIntParamNumberMap.at(paramNumber);

        if(debuglevel >= 5)
        {
            std::cout << "paramNumber=" << paramNumber << " -> " << paramNumberInt << std::endl;
        }

        if(paramEnabled == true)
        {
            if((paramEnabledP1 == true) || (paramEnabledP2 == true))
            {
                // do nothing
            }
            else
            {
                // not enabled
                continue;
            }
        }
        else
        {
            // not enabled
            continue;
        }

        if(paramConstraintMode == MODE_CONSTRAINT_SOFT)
        {
            if(debuglevel >= 6)
            {
                std::cout << "soft" << std::endl;
            }
            // do nothing
        }
        else if(paramConstraintMode == MODE_CONSTRAINT_FREE)
        {
            //-- ndf;
            if(debuglevel >= 6)
            {
                std::cout << "free" << std::endl;
            }
            continue;
        }
        else if(paramConstraintMode == MODE_CONSTRAINT_HARD)
        {
            if(debuglevel >= 6)
            {
                std::cout << "hard" << std::endl;
            }
            continue;
        }
        else
        {
            std::cout << "ERROR: Invalid value for paramNumber=" << paramNumber << ", paramConstraintMode=" << paramConstraintMode << std::endl;
        }

        // this parameter is from minuit internal and thus is in minuit
        // internal units (not Bq)
        // have to convert to Bq units
        double param_value = param.at(paramNumberInt); // TODO this might break if we try to fit P1 seperatly and P2 seperatly
        double activity_value_Bq = paramInitValue;
    
        // convert to Bq
        // multiply by the initial value
        // activity_value_Bq should really be called param_initial_value
        //double activity_value_Bq = 0.0;
        //double tmp_err;
        //get_paramInitValueError(thePhase, i, activity_value_Bq, tmp_err);

        double value = param_value * activity_value_Bq;
        //double penalty = std::pow((param_value - constraint) / error, 2.0);
        double penalty = 0.0;

        // TODO
        /*
        if(EMODE == 1)
        {
            // data
        }
        else if(EMODE == 2)
        {
            // MC
        }   
        else if(EMODE == 3)
        {
            // quadrature
        }
        */
        
        penalty = std::pow((value - paramConstraintValue) / paramConstraintError, 2.0);
        if(debuglevel >= 5)
        {
            //std::cout << "j=" << j << std::endl;
            std::cout << "paramNumber=" << paramNumber
                      << " value=" << value
                      << " param_value=" << param_value
                      << " paramConstraintValue=" << paramConstraintValue
                      << " paramConstraintError=" << paramConstraintError
                      << " penalty=" << penalty
                      << std::endl;
        }


        // TODO: is this the correct error term?
        // error on constraint rather than error on current fit value?

        penalty_term += penalty;
    }

    penalty_term_ret = penalty_term;
}

#endif //MINIMIZEFCNAXIALVECTOR_CHECK_ALLOC_H
